#include <time.h>
#include <fstream>
#include "FreeFunctions.h"
#include "GlobalPrecisionParameters.h"

using namespace Eigen;

#define jprime couplingConstants[2]
#define h couplingConstants[3]

int main()
{
    clock_t start = clock();
    
    std::ifstream filein("Input/Input");
    if(!filein)
    {
        std::cerr << "Couldn't open input file." << std::endl;
        exit(EXIT_FAILURE);
    };
    
    // read in parameters that are constant across all trials:
    int nTrials;
    filein >> nTrials;
    stepData data; // this struct will contain most of the important parameters
    data.compBlock = data.beforeCompBlock = NULL;
    for(int trial = 1; trial <= nTrials; trial++)
    {
        clock_t startTrial = clock();
        std::cout << "Trial " << trial << ":" << std::endl;
        std::ofstream fileout("Output/Trial_" + std::to_string(trial));
        fileout << "Trial " << trial << ":\n" << std::endl;
        
        // read in parameters that vary over trials:
        int lSys;                                              // system length
        filein >> lSys;
        std::vector<double> couplingConstants(nCouplingConstants);
        for(int i = 0; i < nCouplingConstants; i++)
            filein >> couplingConstants[i];
        int wavelength;       // wavelength of interstitial spin pattern ansatz
        filein >> wavelength;
        std::vector<double> intSpinPattern(wavelength);
                                     // repeating pattern of interstitial spins
        for(int i = 0; i < wavelength; i++)
            filein >> intSpinPattern[i];
        int nSweeps;                        // number of sweeps to be performed
        filein >> data.mMax >> nSweeps;
        std::vector<double> lancTolerances(nSweeps + 1);
        for(int sweep = 0; sweep <= nSweeps; sweep++)
            filein >> lancTolerances[sweep];
        
        fileout << "System length: " << lSys << "\nCoupling constants:";
        for(double couplingConstant : couplingConstants)
            fileout << " " << couplingConstant;
        fileout << "\nInitial pattern: ";
        for(double sz : intSpinPattern)
            fileout << sz << " ";
        fileout << "\nMaximum bond dimension: "
                << data.mMax << "\nNumber of sweeps: " << nSweeps
                << "\nLanczos tolerances:";
        for(double lancTolerance : lancTolerances)
            fileout << " " << lancTolerance;
        fileout << std::endl << std::endl;
        data.ham.setParams(couplingConstants, lSys);
        int skips = 0,
            runningKeptStates = d * d;
        for(; runningKeptStates <= data.mMax; skips++)
            runningKeptStates *= d;  // find how many edge sites can be skipped
        bool oddSize = lSys % 2;
        int lSFinal,                        // final length of the system block
            lEFinal;                   // final length of the environment block
        if(oddSize)
        {
            lSFinal = (lSys - 1)/2;
            lEFinal = (lSys - 3)/2;
        }
        else
            lSFinal = lEFinal = lSys / 2 - 1;
        bool completeED = false;
        if(skips + 1 >= lSFinal)
        {
            if(skips + 1 == lSFinal && runningKeptStates == data.mMax * d)
            {
                std::cout << "Note: the maximum bond dimension is large enough "
                          << "to perform exact diagonalization." << std::endl;
                completeED = true;
            }
            else
            {
                std::cout << "Error: the maximum bond dimension is larger than "
                          << "required for exact diagonalization." << std::endl;
                continue;
            };
        };
        std::vector<TheBlock> westBlocks(lSys - 2 - skips),
                              eastBlocks(lSys - 2 - skips);
             // initialize system - the last block is only used for odd-size ED
        TheBlock* eastBlocksStart = eastBlocks.data();
        std::vector<double> intSpins;
                         // initial ansatz for interstitial spin magnetizations
        intSpins.reserve(lSys);
        for(int i = 0; i < lSys; i++)
            intSpins.push_back(intSpinPattern[i % wavelength]);
        data.ham.calcEffectiveH(intSpins);
        westBlocks.front() = TheBlock(data.ham, true);
        eastBlocks.front() = TheBlock(data.ham, false);
                                         // initialize the edge one-site blocks
        std::cout << "Performing iDMRG..." << std::endl;
            // note: this iDMRG code assumes parity symmetry of the Hamiltonian
        data.sweepingEast = true;
        data.exactDiag = true;
        data.compBlock = eastBlocksStart;
        data.infiniteStage = true;
        data.lancTolerance = lancTolerances.front();
        rmMatrixX_t psiGround;                    // seed for Lanczos algorithm
        double cumulativeTruncationError = 0.;
        for(int site = 0; site < skips; site++, data.compBlock++) // initial ED
            westBlocks[site + 1]
                = westBlocks[site].nextBlock(data, psiGround,
                                             cumulativeTruncationError,
                                             eastBlocksStart + site + 1);
        data.exactDiag = completeED;
        for(int site = skips, end = lEFinal - 1; site < end;
            site++, data.compBlock++)                                  // iDMRG
        {
            psiGround = randomSeed(westBlocks[site], eastBlocks[site]);
            westBlocks[site + 1]
                = westBlocks[site].nextBlock(data, psiGround,
                                             cumulativeTruncationError,
                                             eastBlocksStart + site + 1);
        };
        if(oddSize)
        {
            data.compBlock = eastBlocksStart + (lSFinal - 2);
            psiGround = randomSeed(westBlocks[lSFinal - 2],
                                   eastBlocks[lSFinal - 2]);
            westBlocks[lSFinal - 1]
                = westBlocks[lSFinal - 2].nextBlock(data, psiGround,
                                                    cumulativeTruncationError);
        };
        fileout << "iDMRG average truncation error: "
                << cumulativeTruncationError / (lSys - 2) << std::endl
                << std::endl;         // handles both even and odd system sizes
        if(completeED || nSweeps == 0)
            psiGround = randomSeed(westBlocks[lSFinal - 1],
                                   eastBlocks[lEFinal - 1]);
        else
        {
            std::cout << "Performing fDMRG..." << std::endl;
            data.infiniteStage = false;
            int endSweep = lSys - 4 - skips;              // last site of sweep
            psiGround = randomSeed(westBlocks[lSFinal - 1],
                                   eastBlocks[lEFinal - 1]);
            double chainEnergy;
               // contribution to GS energy captured directly by DMRG algorithm
            for(int sweep = 1; sweep <= nSweeps; sweep++)
                                                    // perform the fDMRG sweeps
            {
                fileout << "Sweep " << sweep
                        << ":\nInterstitial spin polarizations:" << std::endl;
                for(double d : intSpins)
                    fileout << d << " ";
                fileout << std::endl << std::endl;
                data.compBlock = eastBlocksStart + (lEFinal - 1);
                data.lancTolerance = lancTolerances[sweep];
                data.beforeCompBlock = data.compBlock - 1;
                cumulativeTruncationError = 0.;
                for(int site = lSFinal - 1; site < endSweep;
                    site++, data.compBlock--, data.beforeCompBlock--)
                    westBlocks[site + 1]
                        = westBlocks[site].nextBlock(data, psiGround,
                                                     cumulativeTruncationError);
                reflectPredictedPsi(psiGround, westBlocks[endSweep],
                                    eastBlocks[skips]);
                               // reflect the system to reverse sweep direction
                data.sweepingEast = false;
                data.compBlock = &westBlocks[endSweep];
                data.beforeCompBlock = data.compBlock - 1;
                for(int site = skips; site < endSweep;
                    site++, data.compBlock--, data.beforeCompBlock--)
                    eastBlocks[site + 1]
                        = eastBlocks[site].nextBlock(data, psiGround,
                                                     cumulativeTruncationError);
                reflectPredictedPsi(psiGround, eastBlocks[endSweep],
                                    westBlocks[skips]);
                data.sweepingEast = true;
                data.compBlock = eastBlocksStart + endSweep;
                data.beforeCompBlock = data.compBlock - 1;
                for(int site = skips, end = lSFinal - 1; site < end;
                    site++, data.compBlock--, data.beforeCompBlock--)
                    westBlocks[site + 1]
                        = westBlocks[site].nextBlock(data, psiGround,
                                                     cumulativeTruncationError);
                std::cout << "Sweep " << sweep << " complete." << std::endl;
                fileout << "Average truncation error: "
                        << cumulativeTruncationError / (2 * lSys - 4)
                        << std::endl;
                data.infiniteStage = false;
                FinalSuperblock hSuperFinal
                    = westBlocks[lSFinal - 1].createHSuperFinal(data, psiGround,
                                                                skips);
                                               // calculate ground-state energy
                chainEnergy = hSuperFinal.gsEnergy;
                fileout << "Chain ground state energy density: "
                        << chainEnergy / lSys << std::endl << std::endl;
                std::cout << "Calculating observables..." << std::endl;
                obsMatrixD_t sz;
                sz << .5, 0.,
                      0., -.5;
                VectorXd oneSiteSzs
                    = oneSiteExpValues(sz, lSys, hSuperFinal, westBlocks,
                                       eastBlocks);
                fileout << "Expectation value of S_z at each chain site:\n"
                        << oneSiteSzs << std::endl << std::endl;
                std::vector<double> nnCorrelations,
                                    nnnCorrelations;
                nnCorrelations.reserve(lSys - 1);
                nnnCorrelations.reserve(lSys - 2);
                for(int i = 0, end = lSys - 2; i < end; i++)
                {
                    nnCorrelations.push_back(expiDotj(i, i + 1, hSuperFinal,
                                                      westBlocks, eastBlocks));
                    nnnCorrelations.push_back(expiDotj(i, i + 2, hSuperFinal,
                                                       westBlocks, eastBlocks));
                };
                nnCorrelations.push_back(expiDotj(lSys - 2, lSys - 1, hSuperFinal,
                                                  westBlocks, eastBlocks));
                fileout << "Expectation values of S_i dot S_{i + 1}:"
                        << std::endl;
                for(double d : nnCorrelations)
                    fileout << d << std::endl;
                fileout << std::endl
                        << "Expectation values of S_i dot S_{i + 2}:\n";
                for(double d : nnnCorrelations)
                    fileout << d << std::endl;
                fileout << std::endl;
                // Determine the polarizations of the next iterations'
                // interstitial spins:
                intSpins.clear();
                RowVector3d hTotal;
                              // induced + applied fields on interstitial spins
                for(int i = 0; i < lSys - 1; i++)
                    intSpins
                        .push_back(((  (-jprime * (  oneSiteSzs(i)
                                                   + oneSiteSzs(i + 1)) + h / 2)
                                     > 0) - .5) * (d - 1));
                intSpins.push_back((((  -jprime * oneSiteSzs(lSys - 1)
                                      + h / 2) > 0) - .5) * (d - 1));
                      // induced + applied field on rightmost interstitial spin
                data.ham.calcEffectiveH(intSpins);
            };
            fileout << "Final interstitial spin polarizations:" << std::endl;
            for(double d : intSpins)
                fileout << d << " ";
            fileout << std::endl << std::endl
                    << "Chain ground state energy density: "
                    << chainEnergy / lSys << std::endl << std::endl;
        };
        clock_t stopTrial = clock();
        fileout << "Elapsed time: "
                << float(stopTrial - startTrial) / CLOCKS_PER_SEC << " s"
                << std::endl;
        fileout.close();
    };
    filein.close();
    
    clock_t stop = clock();
    std::cout << "Done. Elapsed time: " << float(stop - start) / CLOCKS_PER_SEC
              << " s" << std::endl;
    
    return 0;
};
