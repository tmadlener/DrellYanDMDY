
# DrellYanDMDY package

This package runs the analysis of DY->ee and computes dSigma/dM and d2Sigma/dMdY 
cross sections. 

The package is organized by directories that correspond to steps of the analysis flow.
Additionally, there are special directories as described below.

## Package directory structure

The list of directories with brief comments follows:

**FullChain** contains the master script that runs the full analysis

**root_files** is where output fiels of various sorts are written in all
intermediate stages.

**config_files** contains configuration files that define samples, cross sections, 
and other conditions.

**Include** contains most of the include files and some class definitions for code
used from multiple steps.

**Selection, YieldsAndBackgrounds, Unfolding, Acceptance, Efficiency, EventScaleFactors, CrossSection**
directories correspond to scripts responsible for specific analysis steps.

**Skimming** contains scripts needed to reduce first-level ntuples to smaller ones.

**LumiAndPileup** contains scripts related to preparation of pileup weights.

**Theory** has macros that prepare and plot theoretical cross section (the actual
computation does not take place here, but comes from elsewhere).

**EnergyScaleFactors, DataDrivenBackgrounds, Fsr** are other directories for use by experts.

**Logs, logs-1D, logs-2D** are directories that contain logs created when the package runs.

## Most important files

**FullChain/FullChain.sh** is the master script. Inside it one can set which steps of
the analysis flow should be run. One can run this script multiple times switching
on and off difference pieces. The full time required to run the analysis could
be as long as 24 hours (without systematics). Different parts of the analysis sometimes
depend on each other, so a proper sequence is necessary for some steps, as specified below:
 - do_selection, do_prepareYields, do_subtractBackground have to be done in this order
 - do_unfoldingFsr, do_acceptance, do_efficiency, do_efficiencyScaleFactors, 
        do_plotFSRCorrections, do_plotFSRCorrectionsSansAcc, do_setupTheory can be done in any
        order, and in any order with the selection above
 - do_crossSectionFsr, do_plotXSec have to be done after all the above is done

**config_files/data8TeV.conf, config_files/summer12mc.input, config_files/xsecCalc8TeV.conf** set where
to find samples and some parameters like cross sections or directory names. Note that these three
configuration scripts are listed in the very beginning of the FullChain.sh.

**root_files/selected_events/DY_j22_19712pb/** (the last subdirectory is configurable, set in 
the config scripts via parameter called dirTag) is where one finds simple ntuples with only the 
events that pass selection.

## How to run the master script

After the master script is ready and all do_xxx flags are set, one can run it interactively
in a shell simply as
```
./FullChain.sh >& mylog.txt &
```
Then one can use, e.g., "tail -f mylog.txt" to watch the progress. At the end, the script will
print out the summary of whether the requested steps of the analysis succeeded or not.

The script has been tested on lxplus with ROOT 6.06 (e.g. available if one does setenv in CMSSW_8_X_X 
or later).
