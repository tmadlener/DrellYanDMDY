
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

## Output of different steps of the analysis flow

The overall formula for a typical (and our) cross section calculation
```
    sigma = Unfold_in_M_and_Y ( N_selected_data - N_background) / ( Acceptance * Efficiency * Eff_scale_Factors * Luminosity)
```
Note on binning: the 2D measurement with respect to (M, Y) is set up to have equal size bins
    in rapidity, but arbitrary bin sizes in mass. Additionally, the rapidity bins have to be equal
    size, for each mass range, but one can have different number of rapidity bins for different
    mass ranges. For example, for 20<M<30 one can have 24 rapidity bins, and for 30<M<40
    one can have 12 mass bins. 

For the steps "prepare yields" and "subtract backgrounds" the results a matrices of yields
stored in ROOT files with similar structure. The "yields" part is explained below in detail:

```
    root_files/yields/DY_j22_19712pb/yields2D.root
```
      This file contains:
        1) Binning info for consistency checking:
```
                 TVectorT<double> massBinLimits
                 TVectorT<double> rapidityCounts
```
          As the name suggests, the limits of mass bins, and then number of rapidity
          bins for each mass bin.
        2) The yields in TMatrix form for the observed in data event counts yields_data
          as well as yields observed in all MC samples: yields_zee is signal MC, and the
          rest like yields_wjets, yields_qcd, etc, are backgrounds. Note that the MC yields
          are already contain proper pile-up, luminosity and cross section weighting.
          (however there is no scale factor corrections at this stage).
        3) The TMatrix objects for each case above with sum of weights squared 
           (the error on the yield in a given bin is sqrt(sumw^2), so it is the squares of errrors).
           The matrix names are yieldsSumw2_xxx where xxx is the sample type.
       One can open the file and quickly view the content by doing executing .ls first to see
        the list of available objects, and then printing them, like yields_data->Print() interactively.

```
    root_files/yields/DY_j22_19712pb/yields_plots2D.root
```
       This file contains a number of plots. For a 2D cross section, the spectra of rapidity
      are drawn separately for each mass range, data with MC overlaid, in two normalizations:
      in each mass range total MC is normalized to have the same area for that mass range,
      and also in each mass range the same normalization is applied so that the Z peak area
      is the same for data and MC.
          Additionally, there is a "flattened" plot in which both mass and rapidity is represented
      along one axis. There are a few more useful plots there.

```
    YieldsAndBackgrounds/plots2D/
```
        The directory contains plots in png, pdf, root formats, of the type described above.

```
    YieldsAndBackgrounds/tables2D/
```
        The directory contains tables of yields in text/latex format suitable for pasting into
      an Analysis Note.
  
The output of the "subtract backgrounds" is similar to the above, a bit smaller in scope:

```
     root_files/yields/DY_j22_19712pb/yields_bg-subtracted2D.root
```
          This file contains the binning info for consistency checking again massBinLimits
          and rapidityCounts, as well as the background subtracted yields. There are three
          matrices:
```
            TMatrixT<double> YieldsSignal
            TMatrixT<double> YieldsSignalErr
            TMatrixT<double> YieldsSignalSystErr
```
          The first error is statistical, the second is systematic as the name suggests.
          These are errors, nor error squares like it was for the "yields" part.

```
     root_files/yields/DY_j22_19712pb/yields_bg-subtracted2D-plots.root
```
        This file contains a few simple plots, not as useful as the "yields" plots

```
    YieldsAndBackgrounds/plots2D/
    YieldsAndBackgrounds/tables2D/
```
          These are same as above, for the "yields" step of the workflow. The tables and
          plots from both steps are found here.

For the steps "acceptance", "efficiency", "event scale factors", the result is constants
that will be applied in the final cross section calculation. The constants are stored
in ROOT files for later use. There are also plots to look at. The structure is the same
for all of these steps. Let me explain with the "acceptance" step:

```
   root_files/constants/DY_j22_19712pb/acceptance_constants2D.root
```
       This file contains:
          1) Binning info for consistency checking:
```
                 TVectorT<double> massBinLimits
                 TVectorT<double> rapidityCounts
```
              These arrays contain, as the name suggests, the limits of mass bins, and then number of rapidity
              bins for each mass bin.
          2) The acceptance corrections themselves:
```
                TMatrixT<double> acceptanceMatrix
                TMatrixT<double> acceptanceErrMatrix
```
              These are matrices with the acceptance constants and errors on them (errors from MC statistics)
               for all mass/rapidity bins.
           One can easily see the content by simply opening this root file and
           running commands like acceptanceMatrix->Print() interactively.

```
     root_files/constants/DY_j22_19712pb/acceptance_plots2D.root
```
          This file contains a canvas with the acceptance correction plot in 2D. One can open it
          and execute acceptance->Draw() to see the plot.

```
     Acceptance/plots2D/
```
            Here one can find acceptance and related plots in png, pdf, root format. 

```
     Acceptance/tables2D/
```
            In this directory, one can find tables of acceptance constants in latex format, suitable
            for an Analysis Note.

     The steps "efficiency" and "event scale factors" result in a similar outcome, the relevant
     locations/names are:

Efficiency: (no binning arrays saved here, only the corrections themselves)
```
 root_files/constants/DY_j22_19712pb/event_efficiency_constants2D.root
 root_files/constants/DY_j22_19712pb/event_efficiency_plots2D.root 
 Efficiency/plots2D/
 Efficiency/tables2D/
```

Event scale factors:
```
  root_files/constants/DY_j22_19712pb/scale_factors_2D_Full2012_hltEffOld_PU.root
  root_files/constants/DY_j22_19712pb/scale_factors_2D_Full2012_hltEffOld_PU-plots.root
  EventScaleFactors/plots2D_PU/
```
Above in the first file, the ROOT file with constants, there are not only mass and rapidity binnings
for consistency checking, but also single electron pt and eta binning that tells us what was
used when running tag-and-probe to find single electron efficiencies for scale factors.
    