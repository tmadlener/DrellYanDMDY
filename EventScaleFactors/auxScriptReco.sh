#!/bin/bash

# This script runs TnP for RECO efficiencies
# stand-alone, not part of full chain. Only
# single electron efficiencies are computed here.

echo "Start MC"
root -b -q eff_Reco.C+\(\"../config_files/sf_8TeV_mc_et6_eta5.conf\",\"RECO\",\"Full2012_hltEffOld\",1,0\) >& log-reco-pu5minus-mc-1.txt
root -b -q calcEff.C+\(\"../config_files/sf_8TeV_mc_et6_eta5.conf\",\"RECO\",\"Full2012_hltEffOld\",1,0\) >& log-reco-pu5minus-mc-2.txt

echo "Start data"
root -b -q eff_Reco.C+\(\"../config_files/sf_8TeV_data_et6_eta5.conf\",\"RECO\",\"Full2012_hltEffOld\",1,0\) >& log-reco-pu5minus-data-1.txt
root -b -q calcEff.C+\(\"../config_files/sf_8TeV_data_et6_eta5.conf\",\"RECO\",\"Full2012_hltEffOld\",1,0\) >& log-reco-pu5minus-data-2.txt

echo "Done"
