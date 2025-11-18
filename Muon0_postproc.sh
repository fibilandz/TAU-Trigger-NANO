#!/bin/bash
cd /afs/cern.ch/work/f/fbilandz/CMSSW_14_1_0_pre4/src/PhysicsTools/NanoAODTools/TAU-Trigger-NANO && cmsenv
export X509_USER_PROXY=$1
voms-proxy-info -all

python3 postproc_jec.py \
 --input $2 \
 --isMC 0 \
 --era 2024 \
 --output /eos/home-f/fbilandz/Muon02024C_skims/