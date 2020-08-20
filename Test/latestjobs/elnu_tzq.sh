#!/bin/bash
export X509_USER_PROXY=x509up_u527
source ~/.bashrc
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh
./build/dptest calchisto elnu tzq
