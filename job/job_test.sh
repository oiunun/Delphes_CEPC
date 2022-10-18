#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /cefs/higgs/gaoxu/delphes
./DelphesSTDHEP cards/delphes_card_CEPC_4th.tcl /cefs/higgs/gaoxu/delphes/rootfile/e2e2h_X_1_99_1.root /cefs/data/stdhep/CEPC240/higgs/update_from_LiangHao_1M/data/E240.Pe2e2h_X.e0.p0.whizard195/e2e2h_X.e0.p0.000*.stdhep