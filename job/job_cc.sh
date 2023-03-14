#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh
cd /cefs/higgs/gaoxu/delphes
./DelphesHepMC2 cards/delphes_card_CEPC_4th.tcl rootfile/cc.root ../pythia8307/code/cc/cc_*.hepmc