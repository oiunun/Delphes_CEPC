source /cvmfs/cepc.ihep.ac.cn/software/cepcenv/setup.sh
cepcenv use 0.1.0-rc9
export PATH=${PATH}:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Framework/Marlin/01-05/bin
export FSClasser_HOME=/workfs2/bes/lig/higgs/FSClasser
export FastJet_HOME=/besfs5/groups/higgs/users/lig/analysis/LCFIplus/fastjet/install
export FastJet_BIN=${FastJet_HOME}/bin
export FastJet_INCLUDE=${FastJet_HOME}/include
export FastJet_LIB=${FastJet_HOME}/lib
unset  MARLIN_DLL
export MARLIN_DLL=${MARLIN_DLL}:$FastJet_LIB/libfastjet.so
export MARLIN_DLL=${MARLIN_DLL}:$FastJet_LIB/libsiscone.so
export MARLIN_DLL=${MARLIN_DLL}:$FastJet_LIB/libsiscone_spherical.so
export MARLIN_DLL=${MARLIN_DLL}:$FastJet_LIB/libfastjetplugins.so
export MARLIN_DLL=${MARLIN_DLL}:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Tracking/MarlinTrkProcessors/01-10/lib/libMarlinTrkProcessors.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Tracking/ForwardTracking/01-07/lib/libForwardTracking.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Tracking/Clupatra/00-10/lib/libClupatra.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Analysis/SubDetectorStudy/TPC/MarlinTPC/00-16/lib/libMarlinTPC.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Pandora/MarlinPandora/00-11/lib/libMarlinPandora.so:
export MARLIN_DLL=${MARLIN_DLL}:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Garlic/2.10.1/lib/libGarlic.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/HighLevelObjectFinding/Jets/LCFIVertex/00-06-01/lib/libLCFIVertex.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Digitization/MarlinReco/01-09/lib/libMarlinReco.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/EventDisplay/CEDViewer/01-07-02/lib/libCEDViewer.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Pandora/PandoraAnalysis/00-05/lib/libPandoraAnalysis.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Analysis/Overlay/00-13/lib/libOverlay.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Arbor/3.4.1/lib/libRanger.so:
export MARLIN_DLL=${MARLIN_DLL}:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/HighLevelObjectFinding/Jets/MarlinFastJet/00-02/lib/libMarlinFastJet.so:
export MARLIN_DLL=${MARLIN_DLL}:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Framework/LCTuple/01-03/lib/libLCTuple.so:
export MARLIN_DLL=${MARLIN_DLL}:$FSClasser_HOME/lib/libFSClasser.so

#export LD_LIBRARY_PATH="$FSClasser_HOME/lib:$LD_LIBRARY_PATH"
#export MARLIN_DLL=$MARLIN_DLL:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Framework/LCIO/02-04-03/lib/liblcio.so

