step1: initialize your grid certificate by grid-proxy-init -debug -verify

step2: run following command to create the scripts to submit condor jobs
#python  produceSample.py -i /eos/uscms/store/group/lpcgem/SingleMuon_Run2017G_v1_RECO/ -o  /eos/uscms/store/user/mkhurana/GEMCSCBending_2017G/
<<<<<<< Updated upstream
#python  produceSample.py -i /eos/uscms/store/user/mkhurana/2018C_data_files/ -o  /eos/uscms/store/user/mkhurana/GEMCSCBending_2018C/ -c  runSliceTestAnalysis_condor_new_2018C.py
python  produceSample.py -i /eos/uscms/store/group/lpcgem/SingleMuon_Run2018C_v1_RECO/ -o  /eos/uscms/store/user/mkhurana/GEMCSCBending_2017G/ -c  runSliceTestAnalysis_condor_new_2018C.py

=======
python  produceSample.py -i /eos/uscms/store/user/mkhurana/2018C_data_files/ -o  /eos/uscms/store/user/mkhurana/GEMCSCBending_2018C/ -c  runSliceTestAnalysis_condor_new_2018C.py
python  produceSample.py -i /eos/uscms/store/group/lpcgem/SingleMuon_Run2017G_v1_RECO/ -o  /eos/uscms/store/user/tahuang/GEMCSCBending_2017G/
python  produceSample.py -i /eos/uscms/store/group/lpcgem/SingleMuon_Run2018C_v1_RECO/ -o  /eos/uscms/store/user/tahuang/GEMCSCBending_2018C/
>>>>>>> Stashed changes

step3:  run submitalljobs.sh




to test the configuration

#change username accordingly
<<<<<<< Updated upstream

cmsRun runSliceTestAnalysis_condor_new_2018C.py inputFiles=file:/eos/uscms/store/group/lpcgem/SingleMuon_Run2018C_v1_RECO/step3_313.root outputFile=/uscms_data/d3/mkhurana/CMSSW_10_1_5/src/GEMCSCBendingAnalyzer/MuonAnalyser/condor/test.root
#cmsRun runSliceTestAnalysis_condor.py inputFiles=file:/eos/uscms/store/group/lpcgem/SingleMuon_Run2017G_v1_RECO/step3_313.root outputFile=/eos/uscms/store/user/tahuang/GEMCSCBending_2017G/out_ana_0.root
=======
cmsRun runSliceTestAnalysis_condor.py inputFiles=file:/eos/uscms/store/group/lpcgem/SingleMuon_Run2017G_v1_RECO/step3_313.root outputFile=/eos/uscms/store/user/tahuang/GEMCSCBending_2017G/out_ana_0.root
cmsRun runSliceTestAnalysis_condor.py inputFiles=file:/eos/uscms/store/group/lpcgem/SingleMuon_Run2018C_v1_RECO/step3_000.root outputFile=/eos/uscms/store/user/tahuang/GEMCSCBending_2018C/out_ana_0.root
>>>>>>> Stashed changes
