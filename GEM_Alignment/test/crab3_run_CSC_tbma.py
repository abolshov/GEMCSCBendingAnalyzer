#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = 'CSC_tbma_CRAFT2022_050522_ExpressCosmics_Initial'
config.General.workArea = 'RdPhiAna_CRAFT2022_ExpressCosmics_TerukiRuns'#working dir
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_CSC_tbma.py'
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles

misalign = False  #Make sure to change the run_GE11ana.py too!!!
if misalign:
  #config.JobType.inputFiles = ['./UL2018IOV3_RE_CSC_03.db']
  #config.JobType.inputFiles = ['./ME11_UL2018IOV_iter1.csv.db']
  #config.JobType.inputFiles = ['./CRAFT2022_ME11Only_Initial.db']
  config.JobType.inputFiles = ['./CRAFT2022_ME11Only_Iter2.db']

#section Data
#config.Data.inputDataset = '/singleMuonGun_MuAl_pT-30to200_1102_phase1_2021_realistic/hyunyong-crab_singleMuonGun_pT-30to200_1102_phase1_2021_realistic_RAW2DIGI_FullRECOv4-1b4eba2dcd577d6bb642bb3e45609e5f/USER'
#config.Data.inputDataset = '/Cosmics/Commissioning2021-CosmicSP-PromptReco-v1/RAW-RECO'
#config.Data.inputDataset = '/Cosmics/Commissioning2021-CosmicTP-PromptReco-v1/RAW-RECO'
#config.Data.inputDataset = '/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/hyunyong-crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3-ce33467258ba6d6e1b97c4b94c6a3a02/USER'
#config.Data.inputDataset = '/Cosmics/Commissioning2021-PromptReco-v1/AOD' #Antonello comparison
#config.Data.inputDataset = '/Cosmics/Commissioning2021-PromptReco-v1/RAW-RECO' #Antonello comparison *** Maybe an AOD issue?
#config.Data.inputDataset = '/Cosmics/Commissioning2022-CosmicTP-PromptReco-v1/RAW-RECO'
config.Data.inputDataset = '/ExpressCosmics/Commissioning2022-Express-v1/FEVT' #ExpressCosmics is bigger than Cosmics
#config.Data.runRange = '348776,348773,349073'

#config.Data.runRange = '347072' #Antonello comparison

#config.Data.runRange = '342810,342966,343034,343082,343171,343266,343387,344134,344186,344266,344366'
config.Data.runRange = '348773,348776,348777,348838,348908,348955,349016,349073,349078,349079,349084,349146,349147,349260,349263,349348,349422,349433,349435,349436,349437,349527,349528,349611,349703,349758,349833,349834,349839,349840,349893,349963,350010,350107,350142,350166,350174,350254,350294,350361,350424,350425,350431,350450,350459,350460,350462,350463,350464,350490,350491,350561,350619'



#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/daebi/CRAFT2022/'
#config.Data.outLFNDirBase = '/store/group/lpcgem/'
config.Data.publication = False
#import FWCore.PythonUtilities.LumiList as LumiList
##lumiList = LumiList(filename='my_original_lumi_mask.json')
#lumiList = LumiList(filename='320887_13TeV_PromptReco_Collisions18_JSON_MuonPhys.txt')
#lumiList.selectRuns(runs = [321475, 321461,  321457,  321434,  321433,  321432,  321431,  321415,  321414,  321396,  321393,  321313,  321312,  321311, 321310,  321305,  321218,  321178,  321177,  321167,  321166,  321165,  321164,  321162,  321149,  321140,  321138,  321134,  321126, 321123,  321122,  321121,  321119,  321069,  321068,  321067,  321055,  321051,  320996,  320995])
#lumiList.writeJSON('my_lumi_mask.json')
#config.Data.lumiMask = 'my_lumi_mask.json'
#process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()
#config.Data.runRange = '%d-%d'%(runstart, runend)#'315257-315270'#'278820-278820' # '193093-194075'
config.Data.outputDatasetTag = config.General.requestName
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.ignoreGlobalBlacklist = True
#config.Site.whitelist = ["T2_KR_KISTI"]
#config.Site.whitelist = ["T0_CH_CERN_MSS"]
