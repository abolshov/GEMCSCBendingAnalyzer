#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = 'RdPhiAna_CRAFT2022_050522_ExpressCosmics_ME11Iter0_GE11Iter0'
config.General.workArea = 'RdPhiAna_CRAFT2022_ExpressCosmics_TerukiRuns'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_both_analyzers.py'
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles

misalign = False  #Make sure to change the run_both_analyzers.py too!!!
if misalign:
  config.JobType.inputFiles = ['./CRAFT_by_layer_GEMIter2.db']

#section Data
config.Data.inputDataset = '/ExpressCosmics/Commissioning2022-Express-v1/FEVT' #ExpressCosmics is bigger than Cosmics
#config.Data.runRange = '348776,348773,349073'
#config.Data.runRange = '348776,349840,348773,349839,350174,349073,348777,349422,348955,349078,349437,348838,348908,349016,348683,349893,348594,349436,349084,349348,348390,350107,349097,349611,350142,348332,349963,350010'
config.Data.runRange = '348773,348776,348777,348838,348908,348955,349016,349073,349078,349079,349084,349146,349147,349260,349263,349348,349422,349433,349435,349436,349437,349527,349528,349611,349703,349758,349833,349834,349839,349840,349893,349963,350010,350107,350142,350166,350174,350254,350294,350361,350424,350425,350431,350450,350459,350460,350462,350463,350464,350490,350491,350561,350619'

#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/daebi/CRAFT2022/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.ignoreGlobalBlacklist = True
