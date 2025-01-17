import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('analyzer',Phase2C9)
process = cms.Process('analyzer',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_design', '')



### This is the misalignment part
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        connect = cms.string('sqlite_file:GEMposition_from_ME11_3times_firstfit_10k_each.db'),
        record = cms.string('GEMAlignmentRcd'),
        tag = cms.string('GEMAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string('sqlite_file:GEMposition_from_ME11_3times_firstfit_10k_each.db'),
        record = cms.string('GEMAlignmentErrorExtendedRcd'),
        tag = cms.string('GEMAlignmentErrorExtendedRcd')
    ),
    cms.PSet(record=cms.string('GlobalPositionRcd'), tag = cms.string('IdealGeometry'))
)

process.GEMGeometryESModule.applyAlignment = cms.bool(True)
################################



#process.CSCGeometryESModule.applyAlignment = cms.bool(False)





process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
			-1, #Max number of events 
			VarParsing.multiplicity.singleton, 
			VarParsing.varType.int, 
			"Number of events")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.nEvents)
)
process.maxEvents.input = cms.untracked.int32(-1)


process.source = cms.Source("PoolSource", 
				fileNames = cms.untracked.vstring(options.inputFiles), 
				inputCommands = cms.untracked.vstring(
			"keep *", 
			"drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO", 
			"drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
			)
				)

#process.source.fileNames.append('file:A6D8F0A8-462B-0344-91D7-35EC85B5DB6A.root')
#process.source.fileNames.append('file:step2_223.root')
process.source.fileNames.append('file:/eos/cms/store/group/alca_muonalign/singleMuonGun_11_3_4_2021_design/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3/210816_170519/0000/step2_266.root')

process.options = cms.untracked.PSet(
                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string("out_analyzer.root"))

process.analyser = cms.EDAnalyzer('analyzer', 
	process.MuonServiceProxy, 
	gemRecHits = cms.InputTag("gemRecHits"), 
	gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("muons"),
	vertexCollection = cms.InputTag("offlinePrimaryVerticies"),
        tracker_prop = cms.bool(True),
        CSC_prop = cms.bool(True),
        Segment_prop = cms.bool(True),
        debug = cms.bool(False),
        isCosmic = cms.bool(False)
)

#process.p = cms.EndPath(process.analyzer)
process.p = cms.Path(process.analyzer)
