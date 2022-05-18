import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('analyzer',Phase2C9)
process = cms.Process('GEMCSCanalyzer',Run3)

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


### This is the misalignment part
misalign = False
if misalign:
  db_file = 'sqlite_file:CRAFT_by_layer_GEMIter2.db'
  process.GlobalTag.toGet = cms.VPSet(
    #GE11 rec/tag
    cms.PSet(
        connect = cms.string(db_file),
        record = cms.string('GEMAlignmentRcd'),
        tag = cms.string('GEMAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string(db_file),
        record = cms.string('GEMAlignmentErrorExtendedRcd'),
        tag = cms.string('GEMAlignmentErrorExtendedRcd')
    ),
    #ME11 rec/tag
    cms.PSet(
        connect = cms.string(db_file),
        record = cms.string('CSCAlignmentRcd'),
        tag = cms.string('CSCAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string(db_file),
        record = cms.string('CSCAlignmentErrorExtendedRcd'),
        tag = cms.string('CSCAlignmentErrorExtendedRcd')
    ),
    cms.PSet(record=cms.string('GlobalPositionRcd'), tag = cms.string('IdealGeometry'))
  )

  process.GEMGeometryESModule.applyAlignment = cms.bool(True)
  process.CSCGeometryESModule.applyAlignment = cms.bool(True)
################################


#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_design', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_prompt', '') #Antonello Comparison

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

testfile = "/eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/348/776/00000/475b2a2f-673c-4104-a360-72ddee06377f.root"
testfile = 'root://eoscms.cern.ch//eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/349/347/00000/304f8f34-433f-4ea9-9df7-e32f769ad904.root'
testfile = "013a0b3d-c139-4f5d-baaa-d7bc01ef886b.root"
testfile = "step2_400.root"
#testfile = "/eos/cms/store/express/Commissioning2022/ExpressCosmics/FEVT/Express-v1/000/350/107/00000/180b850e-f936-4ae4-b305-da24855658b2.root" #One of teruki's listed runs for CRAFT
outfile = "CSC_tbma.root"
#outfile = "CSC_tbma_COSMIC.root"
process.source.fileNames.append('file:'+testfile)

process.options = cms.untracked.PSet(
                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) #variable name set above

process.CSC_tbma = cms.EDAnalyzer('CSC_tbma',
	process.MuonServiceProxy,
        csc2DRecHits = cms.InputTag("csc2DRecHits"),
        muons = cms.InputTag("muons"),
	vertexCollection = cms.InputTag("offlinePrimaryVerticies"),
        debug = cms.bool(False),
        isCosmic = cms.bool(True)
)




process.p = cms.Path(process.CSC_tbma)
