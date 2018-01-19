import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/Photon/AOD/12Oct2013-v1/10000/AC980FCC-943D-E311-A1D6-003048F11824.root'
    )
)

process.demo = cms.EDAnalyzer('Analyzer1'
    , tracks = cms.untracked.InputTag('generalTracks')
)


process.p = cms.Path(process.demo)
