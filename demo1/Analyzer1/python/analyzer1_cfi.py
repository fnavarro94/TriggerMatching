import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('Analyzer1'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
