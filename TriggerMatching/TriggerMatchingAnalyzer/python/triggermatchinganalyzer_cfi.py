import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('TriggerMatchingAnalyzer'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
