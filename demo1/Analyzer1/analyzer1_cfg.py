import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         #'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/ZZTo2e2mu_7TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S13_START53_LV6-v1/020000/30FDADFE-4C9D-E411-977F-00266CFFBF30.root'
    #'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/ZZTo2mu2tau_7TeV_mll8_mZZ95-160-powheg15-pythia6/AODSIM/PU_S13_START53_LV6-v1/10000/065161D2-5692-E411-A439-0025905964CC.root'
    'root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/VBFToHToZZTo2L2Nu_M-450_7TeV-powheg15-pythia6/AODSIM/PU_S13_START53_LV6-v1/00000/0AE64C55-1FC8-E311-9E77-002481E15008.root'
    )
)

process.demo = cms.EDAnalyzer('Analyzer1'
    , tracks = cms.untracked.InputTag('generalTracks')
)


process.p = cms.Path(process.demo)
