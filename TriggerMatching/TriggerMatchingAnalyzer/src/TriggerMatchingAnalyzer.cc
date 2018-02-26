// -*- C++ -*-
//
// Package:    TriggerMatchingAnalyzer
// Class:      TriggerMatchingAnalyzer
// 
/**\class TriggerMatchingAnalyzer TriggerMatchingAnalyzer.cc TriggerMatching/TriggerMatchingAnalyzer/src/TriggerMatchingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon Feb 26 04:40:10 CLST 2018
// $Id$
//
//

// include root
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
// system include files
#include <memory>

// For trigger object
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//
// class declaration
//

class TriggerMatchingAnalyzer : public edm::EDAnalyzer {
 public:
      explicit TriggerMatchingAnalyzer(const edm::ParameterSet&);
      ~TriggerMatchingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      int countMuons(const edm::Event&, const edm::EventSetup&);
      double dR(const trigger::TriggerObject& ,  const reco::GenParticle &  );
      // Add a Root TFiel and a Histogram
      
     TFile * file;
     TH1F * histo; 
     int  triggerCount, muonCount, correctActivations; // How many matches in total, how many objects passed the filter 
      
      

      // ----------member data ---------------------------
      edm::InputTag trackTags_; //used to select what tracks to read from configuration file
     };

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TriggerMatchingAnalyzer::TriggerMatchingAnalyzer(const edm::ParameterSet& iConfig)
:
 trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))

{
   //now do what ever initialization is needed

}


TriggerMatchingAnalyzer::~TriggerMatchingAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerMatchingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
   using namespace reco;
   
   
   
 muonCount = muonCount + countMuons(iEvent, iSetup);  
   
// Trigger activiation count

edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
//data process=HLT, MC depends, Spring11 is REDIGI311X
iEvent.getByLabel(trigResultsTag,trigResults);
const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);   

//std::string pathName="HLT_DoublePhoton33_v2"; //data
std::string pathName="HLT_L2DoubleMu30_NoVertex_v5";  // simulation


bool passTrig=trigResults->accept(trigNames.triggerIndex(pathName)); 

// count how many tymes the trigger was activated
if (passTrig)
{triggerCount ++;}
  
// *** Trigger Matching   
InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); 
Handle<trigger::TriggerEvent> trigEvent; 
iEvent.getByLabel(trigEventTag,trigEvent);

Handle<GenParticleCollection> genParticles;
iEvent.getByLabel("genParticles", genParticles);
   
std::string filterName("hltL2DoubleMu30NoVertexL2PreFiltered"); // Filter corresponding to trigger HLT_L2DoubleMu30_NoVertex_v5
trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 

int match=0; 
if(passTrig)
{ // 
    if(filterIndex<trigEvent->sizeFilters())
    { 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        using reco::TrackCollection;
    
        for(size_t i = 0; i < genParticles->size(); ++ i) 
        {
            const GenParticle & p = (*genParticles)[i];
     
            // Trigger object loop starts
		    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt)
		    { 
			    const trigger::TriggerObject& obj = trigObjColl[*keyIt];
			    if((dR(obj,p)<0.1))
			    {
			        if(abs(p.pdgId())== 13)
				    {
					    match++;
				    }	
		        }
            }
         }   
     }
}
if (match >= 2){correctActivations ++;} // if at least two muons passed the filter add a count to correct activations
   




#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TriggerMatchingAnalyzer::beginJob()
{
	// Initialize our variables
    file = new TFile("outfile.root","recreate");
    const bool oldAddDir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(true);
    histo = new TH1F("pt","pt",1000,0,100);
    TH1::AddDirectory(oldAddDir); 
	triggerCount=0; // Counts the number of times the trigger was activated
	muonCount = 0;  // Counts the number of events that should activate the trigger
	correctActivations = 0; // Counts the nuber of times the trigger was activated by the correct object
	
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerMatchingAnalyzer::endJob() 
{
	 std::cout<<"Number of events that should activate trigger "<<muonCount<<std::endl;
	std::cout<<"Number of times trigger was activated "<<triggerCount<<std::endl;
	std::cout<<"Number of correct activations "  <<correctActivations<<std::endl;
	std::cout<<"Trigger efficiency 1 (times it was activated/times it should have been activated)"<< (float)triggerCount/muonCount<<std::endl;
	std::cout<<"Trigger efficiency 2 (times it was activated / times it was activated by two muons)"<< (float)triggerCount/correctActivations<<std::endl;
	
	
	
	
    file->Write();
    file->Close(); 
}

double
TriggerMatchingAnalyzer::dR(const trigger::TriggerObject& obj,  const reco::GenParticle & p )
{
	        double dEta2 =pow( p.eta()-obj.eta(),2); 
			double dPhi2 =pow( p.phi()-obj.phi(),2);
			double dR = sqrt(dPhi2+dEta2);
			return dR;
}

int
TriggerMatchingAnalyzer::countMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;
	int numMu = 0;
	int triggeringEvent = 0;
	Handle<GenParticleCollection> genParticles;
iEvent.getByLabel("genParticles", genParticles);
for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     if (abs(p.pdgId()) == 13 && p.pt() > 30){
	  numMu++;
	 
	 }
	
}
 if(numMu >=2)
 {triggeringEvent = 1;}
 return triggeringEvent;
}
// ------------ method called when starting to processes a run  ------------
void 
TriggerMatchingAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TriggerMatchingAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TriggerMatchingAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TriggerMatchingAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerMatchingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

 //Specify that only 'tracks' is allowed
 //To use, remove the default given above and uncomment below
 //ParameterSetDescription desc;
 //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
 //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerMatchingAnalyzer);
