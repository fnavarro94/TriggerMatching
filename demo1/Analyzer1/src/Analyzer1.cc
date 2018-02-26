// -*- C++ -*-
//
// Package:    Analyzer1
// Class:      Analyzer1
// 
/**\class Analyzer1 Analyzer1.cc demo1/Analyzer1/src/Analyzer1.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Fri Jan 19 03:47:24 GMT 2018
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

class Analyzer1 : public edm::EDAnalyzer {
   public:
      explicit Analyzer1(const edm::ParameterSet&);
      ~Analyzer1();

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
     int matchCount, matchCountRep,  filterCount, triggerCount, filterCount2, muonCount, correctMatches, muonMatches; // How many matches in total, how many objects passed the filter 
      
      

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
Analyzer1::Analyzer1(const edm::ParameterSet& iConfig)
:
 trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))

{
   //now do what ever initialization is needed

}


Analyzer1::~Analyzer1()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Analyzer1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   
//std::string filterName("hltDoublePhoton33EgammaLHEDoubleFilter"); // data
std::string filterName("hltL2DoubleMu30NoVertexL2PreFiltered"); // simulation 

//it is important to specify the right HLT process for the filter, not doing this is a common bug
trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
bool primeraVuelta = true;
int match=0; 
if(passTrig){
if(filterIndex<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
   

  using reco::TrackCollection;

   
    
  for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     
    // Trigger object loop starts
		for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
			if (primeraVuelta && keyIt == trigKeys.begin()){filterCount ++;}
			const trigger::TriggerObject& obj = trigObjColl[*keyIt];
			
			if((dR(obj,p)<0.1)){
				matchCount++;
				if(abs(p.pdgId())== 13)
				{muonMatches++;
				 match++;	 }
				
		}
    }
primeraVuelta = false;

}//end filter size check
      
  }
  
      
      
      
      
   }
   if (match >= 2){correctMatches ++;}
   




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
Analyzer1::beginJob()
{
	// Initialize our variables
    file = new TFile("outfile.root","recreate");
    const bool oldAddDir = TH1::AddDirectoryStatus();
    TH1::AddDirectory(true);
    histo = new TH1F("pt","pt",1000,0,100);
    TH1::AddDirectory(oldAddDir); 
	matchCount=0;
	matchCountRep=0;
	filterCount=0;
	triggerCount=0;
	filterCount2=0;
	muonCount = 0;
	correctMatches = 0;
	muonMatches=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer1::endJob() 
{   std::cout<<"Number of events that should activate trigger "<<muonCount<<std::endl;
	std::cout<<"Number of times trigger was activated "<<triggerCount<<std::endl;
    std::cout<<"Number of muon matches "<<muonMatches<<std::endl;
    std::cout<<"Number of matches: "<<matchCount<<std::endl;
	std::cout<<"Number of correct matches"  <<correctMatches<<std::endl;
	std::cout<<"Trigger efficiency 1 (times it should have been activated / times it was activated)"<< (float)muonCount/triggerCount<<std::endl;
	std::cout<<"Trigger efficiency 2 (times it was activated / times it was activated by two muons)"<< (float)triggerCount/correctMatches<<std::endl;
	
	
	
	
    file->Write();
    file->Close(); 
}

// ------------ method called when starting to processes a run  ------------

double
Analyzer1::dR(const trigger::TriggerObject& obj,  const reco::GenParticle & p )
{
	        double dEta2 =pow( p.eta()-obj.eta(),2); 
			double dPhi2 =pow( p.phi()-obj.phi(),2);
			double dR = sqrt(dPhi2+dEta2);
			return dR;
}

int
Analyzer1::countMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;
	int numMu = 0;
	int triggeringEvent = 0;
	Handle<GenParticleCollection> genParticles;
iEvent.getByLabel("genParticles", genParticles);
for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     if (abs(p.pdgId()) == 13 && p.pt() >= 30){
	  numMu++;
	 
	 }
	
}
 if(numMu >=2)
 {triggeringEvent = 1;}
 return triggeringEvent;
}
void 
Analyzer1::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Analyzer1::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analyzer1::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analyzer1::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(Analyzer1);
