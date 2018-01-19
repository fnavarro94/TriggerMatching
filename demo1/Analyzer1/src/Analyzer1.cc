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


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
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
      
      // Add a Root TFiel and a Histogram
      
     TFile * file;
     TH1F * histo; 
     int matchCount, matchCountRep,  filterCount; // How many matches in total, how many objects passed the filter 
      
      

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
   
   
   
   
InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
//data process=HLT, MC depends, Spring11 is REDIGI311X
Handle<trigger::TriggerEvent> trigEvent; 
iEvent.getByLabel(trigEventTag,trigEvent);

std::string filterName("hltDoublePhoton33EgammaLHEDoubleFilter"); 

//it is important to specify the right HLT process for the filter, not doing this is a common bug
trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
bool primeraVuelta = true;
if(filterIndex<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
   

  using reco::TrackCollection;

   Handle<TrackCollection> tracks;
   iEvent.getByLabel(trackTags_,tracks);
   for(TrackCollection::const_iterator itTrack = tracks->begin();
       itTrack != tracks->end();                      
       ++itTrack) {
      
  
     
      
      bool check = true;
    
    // Trigger object loop starts
		for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
			if (primeraVuelta){filterCount ++;}
			const trigger::TriggerObject& obj = trigObjColl[*keyIt];
			double dEta2 =pow( itTrack->eta()-obj.eta(),2); 
			double dPhi2 =pow( itTrack->phi()-obj.phi(),2);
			double dR = sqrt(dPhi2+dEta2);
			if((dR<0.1)&&(abs(itTrack->pt() - obj.pt()) < 3)){
				
				if (check){
				matchCount++;}
				check = false;
				matchCountRep++;
				
		}
    }
primeraVuelta = false;
}//end filter size check
      
      
      
      
      
      
   }

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
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analyzer1::endJob() 
{
	std::cout<<"Number of objects passing filter: "<<filterCount<<std::endl;
	std::cout<<"Number of matches: "<<matchCount<<std::endl;
	std::cout<<"Number of matches with rep: "<<matchCount<<std::endl;
    file->Write();
    file->Close(); 
}

// ------------ method called when starting to processes a run  ------------
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