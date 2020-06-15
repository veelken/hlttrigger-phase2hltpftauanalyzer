#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoTrackAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"                                // JetMCTagUtils::genTauDecayMode()
#include "DataFormats/Math/interface/deltaR.h"                                         // deltaR

#include "HLTTrigger/TallinnHLTPFTauAnalyzer/interface/GenChargedHadronToTrackMatch.h" // GenChargedHadronToRecoTrackMatch

#include <algorithm> // std::sort()

enum { kQualityCuts_disabled, kQualityCuts_enabled };

RecoTrackAnalyzer::RecoTrackAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , debug_(cfg.getParameter<bool>("debug"))
{
  src_genTaus_ = cfg.getParameter<edm::InputTag>("srcGenTaus");
  token_genTaus_ = consumes<reco::GenJetCollection>(src_genTaus_);

  std::string vtxMode_string = cfg.getParameter<std::string>("vtxMode");
  if ( vtxMode_string == "genVtx" ) 
  {
    vtxMode_ = kGenVtx;

    src_genVertex_position_ = cfg.getParameter<edm::InputTag>("srcGenVertex_position");
    token_genVertex_position_ = consumes<reco::TrackBase::Point>(src_genVertex_position_);
  }
  else if ( vtxMode_string == "recVtx" ) 
  {
    vtxMode_ = kRecVtx;

    src_offlineVertices_ = cfg.getParameter<edm::InputTag>("srcOfflineVertices");
    token_offlineVertices_ = consumes<reco::VertexCollection>(src_offlineVertices_);

    src_hltVertices_ = cfg.getParameter<edm::InputTag>("srcHLTVertices");
    token_hltVertices_ = consumes<reco::VertexCollection>(src_hltVertices_);
  }
  else
  {
    throw cms::Exception("RecoTrackAnalyzer") 
      << " Invalid Configuration parameter 'vtxMode' = " << vtxMode_string << " !!\n";;
  }

  src_offlineTracks_ = cfg.getParameter<edm::InputTag>("srcOfflineTracks");
  token_offlineTracks_ = consumes<reco::TrackCollection>(src_offlineTracks_);

  src_offlinePFCands_ = cfg.getParameter<edm::InputTag>("srcOfflinePFCands");
  //token_offlinePFCands_ = consumes<reco::PFCandidateCollection>(src_offlinePFCands_);
  token_offlinePFCands_ = consumes<pat::PackedCandidateCollection>(src_offlinePFCands_);

  src_hltTracks_ = cfg.getParameter<edm::InputTag>("srcHLTTracks");
  token_hltTracks_ = consumes<reco::TrackCollection>(src_hltTracks_);

  src_hltPFCands_ = cfg.getParameter<edm::InputTag>("srcHLTPFCands");
  token_hltPFCands_ = consumes<reco::PFCandidateCollection>(src_hltPFCands_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoTrackAnalyzer::~RecoTrackAnalyzer()
{
  for ( auto efficiencyPlot : efficiencyPlots_offlineTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_offlineTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_offlinePFCandTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_offlinePFCandTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }
  
  for ( auto efficiencyPlot : efficiencyPlots_hltTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_hltTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_hltPFCandTracks_woQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }

  for ( auto efficiencyPlot : efficiencyPlots_hltPFCandTracks_wQualityCuts_ ) 
  {
    delete efficiencyPlot;
  }
}

void RecoTrackAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoTrackAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  std::vector<std::string> genTau_decayModes = { "oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "threeProng0Pi0", "threeProng1Pi0", "all" };
  std::vector<double> genChargedHadron_min_absEtaValues = { -1.,   1.4,   1.4, -1.,    -1.  };
  std::vector<double> genChargedHadron_max_absEtaValues = {  1.4,  2.172, 2.4,  2.172,  2.4 };
  assert(genChargedHadron_min_absEtaValues.size() == genChargedHadron_max_absEtaValues.size());
  size_t numAbsEtaRanges = genChargedHadron_min_absEtaValues.size();
  for ( size_t idxAbsEtaRange = 0; idxAbsEtaRange < numAbsEtaRanges; ++idxAbsEtaRange )
  {
    double genChargedHadron_min_absEta = genChargedHadron_min_absEtaValues[idxAbsEtaRange];
    double genChargedHadron_max_absEta = genChargedHadron_max_absEtaValues[idxAbsEtaRange];
    for ( auto genTau_decayMode : genTau_decayModes )
    {
      TString dqmDirectory = dqmDirectory_.data();
      if ( genChargedHadron_min_absEta >= 0. && genChargedHadron_max_absEta > 0. ) 
      { 
        dqmDirectory.Append(Form("/absEta%1.2fto%1.2f", genChargedHadron_min_absEta, genChargedHadron_max_absEta));
      }
      else if ( genChargedHadron_min_absEta >= 0. ) 
      {
        dqmDirectory.Append(Form("/absEtaGt%1.2f", genChargedHadron_min_absEta));
      }
      else if ( genChargedHadron_max_absEta > 0. ) 
      {
        dqmDirectory.Append(Form("/absEtaLt%1.2f", genChargedHadron_max_absEta));
      }
      std::string genTau_decayMode_capitalized = genTau_decayMode;
      genTau_decayMode_capitalized[0] = toupper(genTau_decayMode_capitalized[0]);	
      dqmDirectory.Append(Form("/gen%sTau", genTau_decayMode_capitalized.data()));
      dqmDirectory = dqmDirectory.ReplaceAll(".", "p");

      dqmStore.setCurrentFolder(Form("%s/offlineTrack_woQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_offlineTracks_woQualityCuts = new efficiencyPlotEntryType(
	"offlineTrack_woQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_offlineTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlineTracks_woQualityCuts_.push_back(efficiencyPlot_offlineTracks_woQualityCuts);

      dqmStore.setCurrentFolder(Form("%s/offlineTrack_wQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_offlineTracks_wQualityCuts = new efficiencyPlotEntryType(
        "offlineTrack_wQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_offlineTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlineTracks_wQualityCuts_.push_back(efficiencyPlot_offlineTracks_wQualityCuts);
      
      dqmStore.setCurrentFolder(Form("%s/offlinePFCandTrack_woQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_offlinePFCandTracks_woQualityCuts = new efficiencyPlotEntryType(
        "offlinePFCandTrack_woQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_offlinePFCandTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlinePFCandTracks_woQualityCuts_.push_back(efficiencyPlot_offlinePFCandTracks_woQualityCuts);
            
      dqmStore.setCurrentFolder(Form("%s/offlinePFCandTrack_wQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_offlinePFCandTracks_wQualityCuts = new efficiencyPlotEntryType(
        "offlinePFCandTrack_wQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_offlinePFCandTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_offlinePFCandTracks_wQualityCuts_.push_back(efficiencyPlot_offlinePFCandTracks_wQualityCuts);
      
      dqmStore.setCurrentFolder(Form("%s/hltTrack_woQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_hltTracks_woQualityCuts = new efficiencyPlotEntryType(
        "hltTrack_woQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_hltTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_hltTracks_woQualityCuts_.push_back(efficiencyPlot_hltTracks_woQualityCuts);
            
      dqmStore.setCurrentFolder(Form("%s/hltTrack_wQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_hltTracks_wQualityCuts = new efficiencyPlotEntryType(
        "hltTrack_wQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_hltTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_hltTracks_wQualityCuts_.push_back(efficiencyPlot_hltTracks_wQualityCuts);
      
      dqmStore.setCurrentFolder(Form("%s/hltPFCandTrack_woQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_hltPFCandTracks_woQualityCuts = new efficiencyPlotEntryType(
        "hltPFCandTrack_woQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_hltPFCandTracks_woQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_hltPFCandTracks_woQualityCuts_.push_back(efficiencyPlot_hltPFCandTracks_woQualityCuts);
            
      dqmStore.setCurrentFolder(Form("%s/hltPFCandTrack_wQualityCuts", dqmDirectory.Data()));
      efficiencyPlotEntryType* efficiencyPlot_hltPFCandTracks_wQualityCuts = new efficiencyPlotEntryType(
        "hltPFCandTrack_wQualityCuts", 1., 1.e+3, genChargedHadron_min_absEta, genChargedHadron_max_absEta, genTau_decayMode); 
      efficiencyPlot_hltPFCandTracks_wQualityCuts->bookHistograms(dqmStore);
      efficiencyPlots_hltPFCandTracks_wQualityCuts_.push_back(efficiencyPlot_hltPFCandTracks_wQualityCuts);
    }
  }
}

namespace
{
  // auxiliary function to select offline reconstructed tracks of "good quality" 
  bool passesRecoTrackQualityCuts(const reco::Track& track, const reco::TrackBase::Point& primaryVertex_position)
  {    
    bool passesQualityCuts = ( TMath::Abs(track.dz(primaryVertex_position))  <   0.4  &&
			       track.hitPattern().numberOfValidPixelHits()   >=  0    &&
			       track.hitPattern().numberOfValidHits()        >=  3    &&
			       TMath::Abs(track.dxy(primaryVertex_position)) <   0.03 &&
			       track.normalizedChi2()                        < 100.   ) ? true : false;
    return passesQualityCuts;
  }

  bool passesRecoTrackQualityCuts(const reco::Track& track, const reco::Vertex* primaryVertex)
  {    
    bool passesQualityCuts = ( primaryVertex && passesRecoTrackQualityCuts(track, primaryVertex->position()) ) ? true : false;
    return passesQualityCuts;
  }

  // auxiliary function for matches between generator-level charged hadrons produced in tau decays and offline reconstructed tracks
  std::vector<GenChargedHadronToRecoTrackMatch_and_genTau_decayMode> 
  getGenChargedHadronToRecoTrackMatches(const std::vector<GenChargedHadron_and_genTau_decayMode>& genTauChargedHadrons, std::vector<const reco::Track*> recTracks, bool debug)
  {
    std::vector<GenChargedHadronToRecoTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches;
    const double dRmatch = 0.05;
    size_t idxGenTauChargedHadron = 0;
    for ( auto genTauChargedHadron : genTauChargedHadrons )
    {
      size_t idxRecTrack = 0;
      for ( auto recTrack : recTracks )
      {
	double dR = reco::deltaR(genTauChargedHadron.genChargedHadron_eta(), genTauChargedHadron.genChargedHadron_phi(), recTrack->eta(), recTrack->phi());
	if ( dR < dRmatch ) 
        {
	  if ( debug ) 
	  {
	    std::cout << "matching genChargedHadron (#" << idxGenTauChargedHadron << "):"
		      << " pT = " << genTauChargedHadron.genChargedHadron_pt() << ","
		      << " eta = " << genTauChargedHadron.genChargedHadron_eta() << ","
		      << " phi = " << genTauChargedHadron.genChargedHadron_phi() << ","
		      << " to recTrack (#" << idxRecTrack << "):" 
		      << " pT = " << recTrack->pt() << ","
		      << " eta = " << recTrack->eta() << ","
		      << " phi = " << recTrack->phi() 
		      << " (dR = " << dR << ")" << std::endl;
	  }
  	  genChargedHadronToTrackMatches.push_back(GenChargedHadronToRecoTrackMatch_and_genTau_decayMode(
	    GenChargedHadronToRecoTrackMatch(genTauChargedHadron.genChargedHadron(), recTrack),
	    genTauChargedHadron.genTau_decayMode(),
	    dR));
	}
	++idxRecTrack;
      }

      // CV: add "empty" match to allow for genParticles without recTrack match (tracking inefficiency)
      genChargedHadronToTrackMatches.push_back(GenChargedHadronToRecoTrackMatch_and_genTau_decayMode(
	GenChargedHadronToRecoTrackMatch(genTauChargedHadron.genChargedHadron(), nullptr), 
	genTauChargedHadron.genTau_decayMode(),
	1.e+3));

      ++idxGenTauChargedHadron;
    }
    return genChargedHadronToTrackMatches;
  }

  // auxiliary function to sort (ranke) matches (needed for cleaning)
  template<class T>
  bool isLowerDeltaR(const T* genChargedHadronToTrackMatch1, const T* genChargedHadronToTrackMatch2)
  {
    return (genChargedHadronToTrackMatch1->dR() < genChargedHadronToTrackMatch2->dR());
  }

  // auxiliary function to clean matches (in order to guarantee that each generator-level charged hadron and each reconstructed track is used in the matching only once)
  template <class T>
  std::vector<const T*> 
  cleanGenChargedHadronToTrackMatches(const std::vector<T>& genChargedHadronToTrackMatches, bool debug)
  {
    std::list<const T*> matches_remaining;
    for ( typename std::vector<T>::const_iterator match = genChargedHadronToTrackMatches.begin();
	  match != genChargedHadronToTrackMatches.end(); ++match )
    {
      matches_remaining.push_back(&(*match));
    }

    std::vector<const T*> matches_selected;
    while ( matches_remaining.size() > 0 ) 
    {
      matches_remaining.sort(isLowerDeltaR<T>);
      const T* bestMatch = matches_remaining.front();
      assert(bestMatch);
      if ( debug )
      {
        std::cout << "bestMatch (#" << matches_selected.size() << "):" 
	  	  << " pT = " << bestMatch->genChargedHadron_pt() << ","
		  << " eta = " << bestMatch->genChargedHadron_eta() << ","
		  << " phi = " << bestMatch->genChargedHadron_phi() << std::endl;
      }
      matches_selected.push_back(bestMatch);
      for ( auto match : matches_remaining ) 
      {
        if ( isOverlap(match->genChargedHadronToTrackMatch(), bestMatch->genChargedHadronToTrackMatch()) )
	{
          matches_remaining.remove(match);
        }
      }      
    }
    return matches_selected;
  }
}
    
void RecoTrackAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<reco::GenJetCollection> genTaus;
  evt.getByToken(token_genTaus_, genTaus);

  std::vector<GenChargedHadron_and_genTau_decayMode> genTauChargedHadrons;
  for ( reco::GenJetCollection::const_iterator genTau = genTaus->begin(); genTau != genTaus->end(); ++genTau ) 
  {
    std::string genTau_decayMode = JetMCTagUtils::genTauDecayMode(*genTau);
    std::vector<const reco::GenParticle*> genTau_daughters = genTau->getGenConstituents();
    for ( auto genTau_daughter : genTau_daughters ) 
    {
      if ( genTau_daughter->charge() != 0 ) 
      {
	genTauChargedHadrons.push_back(GenChargedHadron_and_genTau_decayMode(genTau_daughter, genTau_decayMode));
      }
    }
  }

  const double evtWeight = 1.;

  //-----------------------------------------------------------------------------
  // process generator-level event vertex
  edm::Handle<reco::TrackBase::Point> genPrimaryVertex_position;
  if ( vtxMode_ == kGenVtx ) 
  {
    evt.getByToken(token_genVertex_position_, genPrimaryVertex_position);    
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process offline reconstructed vertices
  const reco::Vertex* offlinePrimaryVertex = nullptr;
  if ( vtxMode_ == kRecVtx && src_offlineVertices_.label() != "" ) 
  {
    edm::Handle<reco::VertexCollection> offlineVertices;
    evt.getByToken(token_offlineVertices_, offlineVertices);    
    if ( offlineVertices->size() > 0 ) offlinePrimaryVertex = &offlineVertices->at(0);
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process offline reconstructed tracks
  if ( src_offlineTracks_.label() != "" )
  {
    edm::Handle<reco::TrackCollection> offlineTracks;
    evt.getByToken(token_offlineTracks_, offlineTracks);

    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const reco::Track*> selectedOfflineTracks;
      for ( reco::TrackCollection::const_iterator recTrack = offlineTracks->begin(); recTrack != offlineTracks->end(); ++recTrack )
      {
        if (  idxQualityCuts == kQualityCuts_disabled                                                   ||
	     (vtxMode_ == kGenVtx && passesRecoTrackQualityCuts(*recTrack, *genPrimaryVertex_position)) ||
	     (vtxMode_ == kRecVtx && passesRecoTrackQualityCuts(*recTrack, offlinePrimaryVertex))       )
	{
	  selectedOfflineTracks.push_back(&(*recTrack));
	}
      }

      if ( debug_ )
      {
	size_t idx = 0;
	for ( auto recTrack : selectedOfflineTracks )
	{
	  std::cout << "recTrack #" << idx << ":" 
		    << " pT = " << recTrack->pt() << ", eta = " << recTrack->eta() << ", phi = " << recTrack->phi() << std::endl;
	  ++idx;
	}
      }
      
      std::vector<GenChargedHadronToRecoTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToRecoTrackMatches(
	genTauChargedHadrons,
	selectedOfflineTracks,
	debug_);   
      if ( debug_ )
      {
	std::cout << "genChargedHadronToTrackMatches BEFORE cleaning:" << std::endl;
	printGenChargedHadronToTrackMatches(genChargedHadronToTrackMatches);
      }

      std::vector<const GenChargedHadronToRecoTrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches,
	debug_); 
      if ( debug_ )
      {
	std::cout << "genChargedHadronToTrackMatches AFTER cleaning:" << std::endl;
	printGenChargedHadronToTrackMatches(cleanedGenChargedHadronToTrackMatches);
      }

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled ) efficiencyPlots = &efficiencyPlots_offlineTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  ) efficiencyPlots = &efficiencyPlots_offlineTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process tracks associated to offline reconstructed PF candidates
  if ( src_offlinePFCands_.label() != "" && src_offlineTracks_.label() != "" )
  {
    //edm::Handle<reco::PFCandidateCollection> offlinePFCands;
    edm::Handle<pat::PackedCandidateCollection> offlinePFCands;
    evt.getByToken(token_offlinePFCands_, offlinePFCands);
  
    edm::Handle<reco::TrackCollection> offlineTracks;
    evt.getByToken(token_offlineTracks_, offlineTracks);

    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const reco::Track*> selectedOfflineTracks;
      //for ( reco::PFCandidateCollection::const_iterator recPFCand = offlinePFCands->begin(); recPFCand != offlinePFCands->end(); ++recPFCand )
      for ( pat::PackedCandidateCollection::const_iterator recPFCand = offlinePFCands->begin(); recPFCand != offlinePFCands->end(); ++recPFCand )
      {
	if ( recPFCand->charge() != 0 )
	{
	  //const reco::Track* recTrack = nullptr;
	  //if      ( recPFCand->trackRef().isNonnull()    ) recTrack = recPFCand->trackRef().get();
	  //else if ( recPFCand->gsfTrackRef().isNonnull() ) recTrack = recPFCand->gsfTrackRef().get();
	  //if ( !recTrack ) continue;
          const reco::Track* recTrack_matched = nullptr;
          double dRmin = 1.e+3;
          for ( reco::TrackCollection::const_iterator recTrack = offlineTracks->begin(); recTrack != offlineTracks->end(); ++recTrack )
          {
            double dR = reco::deltaR(recPFCand->eta(), recPFCand->phi(), recTrack->eta(), recTrack->phi());
            if ( dR < 1.e-2 && dR < dRmin ) 
            {
              recTrack_matched = &(*recTrack);
              dRmin = dR;
            }
          }
          if ( !recTrack_matched ) continue;

          if (  idxQualityCuts == kQualityCuts_disabled                                                           ||
	       (vtxMode_ == kGenVtx && passesRecoTrackQualityCuts(*recTrack_matched, *genPrimaryVertex_position)) ||
	       (vtxMode_ == kRecVtx && passesRecoTrackQualityCuts(*recTrack_matched, offlinePrimaryVertex))       )
  	  {
	    selectedOfflineTracks.push_back(&(*recTrack_matched));
	  }
	}
      }
      
      std::vector<GenChargedHadronToRecoTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToRecoTrackMatches(
	genTauChargedHadrons,
	selectedOfflineTracks,
	debug_);     
      std::vector<const GenChargedHadronToRecoTrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches,
	debug_);     

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled ) efficiencyPlots = &efficiencyPlots_offlinePFCandTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  ) efficiencyPlots = &efficiencyPlots_offlinePFCandTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process tracks reconstructed on HLT level
  if ( src_hltTracks_.label() != "" )
  {
    const reco::Vertex* hltPrimaryVertex = nullptr;
    if ( vtxMode_ == kRecVtx && src_hltVertices_.label() != "" ) 
    {
      edm::Handle<reco::VertexCollection> hltVertices;
      evt.getByToken(token_hltVertices_, hltVertices);    

      if ( hltVertices->size() > 0 ) hltPrimaryVertex = &hltVertices->at(0);
    }

    edm::Handle<reco::TrackCollection> hltTracks;
    evt.getByToken(token_hltTracks_, hltTracks);

    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const reco::Track*> selectedHLTTracks;
      for ( reco::TrackCollection::const_iterator recTrack = hltTracks->begin(); recTrack != hltTracks->end(); ++recTrack )
      {
        if (  idxQualityCuts == kQualityCuts_disabled                                                   ||
	     (vtxMode_ == kGenVtx && passesRecoTrackQualityCuts(*recTrack, *genPrimaryVertex_position)) ||
	     (vtxMode_ == kRecVtx && passesRecoTrackQualityCuts(*recTrack, hltPrimaryVertex))           )
	{
	  selectedHLTTracks.push_back(&(*recTrack));
	}
      }
      
      std::vector<GenChargedHadronToRecoTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToRecoTrackMatches(
        genTauChargedHadrons,
	selectedHLTTracks,
	debug_);     
      std::vector<const GenChargedHadronToRecoTrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches,
	debug_);     

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled ) efficiencyPlots = &efficiencyPlots_hltTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  ) efficiencyPlots = &efficiencyPlots_hltTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // process tracks associated to PF candidates reconstructed on L1 trigger level
  if ( src_hltPFCands_.label() != "" )
  {    
    const reco::Vertex* hltPrimaryVertex = nullptr;
    if ( vtxMode_ == kRecVtx && src_hltVertices_.label() != "" ) 
    {
      edm::Handle<reco::VertexCollection> hltVertices;
      evt.getByToken(token_hltVertices_, hltVertices);    

      if ( hltVertices->size() > 0 ) hltPrimaryVertex = &hltVertices->at(0);
    }

    edm::Handle<reco::PFCandidateCollection> hltPFCands;
    evt.getByToken(token_hltPFCands_, hltPFCands);
  
    for ( int idxQualityCuts = kQualityCuts_disabled; idxQualityCuts <= kQualityCuts_enabled; ++idxQualityCuts )
    {
      std::vector<const reco::Track*> selectedHLTTracks;
      for ( reco::PFCandidateCollection::const_iterator recPFCand = hltPFCands->begin(); recPFCand != hltPFCands->end(); ++recPFCand )
      {
	if ( recPFCand->charge() != 0 )
	{
	  const reco::Track* recTrack = nullptr;
	  if      ( recPFCand->trackRef().isNonnull()    ) recTrack = recPFCand->trackRef().get();
	  else if ( recPFCand->gsfTrackRef().isNonnull() ) recTrack = recPFCand->gsfTrackRef().get();
	  if ( !recTrack ) continue;

          if (  idxQualityCuts == kQualityCuts_disabled                                                   ||
	       (vtxMode_ == kGenVtx && passesRecoTrackQualityCuts(*recTrack, *genPrimaryVertex_position)) ||
	       (vtxMode_ == kRecVtx && passesRecoTrackQualityCuts(*recTrack, hltPrimaryVertex))       )
  	  {
	    selectedHLTTracks.push_back(&(*recTrack));
	  }
	}
      }
      
      std::vector<GenChargedHadronToRecoTrackMatch_and_genTau_decayMode> genChargedHadronToTrackMatches = getGenChargedHadronToRecoTrackMatches(
	genTauChargedHadrons,
	selectedHLTTracks,
	debug_);     
      std::vector<const GenChargedHadronToRecoTrackMatch_and_genTau_decayMode*> cleanedGenChargedHadronToTrackMatches = cleanGenChargedHadronToTrackMatches(
        genChargedHadronToTrackMatches,
	debug_);     

      const std::vector<efficiencyPlotEntryType*>* efficiencyPlots = nullptr;
      if      ( idxQualityCuts == kQualityCuts_disabled ) efficiencyPlots = &efficiencyPlots_hltPFCandTracks_woQualityCuts_;
      else if ( idxQualityCuts == kQualityCuts_enabled  ) efficiencyPlots = &efficiencyPlots_hltPFCandTracks_wQualityCuts_;
      else assert(0);
      for ( auto efficiencyPlot : *efficiencyPlots )
      {
	efficiencyPlot->fillHistograms(cleanedGenChargedHadronToTrackMatches, evtWeight);
      }
    }
  }
  //-----------------------------------------------------------------------------
}

void RecoTrackAnalyzer::endJob()
{}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoTrackAnalyzer);


