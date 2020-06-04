#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/RecoPFCandidateTypeAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

RecoPFCandidateTypeAnalyzer::RecoPFCandidateTypeAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , isolationQualityCuts_dzCut_disabled_(nullptr)
  , isolationQualityCuts_dzCut_enabled_primary_(nullptr)
  , pfChargedHadronPlots_(nullptr)
  , pfChargedHadronPileupPlots_(nullptr)
  , pfElectronPlots_(nullptr)
  , pfNeutralHadronPlots_(nullptr)
  , pfPhotonPlots_(nullptr)
  , pfMuonPlots_(nullptr)
{
  src_pfCands_ = cfg.getParameter<edm::InputTag>("srcPFCands");
  token_pfCands_ = consumes<reco::PFCandidateCollection>(src_pfCands_);

  srcVertices_ = cfg.getParameter<edm::InputTag>("srcVertices");
  if ( srcVertices_.label() != "" ) 
  {
    tokenVertices_ = consumes<reco::VertexCollection>(srcVertices_);
  }

  srcPileupSummaryInfo_ = cfg.getParameter<edm::InputTag>("srcPileupSummaryInfo");
  if ( srcPileupSummaryInfo_.label() != "" ) 
  {
    tokenPileupSummaryInfo_ = consumes<PileupSummaryInfoCollection>(srcPileupSummaryInfo_);
  }

  edm::ParameterSet cfg_isolationQualityCuts = cfg.getParameter<edm::ParameterSet>("isolationQualityCuts");  
  edm::ParameterSet cfg_isolationQualityCuts_dzCut_disabled(cfg_isolationQualityCuts);
  cfg_isolationQualityCuts_dzCut_disabled.addParameter<double>("maxDeltaZ", 1.e+3);
  isolationQualityCuts_dzCut_disabled_ = new reco::tau::RecoTauQualityCuts(cfg_isolationQualityCuts_dzCut_disabled);
  edm::ParameterSet cfg_isolationQualityCuts_dzCut_enabled_primary(cfg_isolationQualityCuts);
  cfg_isolationQualityCuts_dzCut_enabled_primary.addParameter<double>("maxDeltaZ", 0.2);
  isolationQualityCuts_dzCut_enabled_primary_ = new reco::tau::RecoTauQualityCuts(cfg_isolationQualityCuts_dzCut_enabled_primary);

  pfChargedHadronPlots_       = new pfCandTypePlotEntryType(-1., 1.4, "chargedHadron");
  pfChargedHadronPileupPlots_ = new pfCandTypePlotEntryType(-1., 1.4, "chargedHadronPileup");
  pfElectronPlots_            = new pfCandTypePlotEntryType(-1., 1.4, "electron");
  pfNeutralHadronPlots_       = new pfCandTypePlotEntryType(-1., 1.4, "neutralHadron");
  pfPhotonPlots_              = new pfCandTypePlotEntryType(-1., 1.4, "photon");
  pfMuonPlots_                = new pfCandTypePlotEntryType(-1., 1.4, "muon");

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

RecoPFCandidateTypeAnalyzer::~RecoPFCandidateTypeAnalyzer()
{}

void RecoPFCandidateTypeAnalyzer::beginJob()
{
  if ( !edm::Service<dqm::legacy::DQMStore>().isAvailable() ) {
    throw cms::Exception("RecoPFCandidateTypeAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<dqm::legacy::DQMStore>());

  dqmStore.setCurrentFolder(dqmDirectory_);

  pfChargedHadronPlots_->bookHistograms(dqmStore);
  pfChargedHadronPileupPlots_->bookHistograms(dqmStore);
  pfElectronPlots_->bookHistograms(dqmStore);
  pfNeutralHadronPlots_->bookHistograms(dqmStore);
  pfPhotonPlots_->bookHistograms(dqmStore);
  pfMuonPlots_->bookHistograms(dqmStore);

  me_EventCounter_ = dqmStore.book1D("EventCounter", "EventCounter", 1, -0.5, +0.5);
  histogram_EventCounter_ = me_EventCounter_->getTH1();
  assert(histogram_EventCounter_);
}
    
void RecoPFCandidateTypeAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<RecoPFCandidateTypeAnalyzer::analyze>:" << std::endl;

  edm::Handle<reco::PFCandidateCollection> pfCands;
  evt.getByToken(token_pfCands_, pfCands);

  if ( srcVertices_.label() != "" ) 
  {
    edm::Handle<reco::VertexCollection> vertices;
    evt.getByToken(tokenVertices_, vertices);
    if ( vertices->size() > 0 ) 
    {
      reco::VertexRef primaryVertexRef(vertices, 0);
      isolationQualityCuts_dzCut_disabled_->setPV(primaryVertexRef);
      isolationQualityCuts_dzCut_enabled_primary_->setPV(primaryVertexRef);
    }
  }

  int numPileup = -1;
  if ( srcPileupSummaryInfo_.label() != "" ) 
  {
    edm::Handle<PileupSummaryInfoCollection> pileupSummaryInfos;
    evt.getByToken(tokenPileupSummaryInfo_, pileupSummaryInfos);
    for ( PileupSummaryInfoCollection::const_iterator pileupSummaryInfo = pileupSummaryInfos->begin(); pileupSummaryInfo != pileupSummaryInfos->end(); ++pileupSummaryInfo ) 
    {
      if ( pileupSummaryInfo->getBunchCrossing() == 0 ) 
      {
	numPileup = pileupSummaryInfo->getPU_NumInteractions();
      }
    }
  }

  const double evtWeight = 1.;
  
  //int idxPFCand = 0;
  for ( reco::PFCandidateCollection::const_iterator pfCand = pfCands->begin(); pfCand != pfCands->end(); ++pfCand )
  {
    //std::string type_string;
    //if      ( pfCand->particleId() == reco::PFCandidate::h     ) type_string = "PFChargedHadron";
    //else if ( pfCand->particleId() == reco::PFCandidate::e     ) type_string = "PFElectron";
    //else if ( pfCand->particleId() == reco::PFCandidate::h0    ) type_string = "PFNeutralHadron";
    //else if ( pfCand->particleId() == reco::PFCandidate::gamma ) type_string = "PFPhoton";
    //else if ( pfCand->particleId() == reco::PFCandidate::mu    ) type_string = "PFMuon";
    //else                                                 type_string = "N/A";
    //std::cout << "PFCandidate #" << idxPFCand << " (type = " << type_string << "):" 
    //	        << " pT = " << l1PFCand->pt() << ", eta = " << pfCand->eta() << ", phi = " << pfCand->phi() << ","
    //	        << " isSelected (dZcut disabled) = " << isolationQualityCuts_dzCut_disabled_->filterCand(*pfCand) << std::endl;
    //++idxPFCand;
    if ( pfCand->particleId() == reco::PFCandidate::h && isolationQualityCuts_dzCut_disabled_->filterCand(*pfCand) )
    {
      if ( isolationQualityCuts_dzCut_enabled_primary_->filterCand(*pfCand) )
      {
        pfChargedHadronPlots_->fillHistograms(*pfCand, numPileup, evtWeight);
      }
      else
      {
        pfChargedHadronPileupPlots_->fillHistograms(*pfCand, numPileup, evtWeight);
      }
    }
    if ( pfCand->particleId() == reco::PFCandidate::e && isolationQualityCuts_dzCut_disabled_->filterCand(*pfCand) )
    {
      pfElectronPlots_->fillHistograms(*pfCand, numPileup, evtWeight);
    }
    if ( pfCand->particleId() == reco::PFCandidate::h0 && isolationQualityCuts_dzCut_disabled_->filterCand(*pfCand) )
    {
      pfNeutralHadronPlots_->fillHistograms(*pfCand, numPileup, evtWeight);
    }
    if ( pfCand->particleId() == reco::PFCandidate::gamma && isolationQualityCuts_dzCut_disabled_->filterCand(*pfCand) )
    {
      pfPhotonPlots_->fillHistograms(*pfCand, numPileup, evtWeight);
    }
    if ( pfCand->particleId() == reco::PFCandidate::mu && isolationQualityCuts_dzCut_disabled_->filterCand(*pfCand) )
    {
      pfMuonPlots_->fillHistograms(*pfCand, numPileup, evtWeight);
    }
  }

  histogram_EventCounter_->Fill(0., evtWeight);
}

void RecoPFCandidateTypeAnalyzer::endJob()
{
  if ( histogram_EventCounter_->Integral() > 0. )
  {
    std::cout << "<RecoPFCandidateTypeAnalyzer::endJob (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
    pfChargedHadronPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfChargedHadronPileupPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfElectronPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfNeutralHadronPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    std::cout << "scaling histogram = " << pfPhotonPlots_->histogram_ptFraction_vs_absEta_->GetName() << " by factor = " << (1./histogram_EventCounter_->Integral()) << std::endl;    
    pfPhotonPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
    pfMuonPlots_->normalizeHistograms(histogram_EventCounter_->Integral());
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(RecoPFCandidateTypeAnalyzer);


