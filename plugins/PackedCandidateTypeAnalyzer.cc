#include "HLTTrigger/TallinnHLTPFTauAnalyzer/plugins/PackedCandidateTypeAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TMath.h> // TMath::Abs()

PackedCandidateTypeAnalyzer::PackedCandidateTypeAnalyzer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label"))
  , isolationQualityCuts_dzCut_disabled_(nullptr)
  , isolationQualityCuts_dzCut_enabled_primary_(nullptr)
  , applyPuppiWeights_(cfg.getParameter<bool>("applyPuppiWeights"))
  , pfChargedHadronPlots_(nullptr)
  , pfChargedHadronPileupPlots_(nullptr)
  , pfElectronPlots_(nullptr)
  , pfNeutralHadronPlots_(nullptr)
  , pfPhotonPlots_(nullptr)
  , pfMuonPlots_(nullptr)
{
  src_packedCands_ = cfg.getParameter<edm::InputTag>("srcPackedCands");
  token_packedCands_ = consumes<pat::PackedCandidateCollection>(src_packedCands_);

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

  pfChargedHadronPlots_       = new pfCandTypePlotEntryType("chargedHadron",       applyPuppiWeights_);
  pfChargedHadronPileupPlots_ = new pfCandTypePlotEntryType("chargedHadronPileup", applyPuppiWeights_);
  pfElectronPlots_            = new pfCandTypePlotEntryType("electron",            applyPuppiWeights_);
  pfNeutralHadronPlots_       = new pfCandTypePlotEntryType("neutralHadron",       applyPuppiWeights_);
  pfPhotonPlots_              = new pfCandTypePlotEntryType("photon",              applyPuppiWeights_);
  pfMuonPlots_                = new pfCandTypePlotEntryType("muon",                applyPuppiWeights_);

  dqmDirectory_ = cfg.getParameter<std::string>("dqmDirectory");
}

PackedCandidateTypeAnalyzer::~PackedCandidateTypeAnalyzer()
{}

void PackedCandidateTypeAnalyzer::beginJob()
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    throw cms::Exception("PackedCandidateTypeAnalyzer") 
      << " Failed to access dqmStore --> histograms will NEITHER be booked NOR filled !!\n";
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());

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

namespace
{
  const int pdgId_chargedHadron = 211;
  const int pdgId_electron      =  11;
  const int pdgId_neutralHadron = 130;
  const int pdgId_photon        =  22;
  const int pdgId_muon          =  13;
}

void PackedCandidateTypeAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<pat::PackedCandidateCollection> packedCands;
  evt.getByToken(token_packedCands_, packedCands);

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
  
  for ( pat::PackedCandidateCollection::const_iterator packedCand = packedCands->begin(); packedCand != packedCands->end(); ++packedCand )
  {    
    int packedCand_absPdgId = TMath::Abs(packedCand->pdgId());
    if ( packedCand_absPdgId == pdgId_chargedHadron && isolationQualityCuts_dzCut_disabled_->filterCand(*packedCand) )
    {
      if ( isolationQualityCuts_dzCut_enabled_primary_->filterCand(*packedCand) )
      {
        pfChargedHadronPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
      }
      else
      {
        pfChargedHadronPileupPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
      }
    }
    if ( packedCand_absPdgId == pdgId_electron && isolationQualityCuts_dzCut_disabled_->filterCand(*packedCand) )
    {
      pfElectronPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_neutralHadron && isolationQualityCuts_dzCut_disabled_->filterCand(*packedCand) )
    {
      pfNeutralHadronPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_photon && isolationQualityCuts_dzCut_disabled_->filterCand(*packedCand) )
    {
      pfPhotonPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
    if ( packedCand_absPdgId == pdgId_muon && isolationQualityCuts_dzCut_disabled_->filterCand(*packedCand) )
    {
      pfMuonPlots_->fillHistograms(*packedCand, numPileup, evtWeight);
    }
  }

  histogram_EventCounter_->Fill(0., evtWeight);
}

void PackedCandidateTypeAnalyzer::endJob()
{
  if ( histogram_EventCounter_->Integral() > 0. )
  {
    std::cout << "<PackedCandidateTypeAnalyzer::endJob (moduleLabel = " << moduleLabel_ << ")>:" << std::endl;
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

DEFINE_FWK_MODULE(PackedCandidateTypeAnalyzer);


