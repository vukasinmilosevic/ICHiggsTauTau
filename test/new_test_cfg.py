import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")
import sys

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('file:/Volumes/HDD/DYJetsToLL.root')
)
process.GlobalTag.globaltag = cms.string('START53_V22::All')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

##############################################################################
## Selection
##############################################################################
process.selectedVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("ndof >= 4 & abs(z) <= 24 & abs(position.Rho) <= 2")
)

process.selectedElectrons = cms.EDFilter("GsfElectronSelector",
  src = cms.InputTag("gsfElectrons"),
  cut = cms.string("pt > 20 & abs(eta) < 2.5")
)

process.selectedMuons = cms.EDFilter("MuonSelector",
  src = cms.InputTag("muons"),
  cut = cms.string("pt > 20 & abs(eta) < 2.1")
)

process.selectedPFMuons = cms.EDFilter("GenericPFCandidateSelector",
  src = cms.InputTag("particleFlow"),
  cut = cms.string("pt > 20 & abs(eta) < 2.1 & abs(pdgId) == 13")
)

process.selectedCaloJets = cms.EDFilter("CaloJetRefSelector",
  src = cms.InputTag("ak5CaloJets"),
  cut = cms.string("abs(eta) < 4.7")
)

process.selectedPFJets = cms.EDFilter("PFJetRefSelector",
  src = cms.InputTag("ak5PFJets"),
  cut = cms.string("abs(eta) < 4.7")
)

##############################################################################
## PF Lepton Isolation
##############################################################################
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'selectedElectrons')
process.muoIsoSequence = setupPFMuonIso(process, 'selectedPFMuons')

##############################################################################
## Jet Parton Flavour
##############################################################################
# This is the default Run1 configuration "Algorithmic" used by b-tag POG.
# Also what is done by PAT. New recipe and details for Run 2:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#Legacy_jet_flavour_definition
# See PAT implementation here:
# https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/PhysicsTools/PatAlgos/plugins/PATJetProducer.cc#L296-L298
process.jetPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)

process.caloJetPartonMatches = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5CaloJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("jetPartons")
)
process.caloJetFlavourAssociation = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("caloJetPartonMatches"),
    physicsDefinition = cms.bool(False)
)

process.pfJetPartonMatches = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5PFJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("jetPartons")
)
process.pfJetFlavourAssociation = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("pfJetPartonMatches"),
    physicsDefinition = cms.bool(False)
)
##############################################################################
## Jet Extras
##############################################################################
process.icCaloJetFlavourCalculator = cms.EDProducer('ICJetFlavourCalculator',
    input       = cms.InputTag("ak5CaloJets"),
    flavourMap  = cms.InputTag("caloJetFlavourAssociation")
)

process.icPFJetFlavourCalculator = cms.EDProducer('ICJetFlavourCalculator',
    input       = cms.InputTag("ak5PFJets"),
    flavourMap  = cms.InputTag("pfJetFlavourAssociation")
)


process.ak5CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    srcRho = cms.InputTag("kt6CaloJets","rho"),
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet')
)
process.ak5CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)
process.ak5CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)

process.ak5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    srcRho = cms.InputTag("kt6PFJets","rho"),
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet')
)
process.ak5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2Relative')
)
process.ak5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L3Absolute')
)

##############################################################################
## Jet B-tagging
##############################################################################
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorAtVertex
process.load("RecoBTag.Configuration.RecoBTag_cff")
import RecoBTag.Configuration.RecoBTag_cff as btag
process.ak5JetTracksAssociatorAtVertexAK5PF = ak5JetTracksAssociatorAtVertex.clone(
  jets = cms.InputTag("ak5PFJets")
  )
process.impactParameterTagInfosAK5PF = btag.impactParameterTagInfos.clone(
  jetTracks = cms.InputTag('ak5JetTracksAssociatorAtVertexAK5PF')
  )
process.secondaryVertexTagInfosAK5PF = btag.secondaryVertexTagInfos.clone(
  trackIPTagInfos = cms.InputTag('impactParameterTagInfosAK5PF')
  )
process.simpleSecondaryVertexHighEffBJetTagsAK5PF = btag.simpleSecondaryVertexHighEffBJetTags.clone (
  tagInfos = cms.VInputTag('secondaryVertexTagInfosAK5PF')
  )
process.simpleSecondaryVertexHighPurBJetTagsAK5PF = btag.simpleSecondaryVertexHighPurBJetTags.clone (
  tagInfos = cms.VInputTag('secondaryVertexTagInfosAK5PF')
  )
process.combinedSecondaryVertexBJetTagsAK5PF = btag.combinedSecondaryVertexBJetTags.clone (
  tagInfos = cms.VInputTag('impactParameterTagInfosAK5PF','secondaryVertexTagInfosAK5PF')
  )
process.btaggingSequenceAK5PF = cms.Sequence(
  process.ak5JetTracksAssociatorAtVertexAK5PF
  +process.impactParameterTagInfosAK5PF
  +process.secondaryVertexTagInfosAK5PF
  +process.simpleSecondaryVertexHighEffBJetTagsAK5PF
  +process.simpleSecondaryVertexHighPurBJetTagsAK5PF
  +process.combinedSecondaryVertexBJetTagsAK5PF
  )

##############################################################################
## Electron Extras
##############################################################################
process.icElectronR9Calculator = cms.EDProducer('ICElectronR9Calculator',
    input = cms.InputTag("selectedElectrons")
)

process.icElectronHcalDepthCalculator = cms.EDProducer('ICElectronHcalDepthCalculator',
    input = cms.InputTag("selectedElectrons")
)

process.icElectronConversionCalculator = cms.EDProducer('ICElectronConversionCalculator',
    input       = cms.InputTag("selectedElectrons"),
    beamspot    = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("allConversions")
)

process.icHttElecIsoCheck = cms.EDProducer('ICHttElecIsoCheck',
    input         = cms.InputTag("selectedElectrons"),
    pfChargedAll  = cms.InputTag("pfAllChargedParticles")
)

process.icHttMuonOverlapCheck = cms.EDProducer('ICHttMuonOverlapCheck',
    input = cms.InputTag("selectedElectrons"),
    muons = cms.InputTag("muons")
)

##############################################################################
# Electron Module
##############################################################################
process.icElectronProducer = cms.EDProducer('ICElectronProducer',
    branch                    = cms.string("electrons"),
    input                     = cms.InputTag("selectedElectrons"),
    includeR9                 = cms.InputTag("icElectronR9Calculator"),
    includeHcalSum            = cms.InputTag("icElectronHcalDepthCalculator"),
    includeConversionMatches  = cms.InputTag("icElectronConversionCalculator"),
    includeVertexIP           = cms.InputTag("selectedVertices"),
    includeBeamspotIP         = cms.InputTag("offlineBeamSpot"),
    includeFloats = cms.PSet(
      # mvaTrigV0       = cms.InputTag("mvaTrigV0"),
      # mvaNonTrigV0    = cms.InputTag("mvaNonTrigV0"),
      trackInIsoSum   = cms.InputTag("icHttElecIsoCheck"),
      matchedRecoMuon = cms.InputTag("icHttMuonOverlapCheck")
    ),
    pfIso03 = cms.PSet(
      chargedAll  = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
      charged     = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
      neutral     = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
      gamma       = cms.InputTag("elPFIsoValueGamma03PFIdPFIso"),
      pu          = cms.InputTag("elPFIsoValuePU03PFIdPFIso")
    ),
    pfIso04 = cms.PSet(
      chargedAll  = cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
      charged     = cms.InputTag("elPFIsoValueCharged04PFIdPFIso"),
      neutral     = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
      gamma       = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
      pu          = cms.InputTag("elPFIsoValuePU04PFIdPFIso")
    )
)

##############################################################################
# Reco Muon Module
##############################################################################
process.icMuonProducer = cms.EDProducer('ICMuonProducer',
    branch                    = cms.string("muons"),
    input                     = cms.InputTag("selectedMuons"),
    isPF                      = cms.bool(False),
    includeVertexIP           = cms.InputTag("selectedVertices"),
    includeBeamspotIP         = cms.InputTag("offlineBeamSpot"),
    includeFloats = cms.PSet(
    )
)

##############################################################################
# PF Muon Module
##############################################################################
process.icPFMuonProducer = cms.EDProducer('ICMuonProducer',
    branch                    = cms.string("pfMuons"),
    input                     = cms.InputTag("selectedPFMuons"),
    isPF                      = cms.bool(True),
    includeVertexIP           = cms.InputTag("selectedVertices"),
    includeBeamspotIP         = cms.InputTag("offlineBeamSpot"),
    includeFloats = cms.PSet(
    ),
    pfIso03 = cms.PSet(
      chargedAll  = cms.InputTag("muPFIsoValueChargedAll03PFIso"),
      charged     = cms.InputTag("muPFIsoValueCharged03PFIso"),
      neutral     = cms.InputTag("muPFIsoValueNeutral03PFIso"),
      gamma       = cms.InputTag("muPFIsoValueGamma03PFIso"),
      pu          = cms.InputTag("muPFIsoValuePU03PFIso")
    ),
    pfIso04 = cms.PSet(
      chargedAll  = cms.InputTag("muPFIsoValueChargedAll04PFIso"),
      charged     = cms.InputTag("muPFIsoValueCharged04PFIso"),
      neutral     = cms.InputTag("muPFIsoValueNeutral04PFIso"),
      gamma       = cms.InputTag("muPFIsoValueGamma04PFIso"),
      pu          = cms.InputTag("muPFIsoValuePU04PFIso")
    )
)

##############################################################################
# CaloJet Module
##############################################################################
process.icCaloJetProducer = cms.EDProducer('ICNewCaloJetProducer',
    branch                    = cms.string("caloJets"),
    input                     = cms.InputTag("selectedCaloJets"),
    includeJetFlavour         = cms.bool(True),
    inputJetFlavour           = cms.InputTag("icCaloJetFlavourCalculator"),
    applyJECs                 = cms.bool(True),
    includeJECs               = cms.bool(True),
    JECs                      = cms.PSet(
      L1FastJet  = cms.string("ak5CaloL1Fastjet"),
      L2Relative = cms.string("ak5CaloL2Relative"),
      L3Absolute = cms.string("ak5CaloL3Absolute")
    ),
    applyCutAfterJECs         = cms.bool(True),
    cutAfterJECs              = cms.string("pt > 30"),
    inputSVInfo               = cms.InputTag(""),
    requestSVInfo             = cms.bool(False),
    BTagDiscriminators        = cms.PSet(),
    specificConfig = cms.PSet(
      includeJetID    = cms.bool(True),
      inputJetID      = cms.InputTag("ak5JetID")
    )
)

##############################################################################
# PF Jet Module
##############################################################################
process.icPFJetProducer = cms.EDProducer('ICNewPFJetProducer',
    branch                    = cms.string("pfJets"),
    input                     = cms.InputTag("selectedPFJets"),
    includeJetFlavour         = cms.bool(True),
    inputJetFlavour           = cms.InputTag("icPFJetFlavourCalculator"),
    applyJECs                 = cms.bool(True),
    includeJECs               = cms.bool(True),
    JECs                      = cms.PSet(
      L1FastJet  = cms.string("ak5PFL1Fastjet"),
      L2Relative = cms.string("ak5PFL2Relative"),
      L3Absolute = cms.string("ak5PFL3Absolute")
    ),
    applyCutAfterJECs         = cms.bool(True),
    cutAfterJECs              = cms.string("pt > 30"),
    inputSVInfo               = cms.InputTag("secondaryVertexTagInfosAK5PF"),
    requestSVInfo             = cms.bool(True),
    BTagDiscriminators        = cms.PSet(
      simpleSecondaryVertexHighEff = cms.InputTag("simpleSecondaryVertexHighEffBJetTagsAK5PF"),
      simpleSecondaryVertexHighPur = cms.InputTag("simpleSecondaryVertexHighPurBJetTagsAK5PF"),
      combinedSecondaryVertex      = cms.InputTag("combinedSecondaryVertexBJetTagsAK5PF")
    ),
    specificConfig = cms.PSet(
      includePileupID    = cms.bool(False),
      inputPileupID      = cms.InputTag(""),
      includeTrackBasedVars = cms.bool(True),
      inputTracks           = cms.InputTag("generalTracks"),
      inputVertices         = cms.InputTag("offlinePrimaryVertices"),
      requestTracks         = cms.bool(True)
    )
)

##############################################################################
# Vertex Module
##############################################################################
process.icVertexProducer = cms.EDProducer('ICVertexProducer',
  branch  = cms.string("vertices"),
  input   = cms.InputTag("offlinePrimaryVertices"),
  firstVertexOnly = cms.bool(True),
  trackPtThreshold = cms.double(0.0),
  requestTracks = cms.bool(True)
)

process.icMergedTracks = cms.EDProducer('ICTrackMerger',
  merge = cms.VInputTag(
    cms.InputTag("icPFJetProducer", "requestedTracks"),
    cms.InputTag("icVertexProducer", "requestedTracks")
  )
)


##############################################################################
# Track Module
##############################################################################

process.icTrackProducer = cms.EDProducer('ICLightTrackProducer',
  branch  = cms.string("tracks"),
  input   = cms.InputTag("icMergedTracks")
)

process.icEventProducer = cms.EDProducer('ICEventProducer')


process.p = cms.Path(
  process.btaggingSequenceAK5PF+
  process.selectedVertices+
  process.selectedElectrons+
  process.selectedMuons+
  process.selectedCaloJets+
  process.selectedPFJets+
  process.selectedPFMuons+
  process.pfParticleSelectionSequence+
  process.eleIsoSequence+
  process.muoIsoSequence+
  process.icElectronR9Calculator+
  process.icElectronHcalDepthCalculator+
  process.icElectronConversionCalculator+
  process.icHttElecIsoCheck+
  process.icHttMuonOverlapCheck+
  process.jetPartons+
  process.caloJetPartonMatches+
  process.caloJetFlavourAssociation+
  process.icCaloJetFlavourCalculator+
  process.pfJetPartonMatches+
  process.pfJetFlavourAssociation+
  process.icPFJetFlavourCalculator+
  process.icElectronProducer+
  process.icMuonProducer+
  process.icPFMuonProducer+
  process.icCaloJetProducer+
  process.icPFJetProducer+
  process.icVertexProducer+
  process.icMergedTracks+
  process.icTrackProducer+
  process.icEventProducer
  )

#print process.dumpPython()