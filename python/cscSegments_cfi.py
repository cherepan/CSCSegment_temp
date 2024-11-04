import FWCore.ParameterSet.Config as cms

from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmSK_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmTC_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmDF_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmST_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmRU_cfi import *
from RecoLocalMuon.CSCSegment.CSCSegmentAlgorithmUF_cfi import *

cscSegments = cms.EDProducer("CSCSegmentProducer",
    #
    #    Parameters for strip hits
    #

    CSCStripPeakThreshold = cms.double(10.0),
    CSCStripClusterChargeCut = cms.double(25.0),
    CSCStripxtalksOffset = cms.double(0.03),
    UseAverageTime = cms.bool(False),
    UseParabolaFit = cms.bool(False),
    UseFivePoleFit = cms.bool(True),
    CSCWireClusterDeltaT = cms.int32(1),
    CSCUseReducedWireTimeWindow = cms.bool(False),
    CSCWireTimeWindowLow = cms.int32(0),
    CSCWireTimeWindowHigh = cms.int32(15),
    CSCUseCalibrations = cms.bool(True),

    CSCUseStaticPedestals = cms.bool(False),
    CSCNoOfTimeBinsForDynamicPedestal = cms.int32(2),
#    wireDigiTag = cms.InputTag("muonCSCDigis","MuonCSCWireDigi"),
#    stripDigiTag = cms.InputTag("muonCSCDigis","MuonCSCStripDigi"),
    readBadChannels = cms.bool(True),
    readBadChambers = cms.bool(True),
    CSCUseTimingCorrections = cms.bool(True),
    CSCUseGasGainCorrections = cms.bool(True),
    CSCDebug = cms.untracked.bool(False),
    CSCstripWireDeltaTime = cms.int32(8),
    CSCStripClusterSize = cms.untracked.int32(3),

    # Define input
    inputObjects = cms.InputTag("csc2DRecHits"),
    # Choice of the building algo: 1 SK, 2 TC, 3 DF, 4 ST, 5 RU, 6 UF ...
    algo_type = cms.int32(6),
    # std::vector<edm::ParameterSet>
    algo_psets = cms.VPSet(
        cms.PSet(
            CSCSegAlgoSK
        ), 
        cms.PSet(
            CSCSegAlgoTC
        ), 
        cms.PSet(
            CSCSegAlgoDF
        ), 
        cms.PSet(
            CSCSegAlgoST
        ),
        cms.PSet(
            CSCSegAlgoRU
        ),
        cms.PSet(
            CSCSegAlgoUF
        )

     )
)


