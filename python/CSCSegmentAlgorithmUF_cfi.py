import FWCore.ParameterSet.Config as cms


UF_ME1234 = cms.PSet(

    minHitsPerSegment = cms.int32(3),
    prePrun = cms.bool(True),
    prePrunLimit = cms.double(3.17),
    NormChi2Cut3D = cms.double(10.0),
    NormChi2Cut2D = cms.double(20.0),
    CorrectTheErrors = cms.bool(True),
    CSCDebug = cms.untracked.bool(False),
    onlyBestSegment = cms.bool(False)


)

UF_ME1A = cms.PSet(
    minHitsPerSegment = cms.int32(3),
    prePrun = cms.bool(True),
    prePrunLimit = cms.double(3.17),
    NormChi2Cut3D = cms.double(10.0),
    NormChi2Cut2D = cms.double(20.0),
    CorrectTheErrors = cms.bool(True),
    CSCDebug = cms.untracked.bool(False),
    onlyBestSegment = cms.bool(False)


)


CSCSegAlgoUF = cms.PSet(
    algo_name = cms.string('CSCSegAlgoST'),
    algo_psets = cms.VPSet( cms.PSet(UF_ME1234), cms.PSet(UF_ME1A) ),
    chamber_types = cms.vstring('ME1/a','ME1/b','ME1/2','ME1/3','ME2/1','ME2/2','ME3/1','ME3/2','ME4/1','ME4/2'),
    parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

