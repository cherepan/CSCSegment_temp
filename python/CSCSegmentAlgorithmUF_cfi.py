import FWCore.ParameterSet.Config as cms


UF_ME1234 = cms.PSet(

    minHitsPerSegment                 =  cms.int32(3),
    prePrun                           =  cms.bool(True),
    prePrunLimit                      =  cms.double(3.17),
    NormChi2Cut3D                     =  cms.double(10.0),
    NormChi2Cut2D                     =  cms.double(20.0),
    CorrectTheErrors                  =  cms.bool(True),
    CSCDebug                          =  cms.untracked.bool(False),


    recoverMissingWireHits            =  cms.bool(False),   # recent implementation creates wierd segments, set to False (To be reviewed)
    recoverMissingStripHits           =  cms.bool(False),   # recent implementation creates wierd segments, set to False (To be reviewed) 





    
#################################################
## consequence of awful oragnisation
## have to add here for debugging purposes
## it does not belong to Segment Algorithm
## at all!!! 
    #    a XT asymmetry model parameter
    XTasymmetry_ME1a = cms.double(0.023),
    XTasymmetry_ME1b = cms.double(0.01),
    XTasymmetry_ME12 = cms.double(0.015),
    XTasymmetry_ME13 = cms.double(0.02),
    XTasymmetry_ME21 = cms.double(0.023),
    XTasymmetry_ME22 = cms.double(0.023),
    XTasymmetry_ME31 = cms.double(0.023),
    XTasymmetry_ME32 = cms.double(0.023),
    XTasymmetry_ME41 = cms.double(0.023),
    #
    #    constant systematics (in cm)
    ConstSyst_ME1a = cms.double(0.01),
    ConstSyst_ME1b = cms.double(0.02),
    ConstSyst_ME12 = cms.double(0.02),
    ConstSyst_ME13 = cms.double(0.03),
    ConstSyst_ME21 = cms.double(0.03),
    ConstSyst_ME22 = cms.double(0.03),
    ConstSyst_ME31 = cms.double(0.03),
    ConstSyst_ME32 = cms.double(0.03),
    ConstSyst_ME41 = cms.double(0.03),
    #
    #    3 time bins noise (in ADC counts)
    NoiseLevel_ME1a = cms.double(9.0),
    NoiseLevel_ME1b = cms.double(6.0),
    NoiseLevel_ME12 = cms.double(7.0),
    NoiseLevel_ME13 = cms.double(4.0),
    NoiseLevel_ME21 = cms.double(5.0),
    NoiseLevel_ME22 = cms.double(7.0),
    NoiseLevel_ME31 = cms.double(5.0),
    NoiseLevel_ME32 = cms.double(7.0),
    NoiseLevel_ME41 = cms.double(5.0),

    SeedBig = cms.double(0.001500),
    SeedSmall = cms.double(0.000200)
#######
#######
    
)

UF_ME1A = cms.PSet(
    minHitsPerSegment = cms.int32(3),
    prePrun = cms.bool(True),
    prePrunLimit = cms.double(3.17),
    NormChi2Cut3D = cms.double(10.0),
    NormChi2Cut2D = cms.double(20.0),
    CorrectTheErrors = cms.bool(True),
    CSCDebug = cms.untracked.bool(False),
    onlyBestSegment = cms.bool(False),

    recoverMissingWireHits            =  cms.bool(False),   # recent implementation creates wierd segments, set to False (To be reviewed)
    recoverMissingStripHits           =  cms.bool(False),   # recent implementation creates wierd segments, set to False (To be reviewed) 

    
#################################################
## consequence of awful oragnisation
## have to add here for debugging purposes
## it does not belong to Segment Algorithm
## at all!!! 

    XTasymmetry_ME1a = cms.double(0.023),
    XTasymmetry_ME1b = cms.double(0.01),
    XTasymmetry_ME12 = cms.double(0.015),
    XTasymmetry_ME13 = cms.double(0.02),
    XTasymmetry_ME21 = cms.double(0.023),
    XTasymmetry_ME22 = cms.double(0.023),
    XTasymmetry_ME31 = cms.double(0.023),
    XTasymmetry_ME32 = cms.double(0.023),
    XTasymmetry_ME41 = cms.double(0.023),
    #
    #    constant systematics (in cm)
    ConstSyst_ME1a = cms.double(0.01),
    ConstSyst_ME1b = cms.double(0.02),
    ConstSyst_ME12 = cms.double(0.02),
    ConstSyst_ME13 = cms.double(0.03),
    ConstSyst_ME21 = cms.double(0.03),
    ConstSyst_ME22 = cms.double(0.03),
    ConstSyst_ME31 = cms.double(0.03),
    ConstSyst_ME32 = cms.double(0.03),
    ConstSyst_ME41 = cms.double(0.03),
    #
    #    3 time bins noise (in ADC counts)
    NoiseLevel_ME1a = cms.double(9.0),
    NoiseLevel_ME1b = cms.double(6.0),
    NoiseLevel_ME12 = cms.double(7.0),
    NoiseLevel_ME13 = cms.double(4.0),
    NoiseLevel_ME21 = cms.double(5.0),
    NoiseLevel_ME22 = cms.double(7.0),
    NoiseLevel_ME31 = cms.double(5.0),
    NoiseLevel_ME32 = cms.double(7.0),
    NoiseLevel_ME41 = cms.double(5.0),
    
    SeedBig = cms.double(0.001500),
    SeedSmall = cms.double(0.000200)


#######
#######
)


CSCSegAlgoUF = cms.PSet(
    algo_name = cms.string('CSCSegAlgoUF'),
    algo_psets = cms.VPSet( cms.PSet(UF_ME1234), cms.PSet(UF_ME1A) ),
    chamber_types = cms.vstring('ME1/a','ME1/b','ME1/2','ME1/3','ME2/1','ME2/2','ME3/1','ME3/2','ME4/1','ME4/2'),
    parameters_per_chamber_type = cms.vint32(2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
)

