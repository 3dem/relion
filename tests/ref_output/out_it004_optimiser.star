# RELION optimiser
# --o test_output/out --i ap_s_c2d_red.star --particle_diameter 200 --angpix 1.0 --ref ref_model.mrc --ini_high 50 --ctf --ctf_corrected_ref --iter 5 --tau2_fudge 2 --K 1 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --j 1 --memory_per_thread 4 --random_seed 1993 --onthefly_shifts --scale 

data_optimiser_general

_rlnOutputRootName                                    test_output/out
_rlnModelStarFile                                     test_output/out_it004_model.star
_rlnExperimentalDataStarFile                          test_output/out_it004_data.star
_rlnOrientSamplingStarFile                            test_output/out_it004_sampling.star
_rlnCurrentIteration                                             4
_rlnNumberOfIterations                                           5
_rlnDoSplitRandomHalves                                          0
_rlnJoinHalvesUntilThisResolution                        -1.000000
_rlnAdaptiveOversampleOrder                                      1
_rlnAdaptiveOversampleFraction                            0.999000
_rlnRandomSeed                                                1993
_rlnParticleDiameter                                    200.000000
_rlnWidthMaskEdge                                                5
_rlnDoZeroMask                                                   1
_rlnDoSolventFlattening                                          1
_rlnSolventMaskName                                   None
_rlnSolventMask2Name                                  None
_rlnTauSpectrumName                                   None
_rlnCoarseImageSize                                             30
_rlnMaximumCoarseImageSize                                     100
_rlnHighresLimitExpectation                              -1.000000
_rlnIncrementImageSize                                          10
_rlnDoMapEstimation                                              1
_rlnDoAutoRefine                                                 0
_rlnAutoLocalSearchesHealpixOrder                                4
_rlnNumberOfIterWithoutResolutionGain                            4
_rlnBestResolutionThusFar                                 0.019774
_rlnNumberOfIterWithoutChangingAssignments                       3
_rlnDoSkipAlign                                                  0
_rlnDoSkipRotate                                                 0
_rlnOverallAccuracyRotations                             25.655556
_rlnOverallAccuracyTranslations                           2.711111
_rlnChangesOptimalOrientations                        2.521697e-07
_rlnChangesOptimalOffsets                                 0.000000
_rlnChangesOptimalClasses                                 0.000000
_rlnSmallestChangesOrientations                       1.581062e-07
_rlnSmallestChangesOffsets                                0.000000
_rlnSmallestChangesClasses                                       0
_rlnHasConverged                                                 0
_rlnHasHighFscAtResolLimit                                       0
_rlnHasLargeSizeIncreaseIterationsAgo                            0
_rlnDoCorrectNorm                                                1
_rlnDoCorrectScale                                               0
_rlnDoCorrectCtf                                                 1
_rlnDoRealignMovies                                              0
_rlnDoIgnoreCtfUntilFirstPeak                                    0
_rlnCtfDataArePhaseFlipped                                       0
_rlnDoOnlyFlipCtfPhases                                          0
_rlnRefsAreCtfCorrected                                          1
_rlnFixSigmaNoiseEstimates                                       0
_rlnFixSigmaOffsetEstimates                                      0
_rlnMaxNumberOfPooledParticles                                   1
_rlnAvailableMemory                                       4.000000
 
