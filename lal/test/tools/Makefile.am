EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += ComputeTransferTest
test_programs += DetResponseTest
test_programs += DetectorSiteTest
test_programs += FrequencySeriesTest
test_programs += SegmentsTest
test_programs += SequenceTest
test_programs += SkymapTest
test_programs += TimeSeriesInterpTest
test_programs += TimeSeriesTest
test_programs += TriggerInterpolateTest
test_programs += UnitsTest
test_programs += ValueTest
#test_programs += CoherentEstimationTest

# Add shell, Python, etc. test scripts to this variable
test_scripts += IndependentDetResponseTest.sh

# Add any helper programs required by tests to this variable
test_helpers += IndependentDetResponseTest

MOSTLYCLEANFILES = \
	PrintVector.* \
	circ_series.txt \
	cross.txt \
	cross_at_0_0.txt \
	cross_at_2_10.txt \
	cross_at_4_15.txt \
	cross_at_m4_15.txt \
	cross_series.txt \
	cross_sq_avg.txt \
	cross_sq_time_avg.txt \
	ff_diff_cros.txt \
	ff_diff_plus.txt \
	ff_lal_cros.txt \
	ff_lal_plus.txt \
	ff_local_cros.txt \
	ff_local_plus.txt \
	ff_rms_diff_cros_vs_fudge.txt \
	ff_rms_diff_plus_vs_fudge.txt \
	phi.txt \
	plus.txt \
	plus_at_0_0.txt \
	plus_at_2_10.txt \
	plus_at_4_15.txt \
	plus_at_m4_15.txt \
	plus_series.txt \
	plus_sq_avg.txt \
	plus_sq_time_avg.txt \
	scalar_series.txt \
	sum_of_sq_time_avg.txt \
	sum_series.txt \
	sum_sq.txt \
	sum_sq_avg.txt \
	theta.txt \
	transfer.out \
	$(END_OF_LIST)

EXTRA_DIST += \
	ComputeCalibrationFactorsTest.c \
	indTestDetResGEO.cfg \
	indTestDetResLHO.cfg \
	indTestDetResLLO.cfg \
	resamp_figs1.png \
	resamp_figs2.png \
	resamp_figs3.png \
	$(END_OF_LIST)
