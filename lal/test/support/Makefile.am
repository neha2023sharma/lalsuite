EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += ConfigFileTest
test_programs += H5FileIOTest
test_programs += LALMath3DPlotTest
test_programs += LALMathNDPlotTest
test_programs += PrintFTSeriesTest
test_programs += PrintVectorTest
test_programs += ReadFTSeriesTest
test_programs += StreamInputTest
test_programs += StreamSeriesInputTest
test_programs += TranslateAnglesTest
test_programs += TranslateMJDTest
test_programs += UserInputTest
test_programs += UserInputParseTest

# Add shell, Python, etc. test scripts to this variable
test_scripts += GzipTest.sh

# Add any helper programs required by tests to this variable
test_helpers += GzipTest

MOSTLYCLEANFILES = \
	*.dat \
	*.out \
	*PrintVector.00* \
	test.h5 \
	ConfigFile.cfg \
	Math3DNotebook.nb \
	MathNDNotebook.nb \
	$(END_OF_LIST)

EXTRA_DIST += \
	ConfigFileSample.cfg \
	ConfigFileSample2.cfg \
	ConfigFileSample3.cfg.gz \
	StreamInput.data \
	StreamSeriesInput.data \
	hdf5_attr_ascii.h5 \
	hdf5_attr_utf8.h5 \
	$(END_OF_LIST)
