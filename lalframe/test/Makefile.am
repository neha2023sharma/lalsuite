EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

SUBDIRS = \
	$(END_OF_LIST)

# Add compiled test programs to this variable
test_programs += LALFrSeriesTest

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers += MakeFrames

FRAMEFILES = \
	F-TEST-600000000-60.gwf \
	F-TEST-600000060-60.gwf \
	F-TEST-600000120-60.gwf \
	H-CAL_FAC_V03-729273600-5094000.gwf \
	H-CAL_REF_V03-734073939-64.gwf \
	L-CAL_FAC_V03-729273600-5094000.gwf \
	L-CAL_REF_V03-731488397-64.gwf \
	$(END_OF_LIST)

BUILT_SOURCES = \
	catalog \
	catalog.test \
	$(END_OF_LIST)

catalog:
	$(AM_V_GEN)$(SHELL) $(srcdir)/mkcatalog.sh $@ $(srcdir) $(pkgdatadir)

catalog.test:
	$(AM_V_GEN)$(SHELL) $(srcdir)/mkcatalog.sh $@ $(srcdir) $(srcdir); \
	echo "X EMPTY 0 0 file:///dev/null" >> $@

MOSTLYCLEANFILES = \
	*.[0-9][0-9][0-9] \
	*.out \
	H-H1_LSC_AS_Q-600000120-60.gwf \
	Response*.txt \
	catalog \
	catalog.out \
	catalog.test \
	$(END_OF_LIST)

EXTRA_DIST += \
	$(FRAMEFILES) \
	H1-CAL-V03-729273600-734367600.cache \
	L1-CAL-V03-729273600-734367600.cache \
	mkcatalog.sh \
	$(END_OF_LIST)

if HAVE_PYTHON
SUBDIRS += python
endif
