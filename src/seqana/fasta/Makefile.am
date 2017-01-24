AM_CXXFLAGS=$(orpara_CXXFLAGS)
AM_LDFLAGS=$(orpara_LDFLAGS)
LDADD=$(orpara_LIBS)

bin_PROGRAMS = fastabund 

fastabund_SOURCES = fastabund.cpp 
fastabund_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/consensus $(AM_CXXFLAGS) $(plotter_CXXFLAGS)
fastabund_LDFLAGS = -L$(top_srcdir)/src $(AM_LDFLAGS) $(plotter_LDFLALGS)
fastabund_LDADD = $(top_srcdir)/src/consensus/liborparaconsensus.la $(LDADD) $(plotter_LIBS) -lpthread



