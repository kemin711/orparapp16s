#AUTOMAKE_OPTIONS = subdir-objects

#headers = dnaqualstore.h consensus.h boxpile.h piechart.h textpile.h stalactitegraph.h constaxmap.h blastnrow.h 
libheaders = dnaqualstore.h consensus.h blastnrow.h 
headers = stalactitegraph.h 
sources = dnaqualstore.cpp blastnrow.cpp consensus.cpp stalactitegraph.cpp 

pkginclude_HEADERS = $(libheaders)
lib_LTLIBRARIES = liborparaconsensus.la
liborparaconsensus_la_SOURCES = $(libheaders) $(headers) $(sources)
liborparaconsensus_la_CXXFLAGS = $(orpara_CXXFLAGS) $(plotter_CXXFLAGS)
liborparaconsensus_la_LDFLAGS = $(orpara_LDFLAGS) $(plotter_LDFLAGS) -version-info 0:0:0
liborparaconsensus_la_LIBADD = $(orpara_LIBS) $(plotter_LIBS)
