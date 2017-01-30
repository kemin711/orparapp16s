#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([orparapp16s], [0.1.0], [kmzhou4@yahoo.com])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([src/consensus/blastnrow.h])
AC_CONFIG_HEADERS([config.h])
AM_PROG_AR
LT_PREREQ([2.4])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL

AC_LANG([C], [C++])
# Checks for libraries.
AC_LANG_PUSH([C++])
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])
AC_LANG_POP([C++])

if test "$prefix" = "NONE"; then
	prefix=$ac_default_prefix
fi
# set up gtest
save_LDFLAGS=$LDFLAGS
save_LIBS=$LIBS
LDFLAGS="-L/usr/lib/x86_64-linux-gnu $LDFLAGS"
LIBS="-lpthread $LIBS"
#this macro is not setting up PTREAD_XXX
AX_PTHREAD([
  AC_DEFINE([HAVE_PTHREAD], [1], [Define if you have POSIX threads libraries and headers])
  gtest_LIBS="$PTHREAD_LIBS $LIBS"
  gtest_CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
  gtest_CXXFLAGS="$CXXFLAGS $PTHREAD_CXXFLAGS"
  gtest_CC="$PTHREAD_CC"
  gtest_CXX="$PTHREAD_CXX"])
LDFLAGS=$save_LDFLAGS
LIBS=$save_LIBS

AC_ARG_VAR([GTESTINC], [gtest PATH for gtest include dir])
AC_ARG_VAR([GTESTLIB], [gtest PATH for gtest lib dir])
if (test -n "$GTESTINC" && test -z "$GTESTLIB") || \
	(test -z "$GTESTINC" && test -n "$GTESTLIB"); then 
   AC_MSG_ERROR([please specify both lib and include dir for gtest])
fi
AC_LANG_PUSH([C++])
save_LIBS=$LIBS
save_LDFLAGS=$LDFLAGS
save_CXXFLAGS=$CXXFLAGS
LDFLAGS="-L$prefix/lib $LDFLAGS"
CXXFLAGS="-I$prefix/include $CXXFLAGS"
LIBS="$LIBS -lgtest -lgtest_main -lgmock -lgmock_main -lpthread"
if test -n "$GTESTINC"; then
	LDFLAGS="-L$GTESTLIB $LDFLAGS"
	CXXFLAGS="-I$GTESTINC $CXXFLAGS"
fi
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([#include <gmock/gmock.h>], [testing::Cardinality dummy])],
	[gtest_LIBS=$LIBS] [gtest_CXXFLAGS=$CXXFLAGS] [gtest_LDFLAGS=$LDFLAGS] [HAVE_GTEST=1] [HAVE_GMOCK=1],
	[AC_MSG_WARN([libgtest and libgmock not found, will not be able to make check.])])
AC_SUBST([gtest_CFLAGS])
AC_SUBST([gtest_CXXFLAGS])
AC_SUBST([gtest_LIBS])
LIBS=$save_LIBS
LDFLAGS=$save_LDFLAGS
CXXFLAGS=$save_CXXFLAGS
#AC_LANG_POP([C++])
### end of gtest

#AC_LANG_PUSH([C++])
# bamtools lib check, this will check both header and lib
save_LIBS=$LIBS
save_CXXFLAGS=$CXXFLAGS
save_LDFLAGS=$LDFLAGS
bam_LIBS="-lbamtools"
bam_CXXFLAGS="-I$prefix/include/bamtools"
bam_LDFLAGS="-L$prefix/lib/bamtools"
if test "$prefix" != "/usr/local"; then
	bam_CXXFLAGS="$bam_CXXFLAGS -I/usr/local/include/bamtools"
	bam_LDFLAGS="$bam_LDFLAGS -L/usr/local/lib/bamtools"
fi

AC_ARG_VAR([BAMTOOLSINC], [path for bamtools include dir, canb environment variable or confg variable])
AC_ARG_VAR([BAMTOOLSLIB], [path for bamtools lib dir])
if test -n "$BAMTOOLSINC"; then 
	bam_CXXFLAGS=$BAMTOOLSINC
fi
if test -n "$BAMTOOLSLIB"; then 
	bam_LDFLAGS=$BAMTOOLSLIB
fi
CXXFLAGS="$bam_CXXFLAGS $CXXFLAGS"
LDFLAGS="$bam_LDFLAGS $LDFLAGS"
LIBS="$LIBS $bam_LIBS"
AC_CHECK_HEADERS([api/BamReader.h], [bamtools_found=yes], [AC_MSG_ERROR([Unable to find bambools header])])
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([#include <api/BamReader.h>],  [BamTools::BamReader dummy])],
		[HAVE_BAMTOOLS=1],
		[AC_MSG_ERROR([libbamtools is not installed])])
AC_SUBST([bam_LIBS])
AC_SUBST([bam_CXXFLAGS])
AC_SUBST([bam_LDFLAGS])
LIBS=$save_LIBS
CXXFLAGS=$save_CXXFLAGS
LDFLAGS=$save_LDFLAGS

# boost lib check
# each program need to add its own library. there too many
# libraries. Many of boost libs are header only
save_LIBS=$LIBS
save_CXXFLAGS=$CXXFLAGS
save_LDFLAGS=$LDFLAGS
boost_LIBS="-lboost_program_options"
boost_CXXFLAGS="-I$prefix/include/boost"
boost_LDFLAGS="-L$prefix/lib"
if test "$prefix" != "/usr/local"; then
	boost_CXXFLAGS="$boost_CXXFLAGS -I/usr/local/include/boost"
	boost_LDFLAGS="$boost_LDFLAGS -L/usr/local/lib"
fi
AC_ARG_VAR([BOOSTINC], [PATH for boost include dir, canb environment variable or confg variable])
AC_ARG_VAR([BOOSTLIB], [PATH for boost lib dir])
if (test -n "$BOOSTINC"); then 
	boost_CXXFLAGS="-I$BOOSTINC"
fi
if (test -n "$BOOSTLIB"); then 
	boost_LDFLAGS="-L$BOOSTLIB"
fi
CXXFLAGS="$boost_CXXFLAGS $CXXFLAGS"
LDFLAGS="$boost_LDFLAGS $LDFLAGS"
LIBS="$LIBS $boost_LIBS"
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([#include <program_options/parsers.hpp>],  [boost::program_options::option dummy])],
		[HAVE_BOOST=1],
		[AC_MSG_ERROR([boost is not installed])])
AC_SUBST([boost_LIBS])
AC_SUBST([boost_CXXFLAGS])
AC_SUBST([boost_LDFLAGS])
LIBS=$save_LIBS
CXXFLAGS=$save_CXXFLAGS
LDFLAGS=$save_LDFLAGS
if test "$HAVE_BOOST"; then
	AC_MSG_NOTICE([found boost])
fi

# Checks for libraries. `-lkemin':
save_LIBS=$LIBS
save_CXXFLAGS=$CXXFLAGS
save_LDFLAGS=$LDFLAGS
orpara_LIBS="-lkemin16s"
orpara_CXXFLAGS="-I$prefix/include/orpara16s"
orpara_LDFLAGS="-L$prefix/lib"
AC_ARG_VAR([ORPARAIND], [PATH for orpara16s include dir, can be environment variable or confg variable])
AC_ARG_VAR([ORPARALIB], [PATH for orpara16s lib dir])
if (test -n "$ORPARAINC"); then 
	orpara_CXXFLAGS="-I$ORPARAINC"
fi
if (test -n "$ORPARALIB"); then 
	orpara_LDFLAGS="-L$ORPARALIBS"
fi
LIBS="$LIBS $orpara_LIB"
CXXFLAGS="$orpara_CXXFLAGS $CXXFLAGS"
LDFLAGS="$orpara_LDFLAGS $LDFLAGS"
AC_CHECK_HEADERS([bioseq.h], [orpara_found=yes], [AC_MSG_ERROR([Unable to find orpara header])])
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([#include <bioseq.h>],  [orpara::bioseq dummy])],
		[HAVE_ORPARA=1],
		[AC_MSG_ERROR([orpara16s is not installed])])
AC_SUBST([orpara_LIBS])
AC_SUBST([orpara_CXXFLAGS])
AC_SUBST([orpara_LDFLAGS])
LIBS=$save_LIBS
CXXFLAGS=$save_CXXFLAGS
LDFLAGS=$save_LDFLAGS

# Checks for libraries. libplotter 
save_LIBS=$LIBS
save_CXXFLAGS=$CXXFLAGS
save_LDFLAGS=$LDFLAGS
# the x** libs are system level libs that the build host should have
plotter_LIBS="-lplotter -lXaw -lXmu -lXt -lSM -lICE -lXext -lX11 -lpng -lz -lm"
plotter_CXXFLAGS="-I$prefix/include"
plotter_LDFLAGS="-L$prefix/lib"
LIBS="$LIBS $plotter_LIBS"
AC_ARG_VAR([PLOTUTILINC], [PATH for plotutil include dir, can be environment variable or confg variable])
AC_ARG_VAR([PLOTUTILLIB], [PATH for plotutil lib dir])
if (test -n "$PLOTUTILINC"); then 
	plotter_CXXFLAGS="-I$PLOTUTILINC"
fi
if (test -n "$PLOTUTILLIB"); then 
	plotter_LDFLAGS="-L$PLOTUTILLIB"
fi
CXXFLAGS="$plotter_CXXFLAGS $CXXFLAGS"
LDFLAGS="$plotter_LDFLAGS $LDFLAGS"
AC_LINK_IFELSE(
		[AC_LANG_PROGRAM([#include <plotter.h>],  [SVGPlotter dummy])],
		[HAVE_PLOTUTIL=1],
		[AC_MSG_ERROR([plotutil is not installed])])
AC_SUBST([plotter_LIBS])
AC_SUBST([plotter_CXXFLAGS])
AC_SUBST([plotter_LDFLAGS])
LIBS=$save_LIBS
CXXFLAGS=$save_CXXFLAGS
LDFLAGS=$save_LDFLAGS
AC_LANG_POP([C++])

# Checks for header files.
AC_CHECK_HEADERS([string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT8_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor mkdir pow strchr strstr])

AC_CONFIG_FILES([Makefile
                 src/Makefile
					  src/consensus/Makefile
                 src/seqconv/Makefile
                 src/test/Makefile])
AC_OUTPUT
