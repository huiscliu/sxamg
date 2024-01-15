
AC_PREREQ([2.71])

AC_INIT
m4_ifdef([AC_DISABLE_OPTION_CHECKING], [AC_DISABLE_OPTION_CHECKING()], [])

dnl This MUST precede any other macro
AC_CONFIG_AUX_DIR([utils])

AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_FILES([Makefile Makefile.inc])
AMG_DIR="`pwd`"

AC_MSG_NOTICE(============================== compilers ===========================)

#--------------------------------------- config files and compilers
# Checks for programs.
AC_PROG_CC(gcc clang icc)
AC_PROG_CPP

if test "$GCC" = "yes"; then
    if $CC -V 2>&1 | grep -q Intel; then
        CFLAGS="$CFLAGS -Wall"
    elif $CC -v 2>&1 | grep -q clang; then
        CFLAGS="$CFLAGS -Wall -Wextra -Warray-bounds"
    elif $CC -v 2>&1 | grep -q gcc; then
        CFLAGS="$CFLAGS -Wall -Wextra"

        for d in -Wmissing-declarations -Wmissing-prototypes -Wimplicit-function-declaration        \
            -Wunused -Wuninitialized -Wvla -Wformat -Wformat-security -Wmain -Wsequence-point       \
            -Wshadow -Wbad-function-cast -Wcomments -Wdeclaration-after-statement                   \
            -Waddress -Wpointer-sign -Wempty-body -Wnested-externs -Wold-style-declaration          \
            -Wold-style-definition -Wmissing-field-initializers -Wmaybe-uninitialized               \
            -Wno-return-local-addr -Wreturn-type -Wsequence-point -Wundef -Wunknown-pragmas         \
            -Wunused-label; do
            AC_MSG_CHECKING([whether $CC supports $d ])
            AC_LANG(C)
            CFLAGS_bak="$CFLAGS"
            CFLAGS="$CFLAGS $d"
            AC_COMPILE_IFELSE([AC_LANG_SOURCE([main(){}])], tmp=yes, tmp=no)
            if test $tmp = yes; then
                AC_MSG_RESULT(yes)
            else
                AC_MSG_RESULT(no)
                CFLAGS="$CFLAGS_bak"
            fi
        done
    fi

fi

dnl Change relative path in compilers to full path
AC_FULLPATH(CC)
AC_FULLPATH(CPP)

# The default linker
AC_SUBST(LINKER, "$LINKER")

AC_PREFIX_DEFAULT(/usr/local/sxamg)

AC_CHECK_PROGS(AR,[ar], no)
if test $AR = "no" ; then
    AC_MSG_ERROR(* AR is nott found!);
fi

AC_SUBST(AR)
AC_PROG_RANLIB

#--------------------------------------- shared library

AC_SUBST(LIB_SUFFIX, "")
AC_SUBST(BUILD_SHARED, "")
AC_SUBST(BUILD_SHARED_LIBS, "")
AC_SUBST(BUILD_CFLAGS, "")
AC_SUBST(BUILD_CXXFLAGS, "")
dnl AC_SUBST(BUILD_FFLAGS, "")
AC_ARG_ENABLE(shared,
[  --enable-shared         build shared library libsxamg.so (default)
  --disable-shared        build static library libsxamg.a
  --with-build-shared     the BUILD_SHARED macro (default "gcc -shared" or "${CC} -shared")],,enable_shared="auto")

if test "$enable_shared" = "auto"; then
dnl    if test "$target_vendor" = "apple"; then
dnl        # use static link by default on MacOS, before the rpath issue is solved
dnl        enable_shared="no"
dnl    else
        enable_shared="yes"
dnl    fi
fi

if test "$enable_shared" = "yes"; then
    /bin/rm -f src/libsxamg.a
    LIB_SUFFIX=".so"

    if test -n "$with_build_shared"; then
        BUILD_SHARED="$with_build_shared"
    elif gcc -v >/dev/null 2>&1 && \
         test "$cross_compiling" != "yes" && \
         test "$target_vendor" != "apple" # MacOS: force using "${CC} -shared"
    then
        BUILD_SHARED="gcc -shared"
    else
        BUILD_SHARED="\${CC} -shared"
    fi

    # TODO: figure out a proper way to determine when to use the PIC flag
    if test "$target_cpu" != "i686" -a "$target_cpu" != "i586" -a \
            "$target_cpu" != "i486" -a "$target_cpu" != "i386"; then
        BUILD_CFLAGS="-fPIC"
        BUILD_CXXFLAGS="-fPIC"
        dnl BUILD_FFLAGS="-fPIC"
    fi
    dnl TODO: check that the above settings work!
    AC_MSG_NOTICE([*** building shared library libsxamg$LIB_SUFFIX ("$BUILD_SHARED")])
else
    /bin/rm -f src/libsxamg.so
    LIB_SUFFIX=".a"
    BUILD_SHARED=""
    BUILD_CFLAGS=""
    BUILD_CXXFLAGS=""
    dnl BUILD_FFLAGS=""
    AC_MSG_NOTICE([*** building static library libsxamg$LIB_SUFFIX])
fi

#--------------------------------------- rpath
AC_SUBST(RPATH_FLAG, "")
AC_ARG_ENABLE(rpath,
[  --enable-rpath          enable use of rpath (default)
  --disable-rpath         disable use of rpath
  --with-rpath-flag=FLAG  compiler flag for rpath (e.g., "-Wl,-rpath,")],,
enable_rpath=yes)
if test -z "$with_rpath_flag" -a "$enable_rpath" != "no"; then
    AC_MSG_CHECKING([for compiler flag for specifying rpath ])
    LIBS_bak="$LIBS"
    for with_rpath_flag in "-Wl,-rpath," "-Wl,-rpath " "-Wl,--rpath " \
                           "-rpath " "--rpath " "-Wl,-R" "-R" "-R " ; do
        LIBS="$LIBS_bak ${with_rpath_flag}/usr/lib"
        AC_LINK_IFELSE([AC_LANG_PROGRAM()], enable_rpath=yes, enable_rpath=no)
        if test "$enable_rpath" = "yes"; then break; fi
    done
    LIBS="$LIBS_bak"
    if test "$enable_rpath" = "yes"; then
        AC_MSG_RESULT("$with_rpath_flag")
        RPATH_FLAG="$with_rpath_flag"
    else
        AC_MSG_RESULT(unknown)
        RPATH_FLAG="unknown"
    fi
fi

#--------------------------------------- header files
AC_CHECK_HEADERS_ONCE([sys/time.h])
AC_CHECK_HEADERS([sys/time.h sys/resource.h unistd.h])

#--------------------------------------- library functions.
AC_CHECK_FUNCS([gettimeofday getrusage])

# fpe needs -lm
LIBS="$LIBS -lm"

#--------------------------------------- basic flags
# Note: -I`pwd`/include and -L`pwd`/src must be the first term in CPPFLAGS
# and LDFLAGS respectively, which will be changed to AMG_INSTDIR/include by
# 'make install'
CPPFLAGS="-I${AMG_DIR}/include $CPPFLAGS"
LDFLAGS="-L${AMG_DIR}/src $LDFLAGS"

#--------------------------------------- Integer
AC_LANG(C)
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([void *])

AC_SUBST(USE_LONG_INT, "")
AC_ARG_ENABLE(big-int,
[  --enable-big-int        use long int for INT
  --disable-big-int       use int for INT (default),
  --with-int=type         integer type(long|long long)],, enable_big_int=no)
if test "$enable_big_int" = "yes"; then
    if test "$with_int" = "long long"; then
        AC_DEFINE(USE_LONG_LONG, 1)
        AC_DEFINE(USE_LONG, 0)
        AC_MSG_NOTICE([* use 'long long' for INT])
        AC_DEFINE(dFMT, ["lld"])
    elif test "$with_int" = "long"; then
        AC_DEFINE(USE_LONG_LONG, 0)
        AC_DEFINE(USE_LONG, 1)
        AC_MSG_NOTICE([* use 'long' for INT])
        AC_DEFINE(dFMT, ["ld"])
    else
        AC_MSG_NOTICE([* use 'long' for INT])
        AC_DEFINE(USE_LONG_LONG, 0)
        AC_DEFINE(USE_LONG, 1)
        AC_MSG_NOTICE([* use 'long' for INT])
        AC_DEFINE(dFMT, ["ld"])
    fi
else
    AC_MSG_NOTICE(* use 'int' for 'INT')
    AC_DEFINE(USE_LONG_LONG, 0)
    AC_DEFINE(USE_LONG, 0)
    AC_DEFINE(dFMT, ["d"])
fi

#--------------------------------------- Floating point
AC_LANG(C)
AC_SUBST(USE_LONG_DOUBLE, "0")
AC_TYPE_LONG_DOUBLE

AC_ARG_ENABLE(long-double,
[  --enable-long-double    use long double for FLOAT
  --disable-long-double   use double for FLOAT (default)],, enable_long_double=no)
if test "$enable_long_double" = "yes" -a "$ac_cv_type_long_double" = "yes"; then
    AC_CHECK_SIZEOF([long double])
    AC_DEFINE(USE_LONG_DOUBLE, 1)

    AC_MSG_NOTICE([* use 'long double' for FLOAT])
    AC_DEFINE(fFMTe, ["Le"])
    AC_DEFINE(fFMTE, ["LE"])
    AC_DEFINE(fFMTf, ["Lf"])
    AC_DEFINE(fFMTg, ["Lg"])
    AC_DEFINE(fFMTG, ["LG"])
else
    AC_MSG_NOTICE(* use 'double' for 'FLOAT')
    AC_CHECK_SIZEOF([double])

    AC_DEFINE(USE_LONG_DOUBLE, 0)
    AC_DEFINE(fFMTe, ["e"])
    AC_DEFINE(fFMTE, ["E"])
    AC_DEFINE(fFMTf, ["lf"])
    AC_DEFINE(fFMTg, ["g"])
    AC_DEFINE(fFMTG, ["G"])
fi

#--------------------------------------- the end
if test -z "$LINKER"; then
    LINKER="$CC"
fi
LIBS="-lsxamg $LIBS"

AC_OUTPUT
