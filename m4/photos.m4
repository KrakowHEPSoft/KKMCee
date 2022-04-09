#AC_SEARCH_PHOTOS(actionIfFound, actionIfNotFound)
AC_DEFUN([AC_SEARCH_PHOTOS],[

if test x$with_photos != x && test x$with_photos != xyes ; then
  AC_MSG_NOTICE([Adding $with_photos to search path for Photos++])
  if test -d $with_photos/include && test -d $with_photos/lib ; then
    found_photos=yes
    photos_include=$with_photos/include
    photos_lib=$with_photos/lib
  else
    found_photos=no
  fi
fi

# if we failed to find it in the specified path then check some "standard locations"

if test x$found_photos != xyes ; then
  searchPaths="$(echo \"${CMAKE_PREFIX_PATH}\"|grep -o '[[^:]]*/[[pP]][[hH]][[oO]][[tT]][[oO]][[sS]]/[[^:]]*')" # FCCSW
  searchPaths="$searchPaths /usr/local/PHOTOS /usr/local/photos"
  searchPaths="$searchPaths /opt/local/PHOTOS /opt/local/photos"
  searchPaths="$searchPaths /opt/PHOTOS /opt/photos"
  searchPaths="$searchPaths /usr/local /usr /opt/local /opt"
  for ac_photos_path_tmp in $searchPaths ; do
    AC_MSG_NOTICE([Testing $ac_photos_path_tmp for Photos++...])
    if test -d $ac_photos_path_tmp/include && test -d $ac_photos_path_tmp/lib && test x$found_photos != xyes ; then
      photos_include=$ac_photos_path_tmp/include
      photos_lib=$ac_photos_path_tmp/lib
      if test -f $photos_include/Photos/Photos.h ; then
        AC_MSG_NOTICE([...found])
        found_photos=yes
      fi
    fi
  done
fi

# final check that the headers and libraries are actually there

if test x$found_photos = xyes ; then
  if test -f $photos_include/Photos/Photos.h && test -f $photos_lib/libPhotospp${shrext_cmds} && test -f $photos_lib/libPhotosppHEPEVT${shrext_cmds} && test -f $photos_lib/libPhotosppHepMC3${shrext_cmds}; then
    PHOTOS_LDFLAGS="-L$photos_lib -lPhotosppHepMC3 -lPhotosppHEPEVT -lPhotospp"
    PHOTOS_CPPFLAGS="-I$photos_include"
    PHOTOS_LIBDIR="$photos_lib"
    has_namespace=`grep "namespace Photospp" $photos_include/Photos/Photos.h | head -n 1`
    if test "x$has_namespace" != "x" ; then
      PHOTOS_CPPFLAGS="$PHOTOS_CPPFLAGS -DPHOTOS_HAS_NAMESPACE"
    fi
  else
    AC_MSG_NOTICE([Tried $photos_include/Photos/Photos.h , $photos_lib/libPhotospp${shrext_cmds} , $photos_lib/libPhotosppHEPEVT${shrext_cmds} and $photos_lib/libPhotosppHepMC3${shrext_cmds}])
    found_photos=no
  fi
fi

if test x$found_photos = xyes ; then
  export PHOTOS_LDFLAGS
  export PHOTOS_CPPFLAGS
  export PHOTOS_LIBDIR
  AC_SUBST([PHOTOS_LDFLAGS])
  AC_SUBST([PHOTOS_CPPFLAGS])
  AC_SUBST([PHOTOS_LIBDIR])
  AC_MSG_NOTICE([Found Photos++ libraries and headers])
  AC_MSG_NOTICE([PHOTOS_LDFLAGS = $PHOTOS_LDFLAGS])
  AC_MSG_NOTICE([PHOTOS_CPPFLAGS = $PHOTOS_CPPFLAGS])
  $1
else
  AC_MSG_NOTICE([Could not find Photos++ libraries and headers])
  $2
fi

])

