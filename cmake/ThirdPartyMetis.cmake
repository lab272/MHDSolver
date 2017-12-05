########################################################################
#
# ThirdParty configuration for Nektar++
#
# Modified METIS library
#
########################################################################

INCLUDE(ExternalProject)

OPTION(NEKTAR_USE_METIS
    "Use Metis library for performing mesh partitioning." OFF)

IF (NEKTAR_USE_METIS)
    FIND_LIBRARY(METIS_LIBRARY NAMES metis PATHS ${MACPORTS_PREFIX}/lib)
    FIND_PATH   (METIS_INCLUDE_DIR metis.h PATHS ${MACPORTS_PREFIX}/include)

    IF (METIS_LIBRARY AND METIS_INCLUDE_DIR)
        SET(BUILD_METIS OFF)
    ELSE()
        SET(BUILD_METIS ON)
    ENDIF ()

    CMAKE_DEPENDENT_OPTION(THIRDPARTY_BUILD_METIS
        "Build Metis library from ThirdParty" ${BUILD_METIS}
        "NEKTAR_USE_METIS" OFF)

    IF (THIRDPARTY_BUILD_METIS)
        EXTERNALPROJECT_ADD(
            metis-5.1.0
            PREFIX ${TPSRC}
            URL ${TPURL}/metis-5.1.0.tar.gz
            URL_MD5 "5465e67079419a69e0116de24fce58fe"
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/metis-5.1.0
            BINARY_DIR ${TPBUILD}/metis-5.1.0
            TMP_DIR ${TPBUILD}/metis-5.1.0-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST}
                -DCMAKE_C_FLAGS:STRING=-fPIC\ -w
                -DGKLIB_PATH:PATH=${TPSRC}/metis-5.1.0/GKlib
                ${TPSRC}/metis-5.1.0
            )

        IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            # Clang 7.3 has a lovely bug that needs to be patched in order for it to
            # compile.
            IF (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "7.3")
                EXTERNALPROJECT_ADD_STEP(metis-5.1.0 patch-install-path
                    COMMAND sed -i ".bak" "s|#define MAX_JBUFS 128|#define MAX_JBUFS 24|" ${TPSRC}/metis-5.1.0/GKlib/error.c
                    DEPENDERS build
                    DEPENDEES download)
            ENDIF()
        ENDIF()

        THIRDPARTY_LIBRARY(METIS_LIB STATIC metis DESCRIPTION "Metis library")
        MARK_AS_ADVANCED(METIS_LIB)
        MESSAGE(STATUS "Build Metis: ${METIS_LIB}")
        SET(METIS_CONFIG_INCLUDE_DIR ${TPINC})

        INCLUDE_DIRECTORIES(${TPDIST}/include)
    ENDIF()
ENDIF()