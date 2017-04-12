##
## NektarCommon.cmake
##
## Frequently used Nektar++ CMake configuration macros and functions
##

# Attempt to determine architecture for debian/RPM packages
execute_process(COMMAND dpkg --print-architecture
    OUTPUT_VARIABLE DPKG_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND rpm --eval %{_arch}
    OUTPUT_VARIABLE RPM_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

#
# CONSTRUCT_DEBIAN_DEPS(depends outvar)
#
# This macro converts a list of component names to a string containing the
# Debian package dependencies for use by the CPack DEB generator, and stores
# this in an output variable. It assumes that packages will be named
# nektar++-<component>.
#
# Arguments:
#   - `depends`: List of packages.
#   - `outvar`: Name of output variable.
#
MACRO(CONSTRUCT_DEBIAN_DEPS depends outvar)
    SET(${outvar} "")

    FOREACH (pkg ${depends})
        STRING(TOLOWER ${pkg} pkg_lower)

        LIST(FIND NEKTAR_LIBS ${pkg_lower} islib)
        IF(islib EQUAL -1)
            SET(${outvar} "${${outvar}}, nektar++-${pkg_lower} (>= ${NEKTAR_VERSION})")
        ELSE()
            SET(${outvar} "${${outvar}}, libnektar++-${pkg_lower} (>= ${NEKTAR_VERSION})")
        ENDIF()
    ENDFOREACH()

    # Remove starting ", "
    STRING(SUBSTRING ${${outvar}} 2 -1 ${outvar})
    UNSET(pkg)
    UNSET(pkg_lower)
ENDMACRO()

#
# FINALISE_CPACK_COMPONENT(name SUMMARY <summary> DESCRIPTION <description>)
#
# Finalises the variables needed for a component (and only really a component
# containing executables) in order to be packaged by CPack. This should be
# called once all executables and libraries have been added to the
# component. This routine will:
#
# - setup the component's name and description
# - compile a unique list of dependencies
# - construct a Debian dependency string using CONSTRUCT_DEBIAN_DEPS to be used
#   in the resulting `.deb` file.
#
# Arguments:
#   - `name`: component name.
#   - `SUMMARY`: a brief summary of the package
#   - `DESCRIPTION`: a more detailed description of the package
#
MACRO(FINALISE_CPACK_COMPONENT name)
    # Don't both doing anything if we aren't building packages.
    IF (NEKTAR_BUILD_PACKAGES)
        CMAKE_PARSE_ARGUMENTS(COMP "" "DESCRIPTION;SUMMARY" "" ${ARGN})

        # Component names are stored as upper case in the CPack variable names.
        STRING(TOUPPER ${name} COMPVAR)

        # Set the component name to `nektar++-<name>`
        SET(CPACK_COMPONENT_${COMPVAR}_DISPLAY_NAME nektar++-${name}
            CACHE INTERNAL "")
        SET(CPACK_COMPONENT_${COMPVAR}_DESCRIPTION ${COMP_DESCRIPTION}
            CACHE INTERNAL "")
        SET(CPACK_COMPONENT_${COMPVAR}_DESCRIPTION_SUMMARY ${COMP_SUMMARY}
            CACHE INTERNAL "")
        SET(CPACK_RPM_${COMPVAR}_PACKAGE_SUMMARY ${COMP_SUMMARY}
            CACHE INTERNAL "")

        # Remove any duplicates from the existing CPack component dependencies
        # which are set by NEKTAR_ADD_EXECUTABLE and NEKTAR_ADD_LIBRARY
        SET(tmp ${CPACK_COMPONENT_${COMPVAR}_DEPENDS})
        LIST(REMOVE_DUPLICATES tmp)
        SET(CPACK_COMPONENT_${COMPVAR}_DEPENDS ${tmp}
            CACHE INTERNAL "")

        # Construct list of Debian dependencies
        CONSTRUCT_DEBIAN_DEPS("${CPACK_COMPONENT_${COMPVAR}_DEPENDS}" "tmp")
        SET(CPACK_DEBIAN_${COMPVAR}_PACKAGE_DEPENDS ${tmp}
            CACHE INTERNAL "")

        # Other Debian details
        SET(CPACK_DEBIAN_${COMPVAR}_FILE_NAME
            "nektar++-${name}-${NEKTAR_VERSION}-${DPKG_ARCHITECTURE}.deb"
            CACHE INTERNAL "")
        SET(CPACK_RPM_${COMPVAR}_FILE_NAME
            "nektar++-${name}-${NEKTAR_VERSION}-1.${RPM_ARCHITECTURE}.rpm"
            CACHE INTERNAL "")
    ENDIF()
ENDMACRO()

#
# THIRDPARTY_LIBRARY(varname DESCRIPTION <description> [STATIC|SHARED])
#
# Updates a variable containing the name of a third-party shared or static
# library to point to an absolute path defining its location instead of adding
# `-llibname` to the linker flags. This avoids the issue of e.g. linking against
# an outdated system zlib installation.
#
# Arguments:
#   - `varname`: variable name containing the third-party library name. On
#     output will be altered to update the correct path.
#   - `DESCRIPTION`: a brief description of the variable (used in the SET
#     command).
#   - `SHARED`: if the library will be built as a shared library
#   - `STATIC`: if the library will be built as a static library
#
MACRO(THIRDPARTY_LIBRARY varname)
    CMAKE_PARSE_ARGUMENTS(TPLIB "" "DESCRIPTION" "STATIC;SHARED" ${ARGN})

    IF(TPLIB_SHARED)
        SET(LIBTYPE "SHARED")
        SET(TPLIBS ${TPLIB_SHARED})
    ELSEIF(TPLIB_STATIC)
        SET(LIBTYPE "STATIC")
        SET(TPLIBS ${TPLIB_STATIC})
    ENDIF()
    
    FOREACH (lib ${TPLIBS})
        LIST(APPEND tmplist "${TPDIST}/lib/${CMAKE_${LIBTYPE}_LIBRARY_PREFIX}${lib}${CMAKE_${LIBTYPE}_LIBRARY_SUFFIX}")
    ENDFOREACH()

    SET(${varname} ${tmplist} CACHE FILEPATH ${TPLIB_DESCRIPTION} FORCE)
    UNSET(tmplist)
    UNSET(LIBTYPE)
    UNSET(TPLIBS)
    UNSET(lib)
ENDMACRO()

#
# SET_COMMON_PROPERTIES(target)
#
# Sets properties that are common to either library or executable targets. This
# includes:
#
# - Name suffixes: -g for debug, -ms for minsize, -rg for release w/debug.
# - Disable some MSVC compiler warnings
# - Add -pg flag if NEKTAR_ENABLE_PROFILE is switched on and we're using gcc
# - Add compiler definitions and appropriate warning levels to gcc-like
#   compilers (e.g. clang)
# - Define versions for the target
# - Make sure that -fPIC is enabled for library code if building shared
#   libraries.
#
# Arguments:
#   - `target`: target name
#
MACRO(SET_COMMON_PROPERTIES name)
    SET_TARGET_PROPERTIES(${name} PROPERTIES DEBUG_POSTFIX -g)
    SET_TARGET_PROPERTIES(${name} PROPERTIES MINSIZEREL_POSTFIX -ms)
    SET_TARGET_PROPERTIES(${name} PROPERTIES RELWITHDEBINFO_POSTFIX -rg)
    
    IF( MSVC )
        # Disable the warnings about duplicate copy/assignment methods 
        #   (4521, 4522)
        # Disable the warning that arrays are default intialized (4351)	
        # Disable "forcing value to bool 'true' or 'false' (performance
        #   warning)" warning (4800)
        # 4250 - Inheritance via dominance.  Nektar appears to be handling the 
        # diamond correctly.
        # 4373 - Overriding a virtual method with parameters that differ by const
        #        or volatile conforms to the standard.
        # /Za is necessary to prevent temporaries being bound to reference
        #   parameters.
        SET_TARGET_PROPERTIES(${name} PROPERTIES COMPILE_FLAGS 
                                "/wd4521 /wd4522 /wd4351 /wd4018 /wd4800 /wd4250 /wd4373")

        # Enable source level parallel builds.
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    ENDIF( MSVC )	
    
    IF (${CMAKE_COMPILER_IS_GNUCXX})
        IF(NEKTAR_ENABLE_PROFILE)
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
            SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg")
            SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
            SET(LINK_FLAGS "${LINK_FLAGS} -pg")
        ENDIF(NEKTAR_ENABLE_PROFILE)
    ENDIF()

    # Prevent including these common flags multiple times.
    IF (NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES ".*DNEKTAR_DEBUG.*")
        SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DNEKTAR_DEBUG")

        IF ( NEKTAR_FULL_DEBUG )
            SET(CMAKE_CXX_FLAGS_DEBUG 
                    "${CMAKE_CXX_FLAGS_DEBUG} -DNEKTAR_FULLDEBUG")
        ENDIF( NEKTAR_FULL_DEBUG)
   
        IF(NOT MSVC)
            SET(CMAKE_CXX_FLAGS_DEBUG 
                "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-deprecated -Wno-sign-compare")
            SET(CMAKE_CXX_FLAGS_RELEASE 
                "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wno-deprecated -Wno-sign-compare")
            SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO
                "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall -Wno-deprecated -Wno-sign-compare")
            IF (NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                SET(CMAKE_CXX_FLAGS_DEBUG 
                    "${CMAKE_CXX_FLAGS_DEBUG} -fpermissive")
            ENDIF()
        ENDIF (NOT MSVC)

        # Define version
        SET_PROPERTY(TARGET ${name}
            APPEND PROPERTY COMPILE_DEFINITIONS
            NEKTAR_VERSION=\"${NEKTAR_VERSION}\")

        SET(CMAKE_CXX_FLAGS_RELEASE 
                "${CMAKE_CXX_FLAGS_RELEASE} -DNEKTAR_RELEASE")
    ENDIF(NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES ".*DNEKTAR_DEBUG.*")
        
    IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
        # The static libraries must be compiled with position independent
        # code on 64 bit Linux.
        SET_PROPERTY(TARGET ${name} APPEND PROPERTY COMPILE_FLAGS "-fPIC")
    ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
ENDMACRO(SET_COMMON_PROPERTIES name)

#
# ADD_NEKTAR_EXECUTABLE(name COMPONENT <component> [DEPENDS dep1 ...] [SOURCES src1 ...])
#
# Adds a new executable to a component with the supplied component dependencies
# and sources files.
#
# Arguments:
#   - `name`: target name to construct
#   - `COMPONENT`: component name in which this target will live (e.g. demos)
#   - `DEPENDS`: a list of components on which this target depends on
#   - `SOURCES`: a list of source files for this target
#
MACRO(ADD_NEKTAR_EXECUTABLE name)
    CMAKE_PARSE_ARGUMENTS(NEKEXE "" "COMPONENT" "DEPENDS;SOURCES" ${ARGN})
    ADD_EXECUTABLE(${name} ${NEKEXE_SOURCES})
    SET_COMMON_PROPERTIES(${name})

    IF (${CMAKE_SYSTEM} MATCHES "Linux.*")
        SET_PROPERTY(TARGET ${name} APPEND_STRING PROPERTY COMPILE_FLAGS " -pthread")
        SET_PROPERTY(TARGET ${name} APPEND_STRING PROPERTY LINK_FLAGS " -pthread")
    ENDIF()

    STRING(TOLOWER ${NEKEXE_COMPONENT} NEKEXE_COMPONENT)
    STRING(TOUPPER ${NEKEXE_COMPONENT} NEKEXE_COMPVAR)

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${NEKEXE_COMPONENT})
    INSTALL(TARGETS ${name}
        RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL
        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL
        LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKEXE_COMPONENT} OPTIONAL)

    # Add dependencies for executable. We append the dependencies to the CPack
    # component list so that we can resolve these later.
    TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC ${NEKEXE_DEPENDS})
    SET(tmp ${CPACK_COMPONENT_${NEKEXE_COMPVAR}_DEPENDS})
    FOREACH(dep ${NEKEXE_DEPENDS})
        STRING(TOLOWER ${dep} tmp2)
        LIST(APPEND tmp ${tmp2})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES tmp)
    SET(CPACK_COMPONENT_${NEKEXE_COMPVAR}_DEPENDS ${tmp}
        CACHE INTERNAL "")
ENDMACRO()

#
# ADD_NEKTAR_LIBRARY(name
#                    DESCRIPTION <description>
#                    SUMMARY <summary>
#                    DEPENDS dep1 dep2 ...
#                    SOURCES src1 src2 ...
#                    HEADERS head1 head2 ...)
#
# Adds a new library to a component with the supplied component dependencies and
# sources files. A new component will be set up automatically with a lower-case
# name: e.g. if the supplied library name is `LibUtilities` the corresponding
# component is `libutilities`.
#
# Arguments:
#   - `name`: target name to construct
#   - `SUMMARY`: a brief summary of the library
#   - `DESCRIPTION`: a more detailed description of the library
#   - `DEPENDS`: a list of components on which this target depends on
#   - `SOURCES`: a list of source files for this target
#   - `HEADERS`: a list of header files for this target. These will be
#     automatically put into a `dev` package.
#
MACRO(ADD_NEKTAR_LIBRARY name)
    CMAKE_PARSE_ARGUMENTS(NEKLIB "" "DESCRIPTION;SUMMARY" "DEPENDS;SOURCES;HEADERS" ${ARGN})

    ADD_LIBRARY(${name} ${NEKTAR_LIBRARY_TYPE} ${NEKLIB_SOURCES} ${NEKLIB_HEADERS})

    # Infer component name from lower-case library name, variables should use
    # upper-case.
    STRING(TOLOWER ${name} NEKLIB_COMPONENT)
    STRING(TOUPPER ${name} NEKLIB_COMPVAR)

    # Add name to a list so that we know for constructing dependencies.
    SET(NEKTAR_LIBS ${NEKTAR_LIBS} ${NEKLIB_COMPONENT} CACHE INTERNAL "")

    SET_PROPERTY(TARGET ${name} PROPERTY FOLDER ${NEKLIB_COMPONENT})
    SET_PROPERTY(TARGET ${name} PROPERTY VERSION ${NEKTAR_VERSION})

    SET_COMMON_PROPERTIES(${name})

    INSTALL(TARGETS ${name}
        EXPORT Nektar++Libraries
        RUNTIME DESTINATION ${NEKTAR_BIN_DIR} COMPONENT ${NEKLIB_COMPONENT} OPTIONAL
        ARCHIVE DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKLIB_COMPONENT} OPTIONAL
        LIBRARY DESTINATION ${NEKTAR_LIB_DIR} COMPONENT ${NEKLIB_COMPONENT} OPTIONAL)

    FOREACH(HEADER ${NEKLIB_HEADERS})
        STRING(REGEX MATCH "(.*)[/\\]" DIR ${HEADER})
        INSTALL(FILES ${HEADER}
            DESTINATION ${NEKTAR_INCLUDE_DIR}/${name}/${DIR}
            COMPONENT dev)
    ENDFOREACH()

    # Add CPack information
    SET(CPACK_COMPONENT_${NEKLIB_COMPVAR}_DISPLAY_NAME libnektar++-${NEKLIB_COMPONENT}
        CACHE INTERNAL "")
    SET(CPACK_COMPONENT_${NEKLIB_COMPVAR}_DISPLAY_GROUP lib
        CACHE INTERNAL "")
    SET(CPACK_COMPONENT_${NEKLIB_COMPVAR}_DESCRIPTION ${NEKLIB_DESCRIPTION}
        CACHE INTERNAL "")
    SET(CPACK_COMPONENT_${NEKLIB_COMPVAR}_DESCRIPTION_SUMMARY ${NEKLIB_SUMMARY}
        CACHE INTERNAL "")

    # Debian specific information
    SET(CPACK_DEBIAN_${NEKLIB_COMPVAR}_FILE_NAME
        "libnektar++-${NEKLIB_COMPONENT}-${NEKTAR_VERSION}-${DPKG_ARCHITECTURE}.deb"
        CACHE INTERNAL "")
    SET(CPACK_DEBIAN_${NEKLIB_COMPVAR}_PACKAGE_NAME
        "libnektar++-${NEKLIB_COMPONENT}" CACHE INTERNAL "")

    # RPM specific information
    SET(CPACK_RPM_${NEKLIB_COMPVAR}_FILE_NAME
        "libnektar++-${NEKLIB_COMPONENT}-1.${RPM_ARCHITECTURE}.rpm"
        CACHE INTERNAL "")
    SET(CPACK_RPM_${NEKLIB_COMPVAR}_PACKAGE_NAME
        "libnektar++-${NEKLIB_COMPONENT}" CACHE INTERNAL "")

    # If we have dependencies then link against them, and also configure CPack
    # Debian dependencies, which are a special case for some reason. Then set up
    # standard CPack components.
    IF(NEKLIB_DEPENDS)
        TARGET_LINK_LIBRARIES(${name} LINK_PUBLIC ${NEKLIB_DEPENDS})
        CONSTRUCT_DEBIAN_DEPS(${NEKLIB_DEPENDS} "tmp")
        SET(CPACK_DEBIAN_${NEKLIB_COMPVAR}_PACKAGE_DEPENDS ${tmp}
            CACHE INTERNAL "")
        STRING(TOLOWER ${NEKLIB_DEPENDS} tmp)
        SET(CPACK_COMPONENT_${NEKLIB_COMPVAR}_DEPENDS ${tmp}
            CACHE INTERNAL "")
    ENDIF()
ENDMACRO()

#
# ADD_NEKTAR_TEST(name [LENGTHY])
#
# Adds a test with a given name.  The Test Definition File should be in a
# subdirectory called Tests relative to the CMakeLists.txt file calling this
# macros. The test file should be called NAME.tst, where NAME is given as a
# parameter to this macro. If the LENGTHY flag is given, the test will only be
# run if `NEKTAR_TEST_ALL` is enabled.
#
# Arguments:
#   - `name`: name of the test file
#   - `LENGTHY`: denotes a test that requires extended runtime.
#
MACRO(ADD_NEKTAR_TEST name)
    CMAKE_PARSE_ARGUMENTS(NEKTEST "LENGTHY" "" "" ${ARGN})

    IF (NOT NEKTEST_LENGTHY OR NEKTAR_TEST_ALL)
        GET_FILENAME_COMPONENT(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        ADD_TEST(NAME ${dir}_${name}
            COMMAND Tester ${CMAKE_CURRENT_SOURCE_DIR}/Tests/${name}.tst)
    ENDIF()
ENDMACRO(ADD_NEKTAR_TEST)
