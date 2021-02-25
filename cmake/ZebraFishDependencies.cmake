# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(ZebraFishDownloadExternal)

################################################################################
# Required libraries
################################################################################

#Polyfem
zebra_download_polyfem()
add_subdirectory(${ZEBRA_EXTERNAL}/polyfem)

# HighFive
if(NOT TARGET highfive)
  zerba_download_HighFive()
	option(HIGHFIVE_USE_EIGEN ON)
	find_package(HDF5 REQUIRED)
	add_library(highfive INTERFACE)
	target_include_directories(highfive SYSTEM INTERFACE ${ZEBRA_EXTERNAL}/HighFive/include/ ${HDF5_INCLUDE_DIRS})
  target_link_libraries(highfive INTERFACE ${HDF5_LIBRARIES})
endif()

# libigl
if(NOT TARGET igl)
  zebra_download_libigl()
  add_subdirectory(${ZEBRA_EXTERNAL}/libigl EXCLUDE_FROM_ALL)
endif()


if(NOT TARGET TinyTiff)
  # TIFF
  zebra_download_tinytiff()
  add_library(TinyTiff
    ${ZEBRA_EXTERNAL}/TinyTIFF/tinytiffhighrestimer.cpp
    ${ZEBRA_EXTERNAL}/TinyTIFF/tinytiffhighrestimer.h
    ${ZEBRA_EXTERNAL}/TinyTIFF/tinytiffreader.cpp
    ${ZEBRA_EXTERNAL}/TinyTIFF/tinytiffreader.h
    ${ZEBRA_EXTERNAL}/TinyTIFF/tinytiffwriter.cpp
    ${ZEBRA_EXTERNAL}/TinyTIFF/tinytiffwriter.h
    # ${ZEBRA_EXTERNAL}/TinyTIFF/libtiff_tools/libtiff_tools.cpp
    # ${ZEBRA_EXTERNAL}/TinyTIFF/libtiff_tools/libtiff_tools.h
  )
  target_include_directories(TinyTiff SYSTEM INTERFACE ${ZEBRA_EXTERNAL}/TinyTIFF)
endif()

if(NOT TARGET polysolve)
  zebra_download_polysolve()
  add_subdirectory(${ZEBRA_EXTERNAL}/polysolve)
endif()


# spdlog
if(NOT TARGET spdlog::spdlog)
	zebra_download_spdlog()
	add_subdirectory(${ZEBRA_EXTERNAL}/spdlog)
endif()


# LBFGS++
if(NOT TARGET LBFGS)
  zebra_download_LBFGS()
  add_library(LBFGS INTERFACE)
  target_include_directories(LBFGS INTERFACE ${ZEBRA_EXTERNAL}/LBFGS/include)
endif()

# MMG
#zebra_download_mmg()
#option(BUILD_TESTING "Enable/Disable continuous integration" OFF)
#set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
#add_subdirectory(${ZEBRA_EXTERNAL}/mmg)
#add_library(mmg::mmg ALIAS libmmg_a)
#add_library(mmg::mmgs ALIAS libmmgs_a)
#add_library(mmg::mmg2d ALIAS libmmg2d_a)
#add_library(mmg::mmg3d ALIAS libmmg3d_a)


if(NOT TARGET tbb_static)
  zebra_download_tbb()
  set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
  set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)

  add_subdirectory(${ZEBRA_EXTERNAL}/tbb tbb)
  set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")

  target_include_directories(tbb_static SYSTEM PUBLIC ${ZEBRA_EXTERNAL}/tbb/include)
endif()