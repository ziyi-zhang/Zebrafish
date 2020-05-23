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