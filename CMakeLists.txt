cmake_minimum_required(VERSION 3.11)
project(ZebraFish)

################################################################################

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(ZEBRA_TOPLEVEL_PROJECT ON)
else()
    set(ZEBRA_TOPLEVEL_PROJECT OFF)
endif()

# project-options
OPTION(ZEBRA_WITH_UNIT_TESTS "Enable unit test"  ON)

# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG               "Use PNG"                      ON)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            ON)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
################################################################################


option(POLYSOLVE_WITH_CHOLMOD          "Enable Cholmod library"            OFF)
option(POLYSOLVE_WITH_UMFPACK          "Enable UmfPack library"            OFF)
option(POLYSOLVE_WITH_SUPERLU          "Enable SuperLU library"            OFF)
option(POLYSOLVE_WITH_MKL              "Enable MKL library"                OFF)
option(POLYSOLVE_WITH_PARDISO          "Enable Pardiso library"            OFF)
option(POLYSOLVE_WITH_HYPRE            "Enable hypre"                      ON)
option(POLYSOLVE_WITH_AMGCL            "Use AMGCL"                         ON)
option(POLYSOLVE_WITH_SPECTRA          "Enable computing spectrum"         OFF)


# option(POLYFEM_NO_UI "" ON)
option(POLYFEM_WITH_APPS "" OFF)
# option(POLYFEM_WITH_MISC "" OFF)
option(POLYFEM_WITH_AMGCL "" OFF)


option(ZEBRA_ANALYSIS_RUNNABLE         "Compile main analysis"             ON)


### Configuration
set(ZEBRA_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(${PROJECT_NAME}Utils)

################################################################################
# Main Project
################################################################################
add_library(${PROJECT_NAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)

add_subdirectory(Cpp)

# libigl
target_link_libraries(${PROJECT_NAME} PUBLIC igl::png)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)

#Tiny Tiff
target_link_libraries(${PROJECT_NAME} PUBLIC TinyTiff)

#Polysolve
target_link_libraries(${PROJECT_NAME} PUBLIC polysolve)

#spdlog
target_link_libraries(${PROJECT_NAME} PUBLIC spdlog::spdlog)

#LBFGS
target_link_libraries(${PROJECT_NAME} PUBLIC LBFGS)

#TBB
target_link_libraries(${PROJECT_NAME} PUBLIC tbb_static)


################################################################################
# Executable
################################################################################
if(ZEBRA_TOPLEVEL_PROJECT)

  add_executable(${PROJECT_NAME}_bin
    Cpp/main.cpp
  )

  #cli11
  zebra_download_cli11()
  add_subdirectory(${ZEBRA_EXTERNAL}/cli11)

  target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME})
  target_link_libraries(${PROJECT_NAME}_bin PUBLIC CLI11::CLI11)

  ######################################

  if (APPLE)
    set(BUILD_TYPE MACOSX_BUNDLE)
  endif()

  add_executable(${PROJECT_NAME}_gui ${BUILD_TYPE}
    Cpp/main_gui.cpp
    Cpp/gui/GUI.cpp
    Cpp/gui/GUI_stage1.cpp
    Cpp/gui/GUI_stage2.cpp
    Cpp/gui/GUI_stage3.cpp
    Cpp/gui/GUI_stage4.cpp
    Cpp/gui/GUI_stage5.cpp
    Cpp/gui/GUI_stage6.cpp
    Cpp/gui/GUI_stage7.cpp
    Cpp/gui/GUI_stage8.cpp
    Cpp/gui/FileDialog.cpp
    Cpp/gui/GUI_markerDrag.cpp
    Cpp/gui/GUI_mask.cpp
    Cpp/gui/GUI_depthCorrection.cpp
    Cpp/gui/GUI_visualization.cpp
    analysis/zebra-analysis.cpp
  )
  zebra_copy_headers(Cpp/gui/GUI.h Cpp/gui/FileDialog.h)

  target_link_libraries(${PROJECT_NAME}_gui PUBLIC ${PROJECT_NAME})
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC CLI11::CLI11)
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC polyfem igl::tetgen)
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC highfive)
  # target_link_libraries(${PROJECT_NAME}_gui PUBLIC mmg::mmg)

  if(NOT TARGET tinyfiledialogs::tinyfiledialogs)
	  zebra_download_tinyfiledialogs()
	  set(TINYFILEDIALOGS_DIR "${ZEBRA_EXTERNAL}/tinyfiledialogs")
	  add_library(tinyfiledialogs_tinyfiledialogs ${TINYFILEDIALOGS_DIR}/tinyfiledialogs.c)
	  add_library(tinyfiledialogs::tinyfiledialogs ALIAS tinyfiledialogs_tinyfiledialogs)
	  target_include_directories(tinyfiledialogs_tinyfiledialogs SYSTEM INTERFACE ${TINYFILEDIALOGS_DIR})
	  set_target_properties(tinyfiledialogs_tinyfiledialogs PROPERTIES POSITION_INDEPENDENT_CODE ON)
  endif()
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC tinyfiledialogs::tinyfiledialogs)


  if (APPLE)
    # Build an application bundle on OSX
    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)
    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME}_gui)
    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_EXECUTABLE_NAME ${PROJECT_NAME}_gui)
    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_GUI_IDENTIFIER com.im.zebrafish)
    # set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_ICON_FILE im.icns)

    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_SHORT_VERSION_STRING "1.0")
    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_BUNDLE_VERSION "1.0.0")
    set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_COPYRIGHT "© 2022, ZebraFish")


    # set_target_properties(${PROJECT_NAME}_gui PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/MacOSXBundleInfo.plist.in)
    # set_source_files_properties(resources/im.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
  endif()


  # Output directory for binary
  set_target_properties(${PROJECT_NAME}_gui PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
  install(TARGETS ${PROJECT_NAME}_gui DESTINATION .)

  # set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/bundle")

  if(APPLE)
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_DMG_FORMAT "UDBZ")
    set(CPACK_DMG_VOLUME_NAME "${PROJECT_NAME}")
    # set(CPACK_BUNDLE_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)
    # set(CPACK_BUNDLE_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/MacOSXBundleInfo.plist.in)
    # set(CPACK_OSX_PACKAGE_VERSION "10.13")

    set(CPACK_SYSTEM_NAME "OSX")

    # set(CPACK_DMG_DS_STORE "${CMAKE_CURRENT_SOURCE_DIR}/resources/DS_Store")
  elseif(WIN32)
    set(CPACK_GENERATOR "NSIS")
      set(CPACK_SOURCE_GENERATOR "ZIP")
  else()
    set(CPACK_GENERATOR "DEB")
      set(CPACK_SOURCE_GENERATOR "TGZ")
  endif()

  set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};${PROJECT_NAME}_gui;Unspecified;/")

  set(CPACK_PACKAGE_VENDOR "ZebraFish")
  set(CPACK_PACKAGE_NAME "${PROJECT_NAME}_gui")
  set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}")
  set(CPACK_PACKAGE_ICON "${EXTRA_SOURCE}")
  # set(CPACK_DMG_BACKGROUND_IMAGE "${CMAKE_CURRENT_SOURCE_DIR}/resources/splash.png")

  # set(CPACK_PACKAGE_ICON "${EXTRA_SOURCE}")
  # set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/resources/README.md")
  # set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/resources/README.md")
  # set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/resources/README.md")

  # set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/bin/*;${CMAKE_SOURCE_DIR}/lib/*;${CMAKE_SOURCE_DIR}/include/*")
  # set(CPACK_INSTALLED_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/resources/scripts;scripts")

  include(CPack)

endif()

################################################################################
# Tests
################################################################################

if(ZEBRA_WITH_UNIT_TESTS)

  # Catch2
  if(NOT TARGET Catch2::Catch2)
      zebra_download_catch2()
      add_subdirectory(${ZEBRA_EXTERNAL}/Catch2 catch2)
      list(APPEND CMAKE_MODULE_PATH ${ZEBRA_EXTERNAL}/Catch2/contrib)
  endif()

  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()









################################################################################
# Analysis
################################################################################
add_library(ZebraFish-Analysis)
target_include_directories(ZebraFish-Analysis PUBLIC ${PROJECT_BINARY_DIR}/include)
add_subdirectory(analysis)


target_link_libraries(ZebraFish-Analysis PRIVATE warnings::all)
target_link_libraries(ZebraFish-Analysis PUBLIC polyfem igl::tetgen)
target_link_libraries(ZebraFish-Analysis PUBLIC highfive)
# target_link_libraries(ZebraFish-Analysis PUBLIC mmg::mmg)

if(ZEBRA_ANALYSIS_RUNNABLE)

  add_executable(ZebraFish-Analysis_bin
    analysis/main.cpp
  )

  #cli11
  if(NOT TARGET CLI11::CLI11)
    zebra_download_cli11()
    add_subdirectory(${ZEBRA_EXTERNAL}/cli11)
  endif()

  target_link_libraries(ZebraFish-Analysis_bin PUBLIC ZebraFish-Analysis)
  target_link_libraries(ZebraFish-Analysis_bin PUBLIC CLI11::CLI11)
endif()


add_executable(Zebra-Split split/main.cpp Cpp/utils/TiffReader.cpp)
target_link_libraries(Zebra-Split PUBLIC ${PROJECT_NAME})
