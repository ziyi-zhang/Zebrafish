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
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            ON)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         OFF)
################################################################################




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
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::png)

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


  add_executable(${PROJECT_NAME}_gui
    Cpp/main_gui.cpp
    Cpp/gui/GUI.cpp
    Cpp/gui/GUI_stage0.cpp
    Cpp/gui/GUI_stage1.cpp
    Cpp/gui/GUI_stage2.cpp
    Cpp/gui/FileDialog.cpp
  )
  zebra_copy_headers(Cpp/gui/GUI.h Cpp/gui/FileDialog.h)

  target_link_libraries(${PROJECT_NAME}_gui PUBLIC ${PROJECT_NAME})
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC CLI11::CLI11)
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)

  if(NOT TARGET tinyfiledialogs::tinyfiledialogs)
	  zebra_download_tinyfiledialogs()
	  set(TINYFILEDIALOGS_DIR "${ZEBRA_EXTERNAL}/tinyfiledialogs")
	  add_library(tinyfiledialogs_tinyfiledialogs ${TINYFILEDIALOGS_DIR}/tinyfiledialogs.c)
	  add_library(tinyfiledialogs::tinyfiledialogs ALIAS tinyfiledialogs_tinyfiledialogs)
	  target_include_directories(tinyfiledialogs_tinyfiledialogs SYSTEM INTERFACE ${TINYFILEDIALOGS_DIR})
	  set_target_properties(tinyfiledialogs_tinyfiledialogs PROPERTIES POSITION_INDEPENDENT_CODE ON)
  endif()
  target_link_libraries(${PROJECT_NAME}_gui PUBLIC tinyfiledialogs::tinyfiledialogs)

  #

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
