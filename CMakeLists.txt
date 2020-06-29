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
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
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
  # igl::opengl_glfw_imgui

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
