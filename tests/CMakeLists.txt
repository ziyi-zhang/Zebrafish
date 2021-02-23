cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
	main.cpp
	CylinderTest.cpp
	InterpTest.cpp
	PaddingTest.cpp
	WindingTest.cpp
	CageTest.cpp
	PointMeshDistance.cpp
)


add_executable(unit_tests ${test_sources})

foreach(source IN ITEMS ${test_sources})
    source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(unit_tests PUBLIC ZebraFish Catch2::Catch2 warnings::all igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui igl::tetgen)

set(DATA_DIR "$CMAKE_CURRENT_SOURCE_DIR/data/")
target_include_directories(unit_tests PUBLIC ${ZEBRA_EXTERNAL}/Catch2/single_include/catch2)
target_compile_definitions(unit_tests PUBLIC -DZEBRA_DATA_DIR=\"${DATA_DIR}\")


# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
