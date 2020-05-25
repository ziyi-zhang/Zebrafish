set(SOURCES
	Cylinder.cpp
    Cylinder.h
)

prepend_current_path(SOURCES)
zebra_copy_headers(${SOURCES})
zebra_set_source_group(${SOURCES})
target_sources(ZebraFish PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(utils)


# if(POLYFEM_TOPLEVEL_PROJECT AND NOT POLYFEM_NO_UI)
# 	add_subdirectory(viewer)
# endif()
