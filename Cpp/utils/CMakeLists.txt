set(SOURCES
    autodiff.h
	TiffReader.cpp
    TiffReader.h
)

prepend_current_path(SOURCES)
zebra_copy_headers(${SOURCES})
zebra_set_source_group(${SOURCES})
target_sources(ZebraFish PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

# add_subdirectory(xxx)

