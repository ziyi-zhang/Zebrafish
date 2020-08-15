set(SOURCES
    autodiff.h
    base64Layer.cpp
    base64Layer.h
	TiffReader.cpp
    TiffReader.h
    Logger.cpp
    Logger.hpp
    VTUwriter.cpp
    VTUwriter.h
)

prepend_current_path(SOURCES)
zebra_copy_headers(${SOURCES})
zebra_set_source_group(${SOURCES})
target_sources(ZebraFish PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

# add_subdirectory(xxx)
# add_subdirectory(LBFGSpp)
