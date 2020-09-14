set(SOURCES
    zebra-analysis.hpp
    zebra-analysis.cpp
)

prepend_current_path(SOURCES)
zebra_copy_headers(${SOURCES})
zebra_set_source_group(${SOURCES})
target_sources(ZebraFish-Analysis PRIVATE ${SOURCES})
