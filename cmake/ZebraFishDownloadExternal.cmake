include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(ZEBRA_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(ZEBRA_EXTRA_OPTIONS "")
endif()

function(zebra_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${ZEBRA_EXTERNAL}/${name}
        DOWNLOAD_DIR ${ZEBRA_EXTERNAL}/.cache/${name}
        QUIET
        ${ZEBRA_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

# libigl
function(zebra_download_libigl)
  zebra_download_project(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        56f129e4403d3b7b04ddd786745fd3d574e95e04
  )
endfunction()

## tinytiff
function(zebra_download_tinytiff)
    zebra_download_project(TinyTIFF
        GIT_REPOSITORY https://github.com/jkriege2/TinyTIFF.git
        GIT_TAG        e392a3e05bed63686cc536b6365f572a67c2e63f
    )
endfunction()

# Catch2 for testing
function(zebra_download_catch2)
    zebra_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.11.1
    )
endfunction()

## CLI11 3-Clause BSD license optional
function(zebra_download_cli11)
    zebra_download_project(cli11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.8.0.tar.gz
        URL_MD5 5e5470abcb76422360409297bfc446ac
    )
endfunction()

## spdlog MIT
function(zebra_download_spdlog)
    zebra_download_project(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG         v1.3.1
    )
endfunction()

## Polysolve MIT
function(zebra_download_polysolve)
    zebra_download_project(polysolve
        GIT_REPOSITORY     https://github.com/polyfem/polysolve.git
        GIT_TAG            358fa9769e1b67c0e7883eb2b27f171ab3b59b62
    )
endfunction()


