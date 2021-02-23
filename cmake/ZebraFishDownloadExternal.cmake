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
        GIT_TAG        v2.13.4
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
        GIT_TAG            b288fe7d52a758fe3594e711bd801530f440ff53
    )
endfunction()


## LBFGS MIT
function(zebra_download_LBFGS)
    zebra_download_project(LBFGS
        GIT_REPOSITORY     https://github.com/yixuan/LBFGSpp.git
        GIT_TAG            f047ef4586869855f00e72312e7b4d78d11694b1
    )
endfunction()

## tinyfiledialogs
function(zebra_download_tinyfiledialogs)
    zebra_download_project(tinyfiledialogs
        GIT_REPOSITORY https://git.code.sf.net/p/tinyfiledialogs/code
        GIT_TAG        511e6500fa9184923d4859e06ee9a6a4e70820c4
    )
endfunction()



## tbb Apache-2.0
function(zebra_download_tbb)
    zebra_download_project(tbb
        GIT_REPOSITORY https://github.com/nTopology/tbb.git
        GIT_TAG        41adc7a7fbe4e6d37fe57186bd85dde99fa61e66
    )
endfunction()

## Polyfem MIT
function(zebra_download_polyfem)
    zebra_download_project(polyfem
        GIT_REPOSITORY https://github.com/polyfem/polyfem.git
        GIT_TAG        180c89157452756ebd39aa8759583999579ab720
    )
endfunction()

## highfive
function(zerba_download_HighFive)
    zebra_download_project(HighFive
        GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
        GIT_TAG        v2.2.2
    )
endfunction()


## mmg
function(zebra_download_mmg)
    zebra_download_project(mmg
        GIT_REPOSITORY https://github.com/MmgTools/mmg.git
        GIT_TAG        88e2dd6cc773c43141b137fd0972c0eb2f4bbd2a
    )
endfunction()
