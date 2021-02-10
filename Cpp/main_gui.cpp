// main gui
#include <zebrafish/autodiff.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <math.h>
#include <CLI/CLI.hpp>
#include <string>
#include <algorithm>
#include <sstream>

using namespace zebrafish;


// DECLARE_DIFFSCALAR_BASE();


int main(int argc, char **argv) {

    // parse input
    std::string imagePath = "", maskPath = "", analysisInputPath = "";
    unsigned int num_threads = std::min(32u, std::max(1u, std::thread::hardware_concurrency() - 1));
        // At least 1 thread, at most 32 threads
        // prefer (#TTL - 1)
    int lsMethod = 2;
    int debugMode = 0;
    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", imagePath, "Input TIFF image to process")->check(CLI::ExistingFile);
    command_line.add_option("-m", maskPath, "Input mask TIFF image to process")->check(CLI::ExistingFile);
    command_line.add_option("-n", num_threads, "Input number of threads");
    command_line.add_option("-l", lsMethod, "Input least square solver method");
    command_line.add_option("-b", debugMode);
    command_line.add_option("-a", analysisInputPath)->check(CLI::ExistingFile);

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return command_line.exit(e);
    }

    // logger
    GUI gui;
    int log_level = 0;
    Logger::init(gui.oss);
    log_level = std::max(0, std::min(6, log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    spdlog::flush_every(std::chrono::seconds(3));
    logger().info("Zebrafish_gui logger initialized.");

    // TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
    logger().info("Desired #threads = {}", num_threads);

    // Start GUI
    gui.init(imagePath, maskPath, analysisInputPath, debugMode);

    return EXIT_SUCCESS;
}
