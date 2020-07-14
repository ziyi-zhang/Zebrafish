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

using namespace std;
using namespace zebrafish;


DECLARE_DIFFSCALAR_BASE();


int main(int argc, char **argv) {

    // logger
    bool is_quiet = false;
    std::string log_file = "Zebrafish_gui.log";
    int log_level = 0;
    Logger::init(!is_quiet, log_file);
    log_level = std::max(0, std::min(6, log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    spdlog::flush_every(std::chrono::seconds(3));
    logger().info("Zebrafish_gui logger initialized.");

    // parse input
    std::string image_path = "";
    unsigned int num_threads = 4; // std::max(1u, std::thread::hardware_concurrency());
    int lsMethod = 2;
    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", image_path, "Input TIFF image to process")->check(CLI::ExistingFile);
    command_line.add_option("-n", num_threads, "Input number of threads");
    command_line.add_option("-l", lsMethod, "Input least square solver method");

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return command_line.exit(e);
    }

    // TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
    logger().info("Desired #threads = {}", num_threads);

    // Start GUI
    GUI gui;
    gui.init(image_path);

    return EXIT_SUCCESS;
}
