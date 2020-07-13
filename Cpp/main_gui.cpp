// main gui
#include <zebrafish/autodiff.h>
#include <zebrafish/GUI.h>

#include <math.h>
#include <CLI/CLI.hpp>
#include <string>

using namespace std;
using namespace zebrafish;


DECLARE_DIFFSCALAR_BASE();


int main(int argc, char **argv) {

    // parse input
    std::string image_path = "";
    unsigned int num_threads = 32; // std::max(1u, std::thread::hardware_concurrency());
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

    // Start GUI
    GUI gui;
    gui.init(image_path);

    return EXIT_SUCCESS;
}
