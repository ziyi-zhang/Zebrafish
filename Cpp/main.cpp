#include <TiffReader.h>

#include <CLI/CLI.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace std;
using namespace zebrafish;

int main(int argc, char **argv)
{
    std::string image_path = "";

    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", image_path, "Input TIFF image to process")->check(CLI::ExistingFile);

    try
    {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
        return command_line.exit(e);
    }

    std::vector<Eigen::MatrixXd> image;
    read_tif_image(image_path, image);


    cout << image.size() << endl;
    cout << image[0] << endl;

    return EXIT_SUCCESS;
}
