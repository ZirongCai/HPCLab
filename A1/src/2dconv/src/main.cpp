#include <iostream>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "helper.h"
#include "image.h"
#include "filter.h"

#include "compute_kernels_default.h"
#include "compute_kernels.h"

#include "program_title.h"

int main(int argc, char *argv[])
{
    po::variables_map map;
    try
    {

        po::options_description desc(TITLE);
        
        // define command line options
        desc.add_options()
            ("help,h", "print program description")
            ("input,i", po::value<std::string>()->default_value("../images/input.bmp"), "path to an input image")
            ("threshold,t", po::value<float>()->default_value(0.25), "chrome keying sensetivity")
            ("width,m", po::value<unsigned>()->default_value(5), "width of gaussian filter")
            ("height,k", po::value<unsigned>()->default_value(5), "height of gaussian filter")
            ("repeats,r", po::value<unsigned>()->default_value(100), "number of repeats to measure performance")
            ("mode", po::value<std::string>()->default_value("default"), "modes: default, intrinsics, asm")
        ;

        // parse command line
        po::store(po::parse_command_line(argc, argv, desc), map);
        
        if (map.count("help")) {
            std::cout << desc << std::endl;
            return 1;
        }

        // print parameter settings and check whether the user provide correct data type for each parameter
        std::cout << std::string(30, '-') << "PARAMETERS" << std::string(30, '-') << std::endl;
        std::cout << "input image: " << map["input"].as<std::string>() << std::endl;
        std::cout << "chrome-keying threshold: " << map["threshold"].as<float>() << std::endl;
        std::cout << "kernel width: " << map["width"].as<unsigned>() << std::endl;
        std::cout << "kernel height: " << map["height"].as<unsigned>() << std::endl;
        std::cout << "num repeats: " << map["repeats"].as<unsigned>() << std::endl;
        std::cout << "mode: " << map["mode"].as<std::string>() << std::endl;
        std::cout << std::string(70, '-') << std::endl << std::endl;

    } catch(const po::error err) {
        std::cout << "ERROR: during command line parsing. Please, read the documentation (--help)" << std::endl;
        std::cout << err.what() << std::endl;
        exit(EXIT_FAILURE); 
    }
    
    Filter filter(map["width"].as<unsigned>(), map["height"].as<unsigned>());
    filter.init_with_gaussian();


    CommandLineSettings::init(map["repeats"].as<unsigned>(), map["threshold"].as<float>());

    // load images
    ImageRGB image = read_bmp_image(map["input"].as<std::string>(), 
                                    filter.get_half_width(), 
                                    filter.get_half_height());

    Graphics *graphics = Graphics::init_graphics(argc, argv);
    graphics->display_window(image);


    ImageRGB output(image);

    //execute kernels
    std::string mode(map["mode"].as<std::string>());
    if(!mode.compare(std::string("default"))) {
        conv::apply_default(image, filter, output);
    }
    else if(!mode.compare(std::string("intrinsics"))) {
        conv::apply_simd_intrinsics(image, filter, output);
    }
    else if (!mode.compare(std::string("asm"))) {
        conv::apply_asm(image, filter, output);
    }
    else {
        std::cout << "ERROR: unknwon mode. Please, read the documentation (--help)" << std::endl;
        exit(EXIT_FAILURE); 
    }

    graphics->display_window(output);

    graphics->finish_graphics();
	return 0;
}