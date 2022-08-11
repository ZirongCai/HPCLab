#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <errno.h>
#include <cstring>

#include "helper.h"

typedef unsigned char byte; 
const unsigned bmp_header_size = 54;

void report_error(const char* file, size_t string_number) {
    std::cout << "REASON: " <<  strerror(errno) << std::endl;
    std::cout << "OCCURED at: " << file << "::" << string_number << std::endl;
}

CommandLineSettings *CommandLineSettings::m_isntance = nullptr;

/** Read a bmp file and extracts an image
 * 
 * @param path a relative location of an image
 */
ImageRGB read_bmp_image(std::string path, size_t pad_x, size_t pad_y) {

    // open a file
    std::fstream file(path.c_str(), std::ios_base::in | std::ios::binary);
    if (!file.is_open()) {
        std::cout << "ERROR: cannot open a file: " << path << std::endl; REPORT_ERR;
        exit(EXIT_FAILURE);
    } 

    // read header
    byte header[bmp_header_size];
    file.read((char*)header, bmp_header_size);

    // check whether a file is bmp or not
    if ((header[0] != 'B') && (header[1] != 'M')) {
        std::cout << "ERROR: file doesn't belong to the BMP format" << std::endl; REPORT_ERR;
        file.close();
        exit(EXIT_FAILURE);
    }

    // check whether a bmp file has been compressed or not
    const unsigned compression = *reinterpret_cast<unsigned*>(header + 30);
    if (compression) {
        std::cout << "ERROR: file is compressed: " << "cannot work with compressed bmp files"<< std::endl; REPORT_ERR;
        file.close();
        exit(EXIT_FAILURE);
    }

    // check whether an image is encoded according to 8-bit color format
    const unsigned bits_per_pixel = *reinterpret_cast<short*>(header + 28);
    if (bits_per_pixel != (8 * ImageRGB::num_channels)) {
        std::cout << "ERROR: an image is not encoded as 8-bit color format" << std::endl; REPORT_ERR;
        file.close();
        exit(EXIT_FAILURE);
    }

    // retrieve info from the header 

    const unsigned pixel_array_offset = *reinterpret_cast<unsigned*>(header + 10);
    const unsigned image_size = *reinterpret_cast<unsigned*>(header + 34);
    const unsigned height = *reinterpret_cast<unsigned*>(header + 22);
    const unsigned width = *reinterpret_cast<unsigned*>(header + 18);
    const unsigned padded_width = image_size / (3 * height);
    ImageRGB image(width, height, pad_x, pad_y);

    // allocate arrays for each channel 
    float *red = image.get_red();
    float *green = image.get_green();
    float *blue = image.get_blue();
      

    // change the current position to the pixel array within a file
    file.seekg(pixel_array_offset, file.beg);

    // read pixel array
    const size_t pixel_array_size = ImageRGB::num_channels * image.get_total_area();
    byte pixel_array[pixel_array_size];
    file.read((char*)pixel_array, pixel_array_size);

    // close file
    file.close(); 

    // read data from each channel
    for(size_t y = 0; y < height; ++ y) {
        for(size_t x = 0; x < width; ++x) {
            const size_t data_index = ImageRGB::num_channels * (x + y * padded_width);
            const size_t image_index = image.get_index(x, y);

            blue[image_index] = float(pixel_array[data_index]) / 255.0f;
            green[image_index] = float(pixel_array[data_index + 1]) / 255.0f;
            red[image_index] = float(pixel_array[data_index + 2]) / 255.0f;
        }
    }


#ifdef DEBUG
    std::cout << "offset to PixelArray: " << pixel_array_offset << std::endl;
    std::cout << "image width: " << image.get_width() << std::endl;
    std::cout << "padded width: " << padded_width << std::endl;
    std::cout << "image height: " << image.get_height() << std::endl;
    std::cout << "image size: " << image_size << std::endl;
    std::cout << "color planes: " << *reinterpret_cast<short*>(header + 26) << std::endl;
    std::cout << "bits per pixel: " << bits_per_pixel << std::endl;
    std::cout << "image size: " << *reinterpret_cast<unsigned*>(header + 34) << std::endl;
    std::cout << "compression: " << *reinterpret_cast<unsigned*>(header + 30) << std::endl;
#endif

    return image;
}



#include <QApplication>
#include <QPushButton>
#include <QLabel>

Graphics *Graphics::m_instance = nullptr;
Graphics* Graphics::init_graphics(int argc, char **argv) {
    if (m_instance == nullptr) {
        m_instance = new Graphics(argc , argv);
    }
    return m_instance;
}

int Graphics::finish_graphics() {
    if (m_instance != nullptr) {
        delete m_instance;
    }
}

Graphics::Graphics(int argc, char **argv) {
#ifdef GRAPHICS   
    m_app = reinterpret_cast<void *>(new QApplication(argc, argv));
#endif
}

Graphics::~Graphics() {
#ifdef GRAPHICS
    if (m_instance != nullptr) {
        delete reinterpret_cast<QApplication *>(m_app);
    }
#endif
}


void Graphics::display_window(ImageRGB &img) {

#ifdef GRAPHICS    
    QImage myImage(img.get_pad_width(), img.get_pad_height(), QImage::Format_RGB32);
    
    float *red = img.get_red();
    float *green = img.get_green();
    float *blue = img.get_blue();
    
    QColor color;
    for (int y = 0; y < img.get_pad_height(); ++y) {
        for(int x = 0; x < img.get_pad_width(); ++x) {
            size_t coord = img.get_real_index(x, y);
            color.setRgb(int(red[coord] * 255), 
                         int(green[coord] * 255),
                         int(blue[coord] * 255));
            myImage.setPixel(x , img.get_pad_height() - y - 1, color.rgba());
        }
    }


    QLabel myLabel;
    myLabel.setPixmap(QPixmap::fromImage(myImage));

    myLabel.show();

    reinterpret_cast<QApplication *>(m_app)->exec();
#endif
}