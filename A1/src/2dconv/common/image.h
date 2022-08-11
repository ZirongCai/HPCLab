#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <stdlib.h>
#include <cstring>

/** TODO: write documentation
 * 
 */
class ImageRGB{
public:
    ImageRGB(size_t width, size_t height, 
             size_t pad_x = 0, size_t pad_y = 0) : m_red(nullptr), 
                                                   m_green(nullptr), 
                                                   m_blue(nullptr),
                                                   m_width(width),
                                                   m_height(height),
                                                   m_pad_width(0),
                                                   m_pad_height(0),
                                                   m_pad_x(pad_x),
                                                   m_pad_y(pad_y) 
    {
        m_pad_width = m_width + 2 * m_pad_x;
        m_pad_height = m_height + 2 * m_pad_y;

        m_red = new float[this->get_real_total_area()];
        m_green = new float[this->get_real_total_area()];
        m_blue = new float[this->get_real_total_area()];

        memset(m_red, 0, sizeof(float) * this->get_real_total_area());
        memset(m_green, 0, sizeof(float) * this->get_real_total_area());
        memset(m_blue, 0, sizeof(float) * this->get_real_total_area());

    }

    ImageRGB(const ImageRGB &other) : m_red(nullptr), 
                                      m_green(nullptr), 
                                      m_blue(nullptr) {
        *this = other;        
    }

    ImageRGB& operator=(const ImageRGB& other) {
        if (this != &other) {
            m_width = other.m_width;
            m_height = other.m_height;
            m_pad_x = other.m_pad_x;
            m_pad_y = other.m_pad_y;
            m_pad_width = other.m_pad_width;
            m_pad_height = other.m_pad_height;


            if (m_red != nullptr) free(m_red);
            if (m_green != nullptr) free(m_green);
            if (m_blue != nullptr) free(m_blue);

            m_red = new float[this->get_real_total_area()];
            m_green = new float[this->get_real_total_area()];
            m_blue = new float[this->get_real_total_area()];

            std::memcpy((void*)m_red, other.m_red, this->get_real_total_area() * sizeof(float));
            std::memcpy((void*)m_green, other.m_green, this->get_real_total_area() * sizeof(float));
            std::memcpy((void*)m_blue, other.m_blue, this->get_real_total_area() * sizeof(float));
        }
        return *this;
    }


    bool operator==(const ImageRGB& rhs) {
        if ((m_width == rhs.m_width) && (m_height == rhs.m_height) \
                                     && (m_pad_width == rhs.m_pad_width) \
                                     && (m_pad_height == rhs.m_pad_height)) {
            return true;
        }
        else {
            return false;
        }
    }

    bool operator!=(const ImageRGB& rhs) {
        return !((*this) == rhs);

    }

    ~ImageRGB() {
        if (m_red != nullptr) delete [] m_red;
        if (m_green != nullptr) delete [] m_green;
        if (m_blue != nullptr) delete [] m_blue;
    }

    size_t get_width() {return m_width;}
    size_t get_height() {return m_height;}
    size_t get_pad_height() {return m_pad_height;}
    size_t get_pad_width() {return m_pad_width;}

    float* get_red() {return m_red;}
    float* get_green() {return m_green;}
    float* get_blue() {return m_blue;}

    size_t get_index(size_t x, size_t y) {return (m_pad_x + x) + (m_pad_y + y) * m_pad_width;}
    size_t get_real_index(size_t x, size_t y) {return x + y * m_pad_width;}

    size_t get_total_area() {return m_height * m_width;}
    size_t get_real_total_area() {return m_pad_height * m_pad_width;}

    const static short num_channels = 3;


private:

    size_t m_height;
    size_t m_width;
    size_t m_pad_x;
    size_t m_pad_y;
    size_t m_pad_height;
    size_t m_pad_width;
    
    float* m_red;
    float* m_green;
    float* m_blue;
};

#endif  // _IMAGE_H_