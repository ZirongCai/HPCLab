#ifndef _FILTER_H_
#define _FILTER_H_

#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <math.h>

#include "alignment.h"


/**TODO: document the code
 * 
 */
class Filter{
public:
    Filter(size_t width, size_t height) : m_width(width),
                                          m_height(height),
                                          m_data(nullptr) {
        
        assert(((m_width % 2) == 1) && ("ERROR: width must be an odd number"));
        assert(((m_height % 2) == 1) && ("ERROR: height must be an odd number"));

        m_half_width = m_width >> 1;
        m_half_height = m_height >> 1;

        m_data = (float*)aligned_alloc(ALIGNMENT, align(sizeof(float) * this->get_total_size()));
        std::memset(m_data, 0, sizeof(float) * this->get_total_size());
    }

    ~Filter() {
        if (m_data != nullptr) free(m_data);
    }

    Filter(const Filter &other) : m_data(nullptr) {
        (*this) = other;
    }

    Filter& operator=(const Filter &other) {
        if (this != &other) {
            m_width = other.m_width;
            m_height = other.m_height;
            m_half_width = other.m_half_width;
            m_half_height = other.m_half_height;

            if (m_data != nullptr) free(m_data);

            m_data = (float*)aligned_alloc(ALIGNMENT, align(sizeof(float) * this->get_total_size()));
            std::memcpy(m_data, other.m_data, sizeof(float) * this->get_total_size());

        }
    }

    void init_with_gaussian() {
        const float deviation_x = float(m_half_width * m_half_width);
        const float deviation_y = float(m_half_height * m_half_height); 
        
        float sum = 0.0f;
        for (size_t y = 0; y < m_height; ++y) {
            for(size_t x = 0; x < m_width; ++x) {
                const float shift_x = float((m_half_width - x) * (m_half_width - x));
                const float shift_y = float((m_half_height - y) * (m_half_height - y));
               
                float value = exp(-0.5 * ((shift_x / deviation_x) + (shift_y / deviation_y)));
                m_data[get_index(x, y)] = value;
                sum += value;
            }
        }

        // normilization i.e. sum over domain is equal to 1
        float check_sum = 0.0f;
        for (size_t y = 0; y < m_height; ++y) {
            for(size_t x = 0; x < m_width; ++x) {
                m_data[get_index(x, y)] /= sum;
#ifdef DEBUG                
                check_sum += m_data[get_index(x, y)];
                std::cout << m_data[get_index(x, y)] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "check sum: " << check_sum << std::endl; 

#else
            }
        }
#endif
    }

    size_t get_index(size_t x, size_t y) {return x + y * m_width;}
    size_t get_rel_index(int x, int y) {
        // x and y are coordinates relative to the centre of a filter 
        return get_index(m_half_width + x, m_half_height + y);
    }

    size_t get_total_size() {return m_width * m_height;}

    size_t get_width() {return m_width;}
    size_t get_height() {return m_height;}
    size_t get_half_width() {return m_half_width;}
    size_t get_half_height() {return m_half_height;}
    const float *get_data() {return m_data;}

private:
    size_t m_width;
    size_t m_height;
    size_t m_half_width;
    size_t m_half_height;

    float *m_data;
};


#endif  // _FILTER_H_