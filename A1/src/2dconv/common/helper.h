#include <string>
#include <vector>
#include <assert.h>
#include "image.h"

#ifndef _HELPER_H_
#define _HELPER_H_

#define REPORT_ERR report_error(__FILE__, __LINE__)

void report_error(const char* file, size_t string_number);
ImageRGB read_bmp_image(std::string path, size_t pad_x = 0, size_t pad_y = 0);
enum COLORS {RED=0, GREEN=1, BLUE=2};


class Graphics {
    public:
    ~Graphics();
    
    static Graphics* init_graphics(int argc, char **argv);
    static int finish_graphics();
    void display_window(ImageRGB &image);

    private:
    Graphics(int argc, char **argv);
    
    std::vector<void *> m_windows{};
    static Graphics *m_instance;
    void *m_app;
};


class CommandLineSettings {
public:

    ~CommandLineSettings() {
        if (m_isntance != nullptr) {delete m_isntance;}
    }

    size_t get_num_repeats() {return m_num_repeats;}
    float get_ck_threshold() {return m_ck_threshold;}
    std::string get_platform_name() {return m_platform_name;}

    static CommandLineSettings* get_settings() {
        assert((m_isntance != nullptr) && ("the has not been initialized"));
        return m_isntance;
    }

    static void init(size_t repeats, float threshold, std::string platform_name = "") {
        if (m_isntance == nullptr) {
            m_isntance = new CommandLineSettings(repeats, threshold, platform_name);
        }
    }

private:
    CommandLineSettings(size_t repeats,
                        float threshold,
                        std::string platform_name = "") : m_num_repeats(repeats), 
                                                          m_ck_threshold(threshold),
                                                          m_platform_name(platform_name) {}

    size_t m_num_repeats;
    float m_ck_threshold;
    std::string m_platform_name{};  // <-- valid only for OpenCL example
    static CommandLineSettings* m_isntance;
};

#endif  // _HELPER_H_