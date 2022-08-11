#ifndef PRAGRAM_TITLE_H_
#define PRAGRAM_TITLE_H_

#include <string>
const std::string TITLE("This project demonstrates two mainly used multimedia algorithms, namely:\n" 
                        "chroma-keying and 2D convolution. In spite of naive implementations of both,\n"
                        "the project is focused on performance engineering i.e. application of the\n" 
                        "SIMD-extention for accelarating run-time. The program is targeted at x86 \n"
                        "architecture with AVX2 support. Both intrinsics and OpenMP are used.\n\n"
                        "program options");

#endif  //PRAGRAM_TITLE_H_