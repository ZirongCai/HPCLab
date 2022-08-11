
#ifndef _CHROME_KEYING_H_
#define _CHROME_KEYING_H_

#include "helper.h"
#include "image.h"
#include "filter.h"


namespace conv {
    void apply_simd_intrinsics(ImageRGB &input, Filter &filter, ImageRGB &output);
    void apply_asm(ImageRGB &input, Filter &filter, ImageRGB &output);
}
#endif  // _CHROME_KEYING_H_