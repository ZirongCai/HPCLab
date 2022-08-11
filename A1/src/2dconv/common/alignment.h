#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#define ALIGNMENT 64
inline size_t align(size_t num_bytes) {return ALIGNMENT * ((num_bytes + ALIGNMENT - 1) / ALIGNMENT);}

#endif  // _ALIGNMENT_H_