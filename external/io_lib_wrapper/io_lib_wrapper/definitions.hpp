// definitions

#ifndef _IOLIB_DEFINITIONS_LIB
#define _IOLIB_DEFINITIONS_LIB

#include <array>

typedef uint genomic_coordinate_t;
typedef short depth_t; // can not exceed 65K
typedef uint8_t  mutation_t;
typedef std::array<depth_t, 16> panel_t;

#endif