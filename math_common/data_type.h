
#ifndef TEST_DATA_TYPE_H
#define TEST_DATA_TYPE_H


typedef float              float32_t;
typedef double             float64_t;

typedef signed char        int8_t;
typedef unsigned char      uint8_t;
typedef signed short       int16_t;
typedef unsigned short     uint16_t;
typedef signed int         int32_t;
typedef unsigned int       uint32_t;
typedef unsigned long      uint32_lt;
typedef unsigned long      int32_lt;
typedef signed long long   int64_lt;
typedef unsigned long long uint64_lt;
typedef char               char_t;
typedef float              float32_t;
typedef double             float64_t;
typedef unsigned char      uchar;

#define PI                                    (3.141592653589793f)
#ifndef EPSILON
#define EPSILON (1e-6)
#endif

#define CALC_DIST(point1, point2)  sqrtf(((point1).x - (point2).x)*((point1).x - (point2).x) + ((point1).y - (point2).y)*((point1).y - (point2).y))

typedef struct Point2D_f64_ {
    float64_t x;
    float64_t y;
}Point2D_f64;
typedef struct Point2D_f32_ {
    float32_t x;
    float32_t y;
}Point2D_f32;

#endif //TEST_DATA_TYPE_H
