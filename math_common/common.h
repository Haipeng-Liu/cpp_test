
#ifndef __H_MATH_COMMON_
#define __H_MATH_COMMON_

/************************** Include *******************************************/
#include "data_type.h"
#include "iou.h"
#include <vector>
#include <cmath>

namespace MathCommon
{
template <typename T, typename U>
void updateMapBySet(T& t, U& u)
{
    if (t.empty()) {
        return;
    }
    for (auto iter = t.begin(); iter != t.end();) {
        if (u.count(iter->first)) {
            iter++;
        } else {
            iter = t.erase(iter);
        }
    }
}
float64_t getDisSquare(Point2D_f64 a, Point2D_f64 b);

void calcObjBoundary(uint8_t nearSideType, Point2D_f32 objPos, float32_t angle, float32_t length, float32_t width, std::vector<float64_t>& boundary);


bool is_segment_cross(Point2D_f64 a,Point2D_f64 b,Point2D_f64 c,Point2D_f64 d);
bool is_segment_crossF32(Point2D_f32 a,Point2D_f32 b,Point2D_f32 c,Point2D_f32 d);
int32_t infer_sig(float64_t d);
float64_t crossProduct(Point2D_f64_ o, Point2D_f64_ a, Point2D_f64_ b);
float32_t crossProductF32(Point2D_f32 o, Point2D_f32 a, Point2D_f32 b) ;
int32_t lineCross(Point2D_f64_ a, Point2D_f64_ b, Point2D_f64_ c, Point2D_f64_ d, Point2D_f64_ &p);
int32_t lineCrossF32(Point2D_f32 a, Point2D_f32 b, Point2D_f32 c, Point2D_f32 d, Point2D_f32 &p);
// 两个线段ab cd是否相交
bool isSegmentCross(const Point2D_f64& a, const Point2D_f64& b, const Point2D_f64& c, const Point2D_f64& d);
// 线段ab是否与多边形相交
bool isSegmentCrossPolygon(const Point2D_f64& a, const Point2D_f64& b, const std::vector<Point2D_f64>& poly);
// 线段ab在四边形内部的长度
bool lengthSegmentInQuadrangle(const Point2D_f64& a, const Point2D_f64& b, const std::vector<Point2D_f64>& poly, float64_t& length);


float64_t distPointToRectangle(const std::vector<float64_t>& obj, const Point2D_f32& point, Point2D_f64& pointC);
Point2D_f32 getDeltaDistToCenterPoint(uint8_t nearSideType, float32_t angle, float32_t length, float32_t width);

// 线段到矩形距离 return 0：线段与矩形相交，>0：线段到矩形最短距离，<0：线段在矩形内，到矩形的距离, pointC1 线段上的点，pointC2 矩形上的点
float64_t distSegmentToRectangle(const std::vector<float64_t>& rect, const Point2D_f32& segA, const Point2D_f32& segB,
                                 Point2D_f64& pointC1, Point2D_f64& pointC2);

// 点到线段的距离
float64_t minDistPointToSeg(const Point2D_f32& A, const Point2D_f32& B, const Point2D_f32& P);
float64_t minDistPointToSegment(const Point2D_f64& A, const Point2D_f64& B, const Point2D_f64& P, Point2D_f64& C);

// 角度类
// return angle range [-PI, PI]
float32_t getAbsAngle(float32_t angle);
// return angle range [0, PI]
float32_t getVector2DAngle(const Point2D_f32& vec1, const Point2D_f32& vec2);

bool isSamePoint(const Point2D_f32& p1, const Point2D_f32& p2);

void lineFitLeastSquares(const std::vector<Point2D_f64> &points, std::vector<float64_t> &line);
// 点pA到pntStart pntEnd两点所在直线的垂点pFoot
void findFoot(const Point2D_f32& pntStart, const Point2D_f32& pntEnd, const Point2D_f32& pA, Point2D_f32 &pFoot);
void findFootF64(const Point2D_f64& pntSart, const Point2D_f64& pntEnd, const Point2D_f64& pA, Point2D_f64 &pFoot);
// 求两个圆的交点，返回值 0-无交点，1-1个交点，2-2个交点
int32_t findTwoCircleIntersection(const Point2D_f32& circle1, const Point2D_f32& circle2, double r1, double r2, Point2D_f32& pA, Point2D_f32 &pB);

} // namespace MathCommon


#endif //__H_MATH_COMMON_
