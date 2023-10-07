//
// Created by dell on 2023/4/20.
//

#ifndef TEST_MATH_TEST_H
#define TEST_MATH_TEST_H

#endif //TEST_MATH_TEST_H
#include <bits/stdc++.h>

typedef float              float32_t;
typedef double             float64_t;
typedef struct Point2D_f64_ {
    float64_t x;
    float64_t y;
}Point2D_f64;
typedef struct Point2D_f32_ {
    float32_t x;
    float32_t y;
}Point2D_f32;

#define CALC_DIST(point1, point2)  sqrtf(((point1).x - (point2).x)*((point1).x - (point2).x) + ((point1).y - (point2).y)*((point1).y - (point2).y))
#ifndef EPSILON
#define EPSILON (1e-6)
#endif
#define PI                                    (3.141592653589793f)


bool isSegmentCross(Point2D_f64 a,Point2D_f64 b,Point2D_f64 c,Point2D_f64 d)
{
    if (std::max(a.x, b.x) < std::min(c.x, d.x) || std::max(a.y, b.y) < std::min(c.y, d.y) ||
        std::min(a.x, b.x) > std::max(c.x, d.x) || std::min(a.y, b.y) > std::max(c.y, d.y)){
        return false;
    }
    if (((c.x - a.x) * (c.y - d.y) - (c.y - a.y) * (c.x - d.x)) *((c.x - b.x) * (c.y - d.y) - (c.y - b.y) * (c.x - d.x)) <= 0 &&
        ((a.x - c.x) * (a.y - b.y) - (a.y - c.y) * (a.x - b.x)) * ((a.x - d.x) * (a.y - b.y) -(a.y - d.y) * (a.x - b.x)) <= 0){
        return true;
    } else{
        return false;
    }
}

float64_t getDisSquare(Point2D_f64 a, Point2D_f64 b)
{
    return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

bool IsPointInPolygon(const Point2D_f64 &point, const std::vector<Point2D_f64> &vecPoints)
{
    /*
    方法：射线法
        从目标点出发引一条射线，看这条射线和多边形所有边的交点数目。
            如果有奇数个交点，则说明在内部
            如果有偶数个交点，则说明在外部
    */
    /*
    点是否在多边形内的算法逻辑：
        如果从点P作水平向左的射线的话
            （1）若射线与多边形的交点为奇数，则P在多边形内部；
            （1）若射线与多边形的交点为偶数，则P在多边形外部。

        特殊情况（若边为（P1,P2））：
            （1）如果射线正好穿过P1或者P2,那么这个交点会被算作2次，处理办法是如果P的从坐标与P1,P2中较小的纵坐标相同，则直接忽略这种情况；
            （2）如果射线水平，则射线要么与其无交点，要么有无数个，这种情况也直接忽略；
            （3）如果射线竖直，而P的横坐标小于P1,P2的横坐标，则必然相交；
            （4）再判断相交之前，先判断P是否在边(P1,P2)的上面，如果在，则直接得出结论：P再多边形内部。
    */

    bool bResult = false; //判断结果（true；点落在多边形内；false:点未落在多边形内）
    size_t nSize = vecPoints.size();
    size_t j = nSize - 1;//nSize -1 是多边形的最后一个顶点
    for (size_t i = 0; i < nSize; i++) {
        //判断点是否在线段的两侧
        if ((vecPoints[i].y < point.y && vecPoints[j].y >= point.y) ||
            (vecPoints[j].y < point.y && vecPoints[i].y >= point.y)) {
            //根据两点式方程计算出过点P且平行于X轴的直线与线段的交点，两点式方程：x = x1 +  (y - y1) * (x2 - x1) / (y2 - y1);
            if (vecPoints[i].x + (point.y - vecPoints[i].y) * (vecPoints[j].x - vecPoints[i].x) /
                                 (vecPoints[j].y - vecPoints[i].y) < point.x)
                bResult = !bResult;
        }

        //进行下一线段判断
        j = i;
    }

    return bResult;
}

float64_t minDistPointToSegment(Point2D_f64 A, Point2D_f64 B, Point2D_f64 P, Point2D_f64& C)
{
    if (std::abs(A.x - B.x)<EPSILON&&std::abs(A.y - B.y)<EPSILON){
        C = A;
        return CALC_DIST(A,P);
    }
    float64_t r = ((P.x - A.x)*(B.x - A.x) + (P.y - A.y)*(B.y - A.y)) / getDisSquare(A,B);
    if (r <= 0) {           //case 1, returns the length of the AP
        C = A;
        return sqrt(getDisSquare(A, P));
    } else if (r >= 1) {    //case 2, returns the length of the BP
        C = B;
        return sqrt(getDisSquare(B, P));
    } else {                             //case 3, returns the length of the PC
        float64_t AC = r*sqrt(getDisSquare(A,B));  //Find the length of AC, (AC=r*|AB|)
        C.x = A.x + r * (B.x - A.x);
        C.y = A.y + r * (B.y - A.y);
        return sqrt(getDisSquare(A,P)-AC*AC);
    }
}

float64_t distSegmentToRectangle(const std::vector<float64_t>& rect, Point2D_f64 A, Point2D_f64 B)
{
    auto minDist = DBL_MAX;
    if (rect.size() != 8) {
        return minDist;
    }
    std::vector<Point2D_f64_> rectPoints;
    for (auto i = 0; i < rect.size(); i += 2) {
        Point2D_f64 newP;
        newP.x = rect[i];
        newP.y = rect[i + 1];
        rectPoints.emplace_back(newP);
    }

    Point2D_f64 pointC1, pointC2, tempC;   // pointC1 线段上的点，pointC2 矩形上的点，tempC 垂足坐标
    for (auto i = 0; i < 4; ++i) {
        if (isSegmentCross(A, B, rectPoints[i], rectPoints[(i+1)%4])) { // 线段与矩形相交
            return 0;
        }
        auto curDist = minDistPointToSegment(rectPoints[i], rectPoints[(i+1)%4], A, tempC);
        if (curDist < minDist) {
            minDist = curDist;
            pointC1 = A;
            pointC2 = tempC;
        }
        curDist = minDistPointToSegment(rectPoints[i], rectPoints[(i+1)%4], B, tempC);
        if (curDist < minDist) {
            minDist = curDist;
            pointC1 = B;
            pointC2 = tempC;
        }
    }
    for (auto i = 0; i < 4; ++i) {
        if (isSegmentCross(A, B, rectPoints[i], rectPoints[(i+1)%4])) {
            return 0;
        }
        auto curDist = minDistPointToSegment(A, B, rectPoints[i], tempC);
        if (curDist < minDist) {
            minDist = curDist;
            pointC1 = tempC;
            pointC2 = rectPoints[i];
        }
    }
    if (IsPointInPolygon(A, rectPoints) && IsPointInPolygon(B, rectPoints)) {
        minDist *= -1;
    }
    return minDist;
}

void test_distSegmentToRectangle()
{
    std::vector<float64_t> rect1 = {0, 0, 1, 1, 2, 2, 3, 3, 4}; // 最后一个坐标缺失
    Point2D_f64 A1 = {0.5, 0.5};
    Point2D_f64 B1 = {1.5, 1.5};
    auto ret = distSegmentToRectangle(rect1, A1, B1);
    std::cout << ret << std::endl;
    assert(distSegmentToRectangle(rect1, A1, B1) == DBL_MAX);

    std::vector<float64_t> rect2 = {0, 0, 0, 1, 1, 1, 1, 0};
    Point2D_f64 A2 = {-1, 0.5};
    Point2D_f64 B2 = {2, 0.5};
    ret = distSegmentToRectangle(rect2, A2, B2);
    std::cout << ret << std::endl;
    assert(distSegmentToRectangle(rect2, A2, B2) == 0);

    std::vector<float64_t> rect3 = {0, 0, 0, 1, 1, 1, 1, 0};
    Point2D_f64 A3 = {-1, -1};
    Point2D_f64 B3 = {2, -1};
    ret = distSegmentToRectangle(rect3, A3, B3);
    std::cout << ret << std::endl;
    assert(distSegmentToRectangle(rect3, A3, B3) == 1);

    std::vector<float64_t> rect4 = {0, 0, 0, 1, 1, 1, 1, 0};
    Point2D_f64 A4 = {0.5, 0.5};
    Point2D_f64 B4 = {2, 2};
    ret = distSegmentToRectangle(rect4, A4, B4);
    std::cout << ret << std::endl;
    assert(distSegmentToRectangle(rect4, A4, B4) == 0);

    std::vector<float64_t> rect5 = {0, 0, 0, 1, 1, 1, 1, 0};
    Point2D_f64 A5 = {0.5, 0.5};
    Point2D_f64 B5 = {0.7, 0.7};
    ret = distSegmentToRectangle(rect5, A5, B5);
    std::cout << ret << std::endl;
//    assert(distSegmentToRectangle(rect5, A5, B5) == -0.141421);

    std::cout << "test_distSegmentToRectangle" << std::endl;
}


float32_t getAbsAngle(float32_t angle)
{
    if (angle > PI) {
        return angle - 2 * PI;
    } else if (angle < -PI) {
        return 2 * PI + angle;
    } else {
        return angle;
    }
}
// return angle range [0, PI]
float32_t getVector2DAngle(const Point2D_f32& vec1, const Point2D_f32& vec2)
{
    // 计算两个向量的点积
    double dot = vec1.x * vec2.x + vec1.y * vec2.y;
    // 计算两个向量的模长
    double len1 = sqrt(vec1.x * vec1.x + vec1.y * vec1.y);
    double len2 = sqrt(vec2.x * vec2.x + vec2.y * vec2.y);
    // 计算两个向量的夹角（弧度制）
    return acos(fminl(fmaxl(dot / (len1 * len2),-1.0),1.0));
}

double wrapAngle(double angle) {
    double wrappedAngle = fmod(angle + 180.0, 360.0);
    if (wrappedAngle < 0) {
        wrappedAngle += 360.0;
    }
    return wrappedAngle - 180.0;
}

void makePerpendicular(double& angle1, const double angle2) {
    double diff = angle1 - angle2;
    if (diff > M_PI) {
        diff -= 2.0 * M_PI;
    }
    else if (diff < -M_PI) {
        diff += 2.0 * M_PI;
    }
    if (std::abs(diff) > M_PI / 4.0 && std::abs(diff) < 3.0 * M_PI / 4.0) {
        angle1 = angle2 + (diff > 0.0 ? M_PI / 2.0 : -M_PI / 2.0);
    }
}

void test_Angle()
{
    constexpr float32_t  RADIAN_45               = 0.25 * PI;       // unit: rad  45 degrees
    constexpr float32_t  RADIAN_90               = 0.5 * PI;   // unit: rad  45 degrees
    constexpr float32_t  RADIAN_135              = 0.75 * PI;  // unit: rad  45 degrees

    float32_t ussLineAngle = PI * (120) /180;
    float32_t headingAngle = PI * (-76) /180;
    Point2D_f32 ussLineVector, fusObjVector;
    Point2D_f32 ussObjPos1, ussObjPos2;
    ussObjPos1.x = 0;
    ussObjPos1.y = 0;
    ussObjPos2.x = std::cos(ussLineAngle);
    ussObjPos2.y = std::sin(ussLineAngle);

    ussLineVector.x = std::cos(ussLineAngle);
    ussLineVector.y = std::sin(ussLineAngle);
    fusObjVector.x = std::cos(headingAngle);
    fusObjVector.y = std::sin(headingAngle);
    float32_t deltaVectorAngle = getVector2DAngle(ussLineVector, fusObjVector);
    if (deltaVectorAngle < RADIAN_45) {
        headingAngle = ussLineAngle;
    } else if (deltaVectorAngle > RADIAN_135) {
        headingAngle = getAbsAngle(ussLineAngle + PI);
    } else if (deltaVectorAngle < RADIAN_135 && deltaVectorAngle > RADIAN_45) {
        double diff = getAbsAngle(headingAngle - ussLineAngle);
        headingAngle = ussLineAngle + (diff > 0.0 ? PI / 2.0 : -PI / 2.0);
    }

    std::cout << headingAngle / PI * 180 << std::endl;

}