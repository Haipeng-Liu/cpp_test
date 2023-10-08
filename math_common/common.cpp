
#include "common.h"
#include "iou.h"
#include <cfloat>


namespace MathCommon {

Point2D_f32 getDeltaDistToCenterPoint(uint8_t nearSideType, float32_t angle, float32_t length, float32_t width)
{
    Point2D_f32 deltaPos;
    switch (nearSideType) {
        case 1: //车尾
            deltaPos.x = length * 0.5 * std::cos(angle);
            deltaPos.y = length * 0.5 * std::sin(angle);
            break;
        case 2:
            break;
        case 3: //右侧
            deltaPos.x = width * 0.5 * std::cos(angle + 0.5 * PI);
            deltaPos.y = width * 0.5 * std::sin(angle + 0.5 * PI);
            break;
        case 4:
            break;
        case 5: //车头
            deltaPos.x = length * 0.5 * std::cos(angle + PI);
            deltaPos.y = length * 0.5 * std::sin(angle + PI);
            break;
        case 6:
            break;
        case 7: //左侧
            deltaPos.x = width * 0.5 * std::cos(angle - 0.5 * PI);
            deltaPos.y = width * 0.5 * std::sin(angle - 0.5 * PI);
            break;
        case 8:
            break;
        default:
            break;
    }
    return deltaPos;
}

void calcObjBoundary(uint8_t nearSideType, Point2D_f32 objPos, float32_t angle, float32_t length, float32_t width,
                     std::vector<float64_t> &boundary)
{
    std::vector<float64_t>(8, 0.f).swap(boundary);
    switch (nearSideType) {
        case 1: //车尾
        {
            //左后点
            boundary[0] = objPos.x + width * 0.5 * std::cos(angle + 0.5 * PI);
            boundary[1] = objPos.y + width * 0.5 * std::sin(angle + 0.5 * PI);
            //右后点
            boundary[2] = objPos.x + width * 0.5 * std::cos(angle - 0.5 * PI);
            boundary[3] = objPos.y + width * 0.5 * std::sin(angle - 0.5 * PI);
            //右前点
            auto deltaX = length * std::cos(angle);
            auto deltaY = length * std::sin(angle);
            boundary[4] = boundary[2] + deltaX;
            boundary[5] = boundary[3] + deltaY;
            //左前点
            boundary[6] = boundary[0] + deltaX;
            boundary[7] = boundary[1] + deltaY;
            break;
        }
        case 2: //右后
            break;
        case 3: //右侧
        {
            //右前点
            boundary[0] = objPos.x + length * 0.5 * std::cos(angle);
            boundary[1] = objPos.y + length * 0.5 * std::sin(angle);
            //右后点
            boundary[2] = objPos.x + length * 0.5 * std::cos(angle + PI);
            boundary[3] = objPos.y + length * 0.5 * std::sin(angle + PI);
            //左后点
            auto deltaX = width * std::cos(angle + 0.5 * PI);
            auto deltaY = width * std::sin(angle + 0.5 * PI);
            boundary[4] = boundary[2] + deltaX;
            boundary[5] = boundary[3] + deltaY;
            //左前点
            boundary[6] = boundary[0] + deltaX;
            boundary[7] = boundary[1] + deltaY;
            break;
        }
        case 4: //右前
            break;
        case 5: //车头
        {
            //左前点
            boundary[0] = objPos.x + width * 0.5 * std::cos(angle + 0.5 * PI);
            boundary[1] = objPos.y + width * 0.5 * std::sin(angle + 0.5 * PI);
            //右前点
            boundary[2] = objPos.x + width * 0.5 * std::cos(angle - 0.5 * PI);
            boundary[3] = objPos.y + width * 0.5 * std::sin(angle - 0.5 * PI);
            //右前点
            auto deltaX = length * std::cos(angle + PI);
            auto deltaY = length * std::sin(angle + PI);
            boundary[4] = boundary[2] + deltaX;
            boundary[5] = boundary[3] + deltaY;
            //左前点
            boundary[6] = boundary[0] + deltaX;
            boundary[7] = boundary[1] + deltaY;
            break;
        }
        case 6: //左前
            break;
        case 7: //左侧
        {
            //左前点
            boundary[0] = objPos.x + length * 0.5 * std::cos(angle);
            boundary[1] = objPos.y + length * 0.5 * std::sin(angle);
            //左后点
            boundary[2] = objPos.x + length * 0.5 * std::cos(angle + PI);
            boundary[3] = objPos.y + length * 0.5 * std::sin(angle + PI);
            //右后点
            auto deltaX = width * std::cos(angle - 0.5 * PI);
            auto deltaY = width * std::sin(angle - 0.5 * PI);
            boundary[4] = boundary[2] + deltaX;
            boundary[5] = boundary[3] + deltaY;
            //右前点
            boundary[6] = boundary[0] + deltaX;
            boundary[7] = boundary[1] + deltaY;
            break;
        }
        case 8: //左后
            break;
        case 9: //中心
        {
            float32_t theta = std::atan2(width, length);
            float32_t halfDiagonal = sqrtf(width * width + length * length) * 0.5;
            //左前点
            boundary[0] = objPos.x + halfDiagonal * std::cos(angle + theta);
            boundary[1] = objPos.y + halfDiagonal * std::sin(angle + theta);
            //右前点
            boundary[2] = objPos.x + halfDiagonal * std::cos(angle - theta);
            boundary[3] = objPos.y + halfDiagonal * std::sin(angle - theta);
            //右后点
            boundary[4] = objPos.x + halfDiagonal * std::cos(angle + theta + PI);
            boundary[5] = objPos.y + halfDiagonal * std::sin(angle + theta + PI);
            //左后点
            boundary[6] = objPos.x + halfDiagonal * std::cos(angle - theta + PI);
            boundary[7] = objPos.y + halfDiagonal * std::sin(angle - theta + PI);
            break;
        }
        default:
        {
            float32_t theta = std::atan2(width, length);
            float32_t halfDiagonal = sqrtf(width * width + length * length) * 0.5;
            //左前点
            boundary[0] = objPos.x + halfDiagonal * std::cos(angle + theta);
            boundary[1] = objPos.y + halfDiagonal * std::sin(angle + theta);
            //右前点
            boundary[2] = objPos.x + halfDiagonal * std::cos(angle - theta);
            boundary[3] = objPos.y + halfDiagonal * std::sin(angle - theta);
            //右后点
            boundary[4] = objPos.x + halfDiagonal * std::cos(angle + theta + PI);
            boundary[5] = objPos.y + halfDiagonal * std::sin(angle + theta + PI);
            //左后点
            boundary[6] = objPos.x + halfDiagonal * std::cos(angle - theta + PI);
            boundary[7] = objPos.y + halfDiagonal * std::sin(angle - theta + PI);
            break;
        }
    }
}

//Get the square of the distance between two points
float64_t getDisSquare(Point2D_f64 a, Point2D_f64 b)
{
    return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

// The shortest distance from point P to segment AB
float64_t minDistPointToSegment(const Point2D_f64& A, const Point2D_f64& B, const Point2D_f64& P, Point2D_f64& C)
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
        float64_t temp = getDisSquare(A,P) - AC*AC;
        if (temp < EPSILON) {
            return 0;
        }
        return sqrt(temp);
    }
}

float64_t getDisSquaref32(Point2D_f32 a, Point2D_f32 b)
{
    return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y);
}

float64_t minDistPointToSeg(const Point2D_f32& A, const Point2D_f32& B, const Point2D_f32& P)
{
    if (std::abs(A.x - B.x)<EPSILON&&std::abs(A.y - B.y)<EPSILON){
        return CALC_DIST(A,P);
    }
    float64_t r = ((P.x - A.x)*(B.x - A.x) + (P.y - A.y)*(B.y - A.y)) / getDisSquaref32(A,B);
    if (r <= 0) {           //case 1, returns the length of the AP
        return sqrt(getDisSquaref32(A, P));
    } else if (r >= 1) {    //case 2, returns the length of the BP
        return sqrt(getDisSquaref32(B, P));
    } else {                             //case 3, returns the length of the PC
        float64_t AC = r*sqrt(getDisSquaref32(A,B));  //Find the length of AC, (AC=r*|AB|)
        float64_t temp = getDisSquaref32(A,P) - AC*AC;
        if (temp < EPSILON) {
            return 0;
        }
        return sqrt(temp);
    }
}

bool isSegmentCross(const Point2D_f64& a, const Point2D_f64& b, const Point2D_f64& c, const Point2D_f64& d)
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

bool isSegmentCrossPolygon(const Point2D_f64& a, const Point2D_f64& b, const std::vector<Point2D_f64>& poly)
{
    for (auto i = 0; i < poly.size(); ++i) {
        Point2D_f64 startP, endP;
        startP.x = poly[i].x;
        startP.y = poly[i].y;
        endP.x = poly[(i + 1) % poly.size()].x;
        endP.y = poly[(i + 1) % poly.size()].y;
        if (isSegmentCross(a, b, startP, endP)) {
            return true;
        }
    }

    return IsPointInPolygon(a, poly) || IsPointInPolygon(b, poly);
}

bool lengthSegmentInQuadrangle(const Point2D_f64& a, const Point2D_f64& b, const std::vector<Point2D_f64>& poly, float64_t& length)
{
    uint8_t crossPointNum = 0;
    Point2D_f64 crossP1, crossP2;
    for (auto i = 0; i < poly.size(); ++i) {
        Point2D_f64 startP, endP;
        startP.x = poly[i].x;
        startP.y = poly[i].y;
        endP.x = poly[(i + 1) % poly.size()].x;
        endP.y = poly[(i + 1) % poly.size()].y;
        if (!isSegmentCross(a, b, startP, endP)) {
            continue;
        }
        switch (++crossPointNum) {
            case 1:
                if (1 != lineCross(a, b, startP, endP, crossP1)) return false;
                break;
            case 2:
                if (1 != lineCross(a, b, startP, endP, crossP2)) return false;
                break;
            default:
                return false;
        }
    }

    switch (crossPointNum) {
        case 0:
            length = (IsPointInPolygon(a, poly) || IsPointInPolygon(b, poly)) ? CALC_DIST(a, b):0.f;
            break;
        case 1:
            length = IsPointInPolygon(a, poly) ? CALC_DIST(a, crossP1):CALC_DIST(b, crossP1);
            break;
        case 2:
            length = CALC_DIST(crossP1, crossP2);
            break;
        default:
            return false;
    }

    return true;
}

float64_t distSegmentToRectangle(const std::vector<float64_t>& rect, const Point2D_f32& segA, const Point2D_f32& segB,
                                 Point2D_f64& pointC1, Point2D_f64& pointC2)
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

    Point2D_f64 A, B, tempC;   // pointC1 线段上的点，pointC2 矩形上的点，tempC 垂足坐标
    A.x = segA.x; A.y = segA.y;
    B.x = segB.x; B.y = segB.y;
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

bool is_segment_cross(Point2D_f64 a,Point2D_f64 b,Point2D_f64 c,Point2D_f64 d){
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

bool is_segment_crossF32(Point2D_f32 a,Point2D_f32 b,Point2D_f32 c,Point2D_f32 d){
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



int32_t infer_sig(float64_t d) {
    return (d > EPSILON) - (d < -EPSILON);
}

float64_t crossProduct(Point2D_f64_ o, Point2D_f64_ a, Point2D_f64_ b) { //叉积
    return (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
}

float32_t crossProductF32(Point2D_f32 o, Point2D_f32 a, Point2D_f32 b) { //叉积
    return (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
}

int32_t lineCross(Point2D_f64_ a, Point2D_f64_ b, Point2D_f64_ c, Point2D_f64_ d, Point2D_f64_ &p) {
    float64_t s1, s2;
    s1 = crossProduct(a, b, c);
    s2 = crossProduct(a, b, d);
    if (infer_sig(s1) == 0 && infer_sig(s2) == 0)
        return 2;
    if (infer_sig(s2 - s1) == 0)
        return 0;
    p.x = (c.x * s2 - d.x * s1) / (s2 - s1);
    p.y = (c.y * s2 - d.y * s1) / (s2 - s1);
    return 1;
}

int32_t lineCrossF32(Point2D_f32 a, Point2D_f32 b, Point2D_f32 c, Point2D_f32 d, Point2D_f32 &p) {
    float64_t s1, s2;
    s1 = crossProductF32(a, b, c);
    s2 = crossProductF32(a, b, d);
    if (infer_sig(s1) == 0 && infer_sig(s2) == 0)
        return 2;
    if (infer_sig(s2 - s1) == 0)
        return 0;
    p.x = (c.x * s2 - d.x * s1) / (s2 - s1);
    p.y = (c.y * s2 - d.y * s1) / (s2 - s1);
    return 1;
}

void findFootF64(const Point2D_f64& pntSart, const Point2D_f64& pntEnd, const Point2D_f64& pA, Point2D_f64 &pFoot)
{
    if(pntSart.x == pntEnd.x) {
        pFoot.x = pntSart.x;
        pFoot.y = pA.y;
        return;
    }
    float k = (pntEnd.y - pntSart.y) * 1.0 / (pntEnd.x - pntSart.x);
    float A = k;
    float B = -1.0;
    float C = pntSart.y - k  * pntSart.x;

    pFoot.x = (B * B * pA.x -A * B * pA.y-A * C) / (A * A + B * B);
    pFoot.y = (A * A * pA.y - A * B * pA.x-B * C) / (A * A + B * B);
}


float64_t distPointToRectangle(const std::vector<float64_t>& obj, const Point2D_f32& point, Point2D_f64& pointC)
{
    auto minDist = DBL_MAX;
    if (obj.size() != 8) {
        return minDist;
    }
    std::vector<Point2D_f64_> objPoints;
    for (auto i = 0; i < obj.size(); i += 2) {
        Point2D_f64 newP;
        newP.x = obj[i];
        newP.y = obj[i + 1];
        objPoints.emplace_back(newP);
    }

    Point2D_f64 newP, tempC;
    newP.x = point.x;
    newP.y = point.y;
    for (auto i = 0; i < 4; ++i) {
        auto curDist = minDistPointToSegment(objPoints[i], objPoints[(i+1)%4], newP, tempC);
        if (curDist < minDist) {
            minDist = curDist;
            pointC = tempC;
        }
    }
    if (IsPointInPolygon(newP, objPoints)) {
        minDist *= -1;
    }
    return minDist;
}

// return angle range [-PI, PI]
float32_t getAbsAngle(float32_t angle)
{
    if (angle > PI) {
        return angle - 2 * PI;
    } else if (angle < -PI) {
        return angle + 2 * PI;
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

bool isSamePoint(const Point2D_f32& p1, const Point2D_f32& p2)
{
    return (std::abs(p1.x - p2.x) < EPSILON) && (std::abs(p1.y - p2.y) < EPSILON);
}

void lineFitLeastSquares(const std::vector<Point2D_f64> &points, std::vector<float64_t> &line)
{
    if (points.empty()) {
        return;
    }
    float32_t x = 0, y = 0, x2 = 0, y2 = 0, xy = 0;

    // Calculating the average of x and y...
    for (int i = 0; i < points.size(); i += 1) {
        x += points[i].x;
        y += points[i].y;
        x2 += points[i].x * points[i].x;
        y2 += points[i].y * points[i].y;
        xy += points[i].x * points[i].y;
    }
    float64_t w = (float64_t) points.size();

    x /= w;
    y /= w;
    x2 /= w;
    y2 /= w;
    xy /= w;

    float64_t dx2 = x2 - x * x;
    float64_t dy2 = y2 - y * y;
    float64_t dxy = xy - x * y;

    float64_t t = (float64_t) atan2(2 * dxy, dx2 - dy2) / 2;
    std::vector<float64_t>(5).swap(line);
    line[0] = (float64_t)cos(t);   //line[0],line[1]为x，y的单位方向向量
    line[1] = (float64_t)sin(t);   //line[2],line[3]为直线经过某点的X,Y值
    line[2] = (float64_t)x;
    line[3] = (float64_t)y;
    line[4] = t;
}

void findFoot(const Point2D_f32& pntStart, const Point2D_f32& pntEnd, const Point2D_f32& pA, Point2D_f32 &pFoot)
{
    if(pntStart.x == pntEnd.x) {
        pFoot.x = pntStart.x;
        pFoot.y = pA.y;
        return;
    }
    float k = (pntEnd.y - pntStart.y) * 1.0 / (pntEnd.x - pntStart.x);
    float A = k;
    float B = -1.0;
    float C = pntStart.y - k  * pntStart.x;

    pFoot.x = (B * B * pA.x -A * B * pA.y-A * C) / (A * A + B * B);
    pFoot.y = (A * A * pA.y - A * B * pA.x-B * C) / (A * A + B * B);
}

int32_t findTwoCircleIntersection(const Point2D_f32& circle1, const Point2D_f32& circle2, double r1, double r2, Point2D_f32& pA, Point2D_f32 &pB)
{
    float64_t x1 = circle1.x;
    float64_t y1 = circle1.y;
    float64_t x2 = circle2.x;
    float64_t y2 = circle2.y;
    // 计算两个圆心之间的距离
    float64_t d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));

    // 如果两圆相离或包含，则无交点
    if (d > r1 + r2 || d < abs(r1 - r2)) {
        return 0;
    }

    // 计算交点的坐标
    float64_t a = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
    float64_t h = sqrt(r1 * r1 - a * a);
    float64_t x3 = x1 + a * (x2 - x1) / d;
    float64_t y3 = y1 + a * (y2 - y1) / d;
    pA.x = x3 + h * (y2 - y1) / d;
    pA.y = y3 - h * (x2 - x1) / d;
    pB.x = x3 - h * (y2 - y1) / d;
    pB.y = y3 + h * (x2 - x1) / d;

    if (std::abs(pA.x - pB.x) < EPSILON && std::abs(pA.y - pB.y) < EPSILON) {
        return 1;
    }
    return 2;
}

} // namespace MathCommon

