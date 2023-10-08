
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstring>
#include "iou.h"


namespace MathCommon
{

#define maxn 51

const float64_t eps = 1E-8;

int sig(float64_t d)
{
    // d < 0 return -1
    // d = 0 return 0
    // d > 0 return 1
    return (d > eps) - (d < -eps);
}

struct Point
{
    float64_t x, y;
    Point() {}
    Point(float64_t x, float64_t y) : x(x), y(y) {}
    bool operator==(const Point &p) const
    {
        return sig(x - p.x) == 0 && sig(y - p.y) == 0;
    }
};

float64_t cross(Point o, Point a, Point b)
{ //叉积
    return (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
}

float64_t area(Point *ps, int n)
{
    ps[n] = ps[0];
    float64_t res = 0;
    for (int i = 0; i < n; i++)
    {
        res += ps[i].x * ps[i + 1].y - ps[i].y * ps[i + 1].x;
    }
    return res / 2.0;
}

int lineCross(Point a, Point b, Point c, Point d, Point &p)
{
    float64_t s1, s2;
    s1 = cross(a, b, c);
    s2 = cross(a, b, d);
    if (sig(s1) == 0 && sig(s2) == 0)   // a b c d 在一条直线上
        return 2;
    if (sig(s2 - s1) == 0)  // ab与cd平行
        return 0;
    p.x = (c.x * s2 - d.x * s1) / (s2 - s1);
    p.y = (c.y * s2 - d.y * s1) / (s2 - s1);
    return 1;
}

//多边形切割
//用直线ab切割多边形p，切割后的在向量(a,b)的左侧，并原地保存切割结果
//如果退化为一个点，也会返回去,此时n为1
//void polygon_cut(Point*p,int&n,Point a,Point b){
//    static Point pp[maxn];
//    int m=0;p[n]=p[0];
//    for(int i=0;i<n;i++){
//        if(sig(cross(a,b,p[i]))>0) pp[m++]=p[i];
//        if(sig(cross(a,b,p[i]))!=sig(cross(a,b,p[i+1])))
//            lineCross(a,b,p[i],p[i+1],pp[m++]);
//    }
//    n=0;
//    for(int i=0;i<m;i++)
//        if(!i||!(pp[i]==pp[i-1]))
//            p[n++]=pp[i];
//    while(n>1&&p[n-1]==p[0])n--;
//}
void polygon_cut(Point *p, int &n, Point a, Point b, Point *pp)
{
    //    static Point pp[maxn];
    int m = 0;
    p[n] = p[0];
    for (int i = 0; i < n; i++)
    {
        if (sig(cross(a, b, p[i])) > 0)
            pp[m++] = p[i];
        if (sig(cross(a, b, p[i])) != sig(cross(a, b, p[i + 1])))
            lineCross(a, b, p[i], p[i + 1], pp[m++]);
    }
    n = 0;
    for (int i = 0; i < m; i++)
        if (!i || !(pp[i] == pp[i - 1]))
            p[n++] = pp[i];
    while (n > 1 && p[n - 1] == p[0])
        n--;
}

//---------------华丽的分隔线-----------------//
//返回三角形oab和三角形ocd的有向交面积,o是原点//
float64_t intersectArea(Point a, Point b, Point c, Point d)
{
    Point o(0, 0);
    int s1 = sig(cross(o, a, b));
    int s2 = sig(cross(o, c, d));
    if (s1 == 0 || s2 == 0)
        return 0.0; //退化，面积为0
    if (s1 == -1)
        std::swap(a, b);
    if (s2 == -1)
        std::swap(c, d);
    Point p[10] = {o, a, b};
    int n = 3;
    Point pp[maxn];
    polygon_cut(p, n, o, c, pp);
    polygon_cut(p, n, c, d, pp);
    polygon_cut(p, n, d, o, pp);
    float64_t res = fabs(area(p, n));
    if (s1 * s2 == -1)
        res = -res;
    return res;
}

//求两多边形的交面积
float64_t intersectArea(Point *ps1, int n1, Point *ps2, int n2)
{
    if (area(ps1, n1) < 0)
        std::reverse(ps1, ps1 + n1);
    if (area(ps2, n2) < 0)
        std::reverse(ps2, ps2 + n2);
    ps1[n1] = ps1[0];
    ps2[n2] = ps2[0];
    float64_t res = 0;
    for (int i = 0; i < n1; i++)
    {
        for (int j = 0; j < n2; j++)
        {
            res += intersectArea(ps1[i], ps1[i + 1], ps2[j], ps2[j + 1]);
        }
    }
    return res; //assumeresispositive!
}

//float64_t iou_poly(std::vector<float64_t> p, std::vector<float64_t> q)
//{
//    Point ps1[maxn], ps2[maxn];
//    int n1 = 4;
//    int n2 = 4;
//    for (int i = 0; i < 4; i++)
//    {
//        ps1[i].x = p[i * 2];
//        ps1[i].y = p[i * 2 + 1];
//
//        ps2[i].x = q[i * 2];
//        ps2[i].y = q[i * 2 + 1];
//    }
//    float64_t inter_area = intersectArea(ps1, n1, ps2, n2);
//    float64_t union_area = fabs(area(ps1, n1)) + fabs(area(ps2, n2)) - inter_area;
//    float64_t iou = inter_area / union_area;
//
//    //    cout << "inter_area:" << inter_area << endl;
//    //    cout << "union_area:" << union_area << endl;
//    //    cout << "iou:" << iou << endl;
//
//    return iou;
//}

//https://blog.csdn.net/qq_14845119/article/details/91876768
struct MyPoint
{
    float64_t x, y;
};

int dcmp(float64_t x)
{
    if(x > eps) return 1;
    return x < -eps ? -1 : 0;
}
float64_t cross(MyPoint a,MyPoint b,MyPoint c) ///叉积
{
    return (a.x-c.x)*(b.y-c.y)-(b.x-c.x)*(a.y-c.y);
}
MyPoint intersection(MyPoint a,MyPoint b,MyPoint c,MyPoint d)
{
    MyPoint p = a;
    float64_t t =((a.x-c.x)*(c.y-d.y)-(a.y-c.y)*(c.x-d.x))/((a.x-b.x)*(c.y-d.y)-(a.y-b.y)*(c.x-d.x));
    p.x +=(b.x-a.x)*t;
    p.y +=(b.y-a.y)*t;
    return p;
}
//计算多边形面积
float64_t PolygonArea(MyPoint p[], int n)
{
    if(n < 3) return 0.0;
    float64_t s = p[0].y * (p[n - 1].x - p[1].x);
    p[n] = p[0];
    for(int i = 1; i < n; ++ i)
        s += p[i].y * (p[i - 1].x - p[i + 1].x);
    return fabs(s * 0.5);
}
float64_t CPIA(MyPoint a[], MyPoint b[], int na, int nb)//ConvexPolygonIntersectArea
{
    MyPoint p[20], tmp[20];
    int tn, eflag;
    a[na] = a[0], b[nb] = b[0];
    memcpy(p,b,sizeof(MyPoint)*(nb + 1));
    for(int i = 0; i < na && nb > 2; i++)
    {
        int sflag = dcmp(cross(a[i + 1], p[0],a[i]));
        for(int j = tn = 0; j < nb; j++, sflag = eflag)
        {
            if(sflag>=0) tmp[tn++] = p[j];
            eflag = dcmp(cross(a[i + 1], p[j + 1],a[i]));
            if((sflag ^ eflag) == -2)
                tmp[tn++] = intersection(a[i], a[i + 1], p[j], p[j + 1]); ///求交点
        }
        memcpy(p, tmp, sizeof(MyPoint) * tn);
        nb = tn, p[nb] = p[0];
    }
    if(nb < 3) return 0.0;
    return PolygonArea(p, nb);
}

float64_t SPIA(MyPoint a[], MyPoint b[], int na, int nb)
{
    int i, j;
    MyPoint t1[4], t2[4];
    float64_t res = 0, num1, num2;
    a[na] = t1[0] = a[0], b[nb] = t2[0] = b[0];

    for (i = 2; i < na; i++)
    {
        t1[1] = a[i-1], t1[2] = a[i];
        num1 = dcmp(cross(t1[1], t1[2],t1[0]));
        if(num1 < 0) std::swap(t1[1], t1[2]);

        for (j = 2; j < nb; j++)
        {

            t2[1] = b[j - 1], t2[2] = b[j];
            num2 = dcmp(cross(t2[1], t2[2],t2[0]));
            if(num2 < 0) std::swap(t2[1], t2[2]);
            res += CPIA(t1, t2, 3, 3) * num1 * num2;
        }
    }
    return std::abs(res);
}

float64_t calcularea(MyPoint r[])
{
    float64_t d12 = sqrt((r[1].x-r[0].x)*(r[1].x-r[0].x) + (r[1].y-r[0].y)*(r[1].y-r[0].y));
    float64_t d14 = sqrt((r[3].x-r[0].x)*(r[3].x-r[0].x) + (r[3].y-r[0].y)*(r[3].y-r[0].y));
    float64_t d24 = sqrt((r[1].x-r[3].x)*(r[1].x-r[3].x) + (r[1].y-r[3].y)*(r[1].y-r[3].y));
    float64_t d32 = sqrt((r[1].x-r[2].x)*(r[1].x-r[2].x) + (r[1].y-r[2].y)*(r[1].y-r[2].y));
    float64_t d34 = sqrt((r[2].x-r[3].x)*(r[2].x-r[3].x) + (r[2].y-r[3].y)*(r[2].y-r[3].y));

    float64_t p1 = (d12+d14+d24) / 2;
    float64_t p2 = (d24+d32+d34) / 2;
    float64_t s1 = sqrt(p1*(p1-d12)*(p1-d14)*(p1-d24));
    float64_t s2 = sqrt(p2*(p2-d32)*(p2-d34)*(p2-d24));
    return s1 + s2;
}

//多边形求交集
float64_t iou_poly(const std::vector<float64_t>& p, const std::vector<float64_t>& q)
{
    if (p.size() < 8 || q.size() < 8) {
        return 0.0;
    }
    MyPoint p1[10], p2[10];
    for (uint32_t i = 0; i < 4; i++)
    {
        p1[i].x = p[i * 2];
        p1[i].y = p[i * 2 + 1];
        p2[i].x = q[i * 2];
        p2[i].y = q[i * 2 + 1];
    }
    float64_t inter = SPIA(p1, p2, 4, 4);
    float64_t unionArea = calcularea(p1) + calcularea(p2) - inter;
    if (fabs(inter) < eps || fabs(unionArea) < eps) {
        return 0;
    }

    return fabs(inter / unionArea);
}

// 计算重叠面积/自身面积
float64_t overlap_poly(const std::vector<float64_t>& p, const std::vector<float64_t>& q)
{
    if (p.size() < 8 || q.size() < 8) {
        return 0.0;
    }

    MyPoint p1[10],p2[10];
    for (uint32_t i = 0; i < 4; i++)
    {
        p1[i].x = p[i * 2];
        p1[i].y = p[i * 2 + 1];
        p2[i].x = q[i * 2];
        p2[i].y = q[i * 2 + 1];
    }
    float64_t inter = SPIA(p1, p2, 4, 4);
    float64_t area_poly = std::min(calcularea(p1), calcularea(p2));
    if (fabs(inter) < eps || fabs(area_poly) < eps) {
        return 0;
    }
    return fabs(inter / area_poly);
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
bool IsPointInPolygonF32(const Point2D_f32 &point, const std::vector<Point2D_f32> &vecPoints)
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
}
