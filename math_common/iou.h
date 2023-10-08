
#ifndef POLYIOU_POLYIOU_H
#define POLYIOU_POLYIOU_H

#include "data_type.h"
#include <vector>

namespace MathCommon
{

float64_t iou_poly(const std::vector<float64_t>& p, const std::vector<float64_t>& q);
float64_t overlap_poly(const std::vector<float64_t>& p, const std::vector<float64_t>& q);
bool IsPointInPolygon(const Point2D_f64 &point, const std::vector<Point2D_f64> &vecPoints);
bool IsPointInPolygonF32(const Point2D_f32 &point, const std::vector<Point2D_f32> &vecPoints);

} // namespace perception


#endif //POLYIOU_POLYIOU_H
