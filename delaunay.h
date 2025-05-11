#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <vector>
#include "frontier.h"
#include "geometry.h"

std::vector<PolarPoint> calculate_polar_coords(const std::vector<Point> &points, const Point &origin);
void sort_points(std::vector<PolarPoint> &points);
void legalize(const size_t &tri_index, const std::vector<Point> &input_points,
              std::vector<Triangle> &triangles, Frontier &frontier);
std::vector<Triangle> triangulate(const std::vector<Point> &input_points);


#endif //DELAUNAY_H
