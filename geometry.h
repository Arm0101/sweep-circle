#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <cmath>

struct Point {
    double x, y;
};

struct PolarPoint {
    double r;
    double theta;
    std::size_t index;
};

struct Triangle {
    std::size_t v1;
    std::size_t v2;
    std::size_t v3;
    std::vector<std::size_t> neighbors;

    Triangle(const std::size_t a, const std::size_t b, const std::size_t c, const std::vector<std::size_t> &n)
        : v1(a), v2(b), v3(c), neighbors(n) {
    }
};


inline bool ccw(const Point &a, const Point &b, const Point &c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0;
}

inline bool intersect(const Point &A, const Point &B, const Point &C, const Point &D) {
    return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
}

inline bool in_circle(const Point &a, const Point &b, const Point &c, const Point &p) {
    const double orient = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);

    const double dx = a.x - p.x;
    const double dy = a.y - p.y;
    const double ex = b.x - p.x;
    const double ey = b.y - p.y;
    const double fx = c.x - p.x;
    const double fy = c.y - p.y;

    const double ap = dx * dx + dy * dy;
    const double bp = ex * ex + ey * ey;
    const double cp = fx * fx + fy * fy;

    const double det = (dx * (ey * cp - bp * fy)
                        - dy * (ex * cp - bp * fx)
                        + ap * (ex * fy - ey * fx));

    return (orient > 0) ? (det > 0) : (det < 0);
}

inline double orientation(const Point &a, const Point &b, const Point &c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline double distance(const Point &a, const Point &b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

inline double circumradius(const Point &a, const Point &b, const Point &c) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double ex = c.x - a.x;
    const double ey = c.y - a.y;

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;
    const double d = dx * ey - dy * ex;

    if (std::abs(d) < 1e-12) return std::numeric_limits<double>::max();

    const double x = (ey * bl - dy * cl) * 0.5 / d;
    const double y = (dx * cl - ex * bl) * 0.5 / d;
    return x * x + y * y;
}

inline Point circumcenter(const Point &a, const Point &b, const Point &c) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double ex = c.x - a.x;
    const double ey = c.y - a.y;

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;
    const double d = dx * ey - dy * ex;

    const double x = a.x + (ey * bl - dy * cl) * 0.5 / d;
    const double y = a.y + (dx * cl - ex * bl) * 0.5 / d;

    return {x, y};
}

inline Point centroid(const Point &a, const Point &b, const Point &c) {
    return {
        (a.x + b.x + c.x) / 3.0,
        (a.y + b.y + c.y) / 3.0
    };
}

inline double point_segment_distance(const Point &P, const Point &A, const Point &B) {
    double dx = B.x - A.x;
    double dy = B.y - A.y;
    if (dx == 0 && dy == 0) {
        dx = P.x - A.x;
        dy = P.y - A.y;
        return std::sqrt(dx * dx + dy * dy);
    }

    double t = ((P.x - A.x) * dx + (P.y - A.y) * dy) / (dx * dx + dy * dy);
    t = std::max(0.0, std::min(1.0, t));
    const double px = A.x + t * dx;
    const double py = A.y + t * dy;
    dx = P.x - px;
    dy = P.y - py;
    return std::sqrt(dx * dx + dy * dy);
}

#endif //GEOMETRY_H
