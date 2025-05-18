#include "delaunay.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <stack>
#include "frontier.h"
#include "geometry.h"

using std::cout;
using std::endl;


constexpr double EPS = 1e-6f;


std::vector<PolarPoint> calculate_polar_coords(const std::vector<Point> &points, const Point &origin) {
    std::vector<PolarPoint> polar_points;
    for (size_t i = 0; i < points.size(); ++i) {
        const double dx = points[i].x - origin.x;
        const double dy = points[i].y - origin.y;
        const double r = std::sqrt(dx * dx + dy * dy);
        double theta = std::atan2(dy, dx);
        if (theta < 0) theta += 2 * M_PI;
        theta = fmod(theta, 2 * M_PI);
        polar_points.emplace_back(r, theta, i);
    }
    return polar_points;
}

void sort_points(std::vector<PolarPoint> &points) {
    std::sort(points.begin(), points.end(), [](const PolarPoint &a, const PolarPoint &b) {
        if (std::abs(a.r - b.r) > EPS)
            return a.r < b.r;
        return a.theta > b.theta;
    });
}

void get_common_edge(const Triangle &t, const Triangle &neighbor, size_t &va, size_t &vb, size_t &vc, size_t &vn) {
    for (int i = 0; i < 3; ++i) {
        if (i == 0) {
            va = t.v1;
            vb = t.v2;
            vc = t.v3;
        } else if (i == 1) {
            va = t.v2;
            vb = t.v3;
            vc = t.v1;
        } else {
            va = t.v3;
            vb = t.v1;
            vc = t.v2;
        }


        // va-vc opposite in neighbor triangle
        if ((neighbor.v1 != vb) && ((neighbor.v2 == vc && neighbor.v3 == va) || (
                                        neighbor.v3 == vc && neighbor.v2 == va))) {
            vn = neighbor.v1;
            return;
        }
        if ((neighbor.v2 != vb) && ((neighbor.v1 == vc && neighbor.v3 == va) || (
                                        neighbor.v3 == vc && neighbor.v1 == va))) {
            vn = neighbor.v2;
            return;
        }
        if ((neighbor.v3 != vb) && ((neighbor.v1 == vc && neighbor.v2 == va) || (
                                        neighbor.v2 == vc && neighbor.v1 == va))) {
            vn = neighbor.v3;
            return;
        }
    }
}

void legalize(const size_t &tri_index, const std::vector<Point> &input_points,
              std::vector<Triangle> &triangles, Frontier &frontier) {
    std::stack<size_t> to_process;
    std::unordered_set<size_t> visited;

    to_process.push(tri_index);

    while (!to_process.empty()) {
        size_t current_index = to_process.top();
        to_process.pop();

        if (visited.contains(current_index)) continue;
        visited.insert(current_index);

        Triangle &t = triangles[current_index];

        for (size_t n: t.neighbors) {
            Triangle &neighbor = triangles[n];
            bool flipped = false;
            size_t va, vb, vc, vn;
            get_common_edge(t, neighbor, va, vb, vc, vn);

            // check vn -> va-vb-vc
            if (in_circle(input_points[t.v1], input_points[t.v3], input_points[t.v2], input_points[vn])) {
                frontier.unmark_edge(va, vc);
                const size_t tri_index1 = current_index;
                const size_t tri_index2 = n;

                std::vector _n1 = {tri_index1};
                if (auto t_1 = frontier.get_edge(vn, va).value(); t_1 != tri_index1)
                    _n1.emplace_back(t_1);
                if (auto t_2 = frontier.get_edge(va, vb).value(); t_2 != tri_index1)
                    _n1.emplace_back(t_2);

                const Triangle new1 = {vn, va, vb, _n1};

                std::vector _n2 = {tri_index2};
                if (auto t_3 = frontier.get_edge(vn, vc).value(); t_3 != tri_index2)
                    _n1.emplace_back(t_3);
                if (auto t_4 = frontier.get_edge(vc, vb).value(); t_4 != tri_index2)
                    _n1.emplace_back(t_4);
                const Triangle new2 = {vn, vb, vc, _n2};

                frontier.insert_edge(vn, vb, tri_index);

                triangles[tri_index2] = new1;
                frontier.unmark_edge(vn, va);
                frontier.unmark_edge(va, vb);
                frontier.insert_edge(vn, va, tri_index2);
                frontier.insert_edge(va, vb, tri_index2);
                triangles[tri_index1] = new2;
                frontier.unmark_edge(vn, vc);
                frontier.unmark_edge(vb, vc);
                frontier.insert_edge(vn, vc, tri_index1);
                frontier.insert_edge(vb, vc, tri_index1);
                flipped = true;
            }

            if (!flipped && !visited.contains(n)) {
                to_process.push(n);
            }
        }
    }
}


std::vector<Triangle> triangulate(const std::vector<Point> &input_points) {
    Frontier frontier;
    std::vector<Triangle> triangles = {};
    std::vector<PolarPoint> points = {};
    // 1. selecting the origin of the polar coordinates O
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    for (const auto &[x, y]: input_points) {
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
    }

    const double cx = (min_x + max_x) / 2.0f;
    const double cy = (min_y + max_y) / 2.0f;

    double min_dist = std::numeric_limits<double>::max();
    Point i0{}, i1{}, i2{};

    size_t i0_index = std::numeric_limits<size_t>::max();
    size_t i1_index = std::numeric_limits<size_t>::max();
    size_t i2_index = std::numeric_limits<size_t>::max();

    // pick a seed point close to the centroid
    for (std::size_t i = 0; i < input_points.size(); i++) {
        if (const double d = distance({cx, cy}, input_points[i]); d < min_dist) {
            i0_index = i;
            min_dist = d;
        }
    }
    i0 = input_points[i0_index];
    min_dist = std::numeric_limits<double>::max();

    // find the point closest to the seed
    for (std::size_t i = 0; i < input_points.size(); i++) {
        if (i == i0_index) continue;
        if (const double d = distance(i0, input_points[i]); d < min_dist && d > 0.0) {
            i1_index = i;
            min_dist = d;
        }
    }
    i1 = input_points[i1_index];

    // find the third point which forms the smallest circumcircle with the first two
    double min_radius = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < input_points.size(); i++) {
        if (i == i0_index || i == i1_index) continue;
        if (const double r = circumradius(i0, i1, input_points[i]); r < min_radius) {
            i2_index = i;
            min_radius = r;
        }
    }
    i2 = input_points[i2_index];

    if (ccw(i0, i1, i2)) {
        std::swap(i1, i2);
        std::swap(i1_index, i2_index);
    }

    // 4. construction of the first triangle
    triangles.emplace_back(Triangle{i0_index, i1_index, i2_index, {}});

    // center O
    const Point O = centroid(input_points[i0_index], input_points[i1_index], input_points[i2_index]);

    //2. calculating the polar coordinates of input points
    points = calculate_polar_coords(input_points, O);

    const double i0_theta = points[i0_index].theta;
    const double i1_theta = points[i1_index].theta;
    const double i2_theta = points[i2_index].theta;

    //3. sorting input points in increasing distances from O
    sort_points(points);

    //5. determination of the initial frontier.

    frontier.insert_sorted(new FrontierNode(i0_index, 0, i0_theta));
    frontier.insert_sorted(new FrontierNode(i1_index, 0, i1_theta));
    frontier.insert_sorted(new FrontierNode(i2_index, 0, i2_theta));

    frontier.insert_edge(i0_index, i1_index, 0);
    frontier.insert_edge(i1_index, i2_index, 0);
    frontier.insert_edge(i2_index, i0_index, 0);

    std::vector<PolarPoint> points_to_process;
    for (const auto &p: points) {
        if (const size_t index = p.index; index == i0_index || index == i1_index || index == i2_index) continue;
        points_to_process.emplace_back(p);
    }
    points.clear();

    // Triangulation
    for (const auto &pp: points_to_process) {
        auto [edg, tri] = frontier.find_edge(pp);
        auto [vl, vr] = edg;
        // create triangle i,L,R and legalize it recursively
        size_t _i0 = pp.index;
        size_t _i1 = vl->vertex_index;
        size_t _i2 = vr->vertex_index;

        const size_t tri_index = triangles.size();

        if (ccw(input_points[_i0], input_points[_i1], input_points[_i2])) {
            std::swap(_i1, _i2);
        }

        triangles.emplace_back(Triangle{_i0, _i1, _i2, {tri}});
        frontier.insert_edge(_i0, _i1, tri_index);
        frontier.insert_edge(_i1, _i2, tri_index);
        frontier.insert_edge(_i2, _i0, tri_index);

        triangles[tri].neighbors.emplace_back(tri_index);

        frontier.insert_between(new FrontierNode(pp.index, tri_index, pp.theta), vl, vr);
        legalize(tri_index, input_points, triangles, frontier);
    }

    // Finalization
    FrontierNode *current = frontier.head;

    bool changed = true;

    while (changed) {
        changed = false;
        if (!current) break;

        FrontierNode *start = current;
        do {
            FrontierNode *a = current;
            FrontierNode *b = a->next;
            FrontierNode *c = b->next;

            if (a->vertex_index == b->vertex_index || b->vertex_index == c->vertex_index || c->vertex_index == a->
                vertex_index) {
                break;
            }

            const Point &pa = input_points[a->vertex_index];
            const Point &pb = input_points[b->vertex_index];
            const Point &pc = input_points[c->vertex_index];

            // const auto _e1 = frontier(pa, pb)
            auto tri1 = frontier.get_edge(a->vertex_index, b->vertex_index);
            auto tri2 = frontier.get_edge(b->vertex_index, c->vertex_index);

            if (ccw(pa, pb, pc) && !frontier.is_edge(a->vertex_index, c->vertex_index)) {
                auto tri = Triangle{
                    a->vertex_index, b->vertex_index, c->vertex_index, {}
                };
                if (tri1.has_value())
                    tri.neighbors.emplace_back(tri1.value());
                if (tri2.has_value())
                    tri.neighbors.emplace_back(tri2.value());

                triangles.emplace_back(tri);

                const size_t tri_index = triangles.size() - 1;
                frontier.insert_edge(a->vertex_index, b->vertex_index, tri_index);
                frontier.insert_edge(b->vertex_index, c->vertex_index, tri_index);
                frontier.insert_edge(c->vertex_index, a->vertex_index, tri_index);

                legalize(tri_index, input_points, triangles, frontier);

                frontier.remove(b);
                current = a;
                changed = true;
                break;
            }
            if (triangles.size() == 5)
                break;
            current = current->next;
        } while (current != start);
    }
    frontier.clear();

    return triangles;
}
