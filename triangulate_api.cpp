#include <vector>
#include "geometry.h"
#include "delaunay.h"

extern std::vector<Triangle> triangulate(const std::vector<Point>& input_points);

extern "C" {
    void triangulate_api(const double* coords, int n_points, int* output, int* n_triangles) {
        std::vector<Point> input;
        input.reserve(n_points);

        for (int i = 0; i < n_points; ++i) {
            double x = coords[2 * i];
            double y = coords[2 * i + 1];
            input.emplace_back(x, y);
        }

        const std::vector<Triangle> tris = triangulate(input);

        *n_triangles = static_cast<int>(tris.size());

        for (int i = 0; i < *n_triangles; ++i) {
            output[3 * i + 0] = static_cast<int>(tris[i].v1);
            output[3 * i + 1] = static_cast<int>(tris[i].v2);
            output[3 * i + 2] = static_cast<int>(tris[i].v3);
        }
    }

}
