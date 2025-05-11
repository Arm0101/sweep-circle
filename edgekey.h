#ifndef EDGE_KEY_H
#define EDGE_KEY_H

#include <functional>

struct EdgeKey {
    size_t a, b;
    size_t triangle_index;


    EdgeKey(const size_t _a, const size_t _b, size_t _triangle_index = static_cast<size_t>(-1)): triangle_index(
        _triangle_index) {
        if (_a < _b) {
            a = _a;
            b = _b;
        } else {
            a = _b;
            b = _a;
        }
    }

    bool operator==(const EdgeKey &other) const {
        return a == other.a && b == other.b;
    }
};

template<>
struct std::hash<EdgeKey> {
    size_t operator()(const EdgeKey &e) const noexcept {
        return hash<size_t>()(e.a) ^ (hash<size_t>()(e.b) << 1);
    }
};

#endif // EDGE_KEY_H
