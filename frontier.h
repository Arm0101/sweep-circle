#ifndef FRONTIER_H
#define FRONTIER_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "geometry.h"
#include "edgekey.h"

struct FrontierNode {
    size_t vertex_index;
    size_t triangle_index;
    double theta;
    FrontierNode *prev = nullptr;
    FrontierNode *next = nullptr;

    FrontierNode(const size_t vi, const size_t ti, const double th)
        : vertex_index(vi), triangle_index(ti), theta(th) {
    }
};

class Frontier {
public:
    Frontier();

    ~Frontier();

    FrontierNode *head = nullptr;

    void insert_sorted(FrontierNode *node);

    void insert_between(FrontierNode *node, FrontierNode *left, FrontierNode *right, const std::vector<Point> &points);

    void remove(const FrontierNode *node);

    void clear();

    FrontierNode *get_node(const size_t &index);

    void print() const;

    [[nodiscard]] std::pair<std::pair<FrontierNode *, FrontierNode *>, size_t> find_edge(const Point &P, const Point &O,
        const std::vector<Point> &points);

    [[nodiscard]] std::optional<size_t> get_edge(size_t, size_t) const;

    void insert_edge(size_t u, size_t v, size_t t);

    void unmark_edge(size_t u, size_t v);

    [[nodiscard]] bool is_edge(size_t u, size_t v) const;

private:
    std::unordered_set<EdgeKey> edges;
    std::unordered_map<size_t, FrontierNode *> hash_table;

    void insert_to_hash(FrontierNode *node);

    void remove_from_hash(const FrontierNode *node);
};

#endif // FRONTIER_H
