#include <iostream>
#include "frontier.h"
#include "geometry.h"

constexpr double EPS = 1e-6f;

Frontier::Frontier() = default;

Frontier::~Frontier() {
    clear();
}

void Frontier::clear() {
    if (!head) return;
    FrontierNode *current = head->next;
    while (current != head) {
        FrontierNode *temp = current;
        current = current->next;
        remove_from_hash(temp);
        delete temp;
    }
    remove_from_hash(head);
    delete head;
    head = nullptr;
    hash_table.clear();
    edges.clear();
}

FrontierNode *Frontier::get_node(const size_t &index) {
    const auto it = hash_table.find(index);
    return it != hash_table.end() ? it->second : nullptr;
}

void Frontier::insert_to_hash(FrontierNode *node) {
    hash_table[node->vertex_index] = node;
}

void Frontier::remove_from_hash(const FrontierNode *node) {
    hash_table.erase(node->vertex_index);
}

void Frontier::insert_sorted(FrontierNode *node) {
    if (!head) {
        node->next = node;
        node->prev = node;
        head = node;
        insert_to_hash(node);
        return;
    }

    FrontierNode *current = head;
    do {
        if (node->theta > current->theta) {
            node->prev = current->prev;
            node->next = current;
            current->prev->next = node;
            current->prev = node;

            if (current == head)
                head = node;

            insert_to_hash(node);
            return;
        }
        current = current->next;
    } while (current != head);

    node->prev = head->prev;
    node->next = head;
    head->prev->next = node;
    head->prev = node;
    insert_to_hash(node);
}

void Frontier::insert_between(FrontierNode *node, FrontierNode *left, FrontierNode *right,
                              const std::vector<Point> &points) {
    if (ccw(points[left->vertex_index], points[node->vertex_index], points[right->vertex_index])) {
        // after r
        node->prev = right;
        node->next = right->next;
        right->next->prev = node;
        right->next = node;
    } else {
        node->prev = left;
        node->next = left->next;
        left->next->prev = node;
        left->next = node;
    }

    if (node->theta > head->theta) {
        head = node;
    }
    insert_to_hash(node);
}

void Frontier::remove(const FrontierNode *node) {
    if (node->next == node) {
        head = nullptr;
    } else {
        node->prev->next = node->next;
        node->next->prev = node->prev;
        if (node == head)
            head = node->next;
    }
    remove_from_hash(node);
    delete node;
}

std::pair<std::pair<FrontierNode *, FrontierNode *>, size_t> Frontier::find_edge(const Point &P, const Point &O,
    const std::vector<Point> &points) {
    double distance = -1;
    size_t triangle_index = 0;
    std::pair<std::pair<FrontierNode *, FrontierNode *>, size_t> result = {
        {nullptr, nullptr},
        std::numeric_limits<size_t>::max()
    };
    for (auto const edge: edges) {
        Point A = points[edge.a];
        Point B = points[edge.b];

        if (intersect(P, O, A, B)) {
            FrontierNode *node_a = get_node(edge.a);
            FrontierNode *node_b = get_node(edge.b);
            const double current_dist = point_segment_distance(O, A, B);
            if ((current_dist - distance) > EPS || (current_dist - distance) < EPS && edge.triangle_index >=
                triangle_index) {
                triangle_index = edge.triangle_index;
                distance = current_dist;
                if (node_a->theta < node_b->theta)
                    result = std::make_pair(std::make_pair(node_b, node_a), edge.triangle_index);
                else result = std::make_pair(std::make_pair(node_a, node_b), edge.triangle_index);
            }
        }
    }
    return result;
}

std::optional<size_t> Frontier::get_edge(const size_t a, const size_t b) const {
    if (const auto it = edges.find(EdgeKey(a, b)); it != edges.end()) {
        return it->triangle_index;
    }
    return std::nullopt;
}


void Frontier::insert_edge(const size_t u, const size_t v, const size_t t) {
    edges.insert(EdgeKey(u, v, t));
}

void Frontier::unmark_edge(size_t u, size_t v) {
    edges.erase(EdgeKey(u, v));
}


bool Frontier::is_edge(const size_t u, const size_t v) const {
    return edges.contains(EdgeKey(u, v));
}

void Frontier::print() const {
    if (!head) {
        std::cout << "Frontier is empty" << std::endl;
        return;
    }

    std::cout << "Frontier:" << std::endl;
    const FrontierNode *current = head;
    do {
        std::cout << "Vi: " << current->vertex_index
                << " | θ: " << current->theta
                << " | Ti: " << current->triangle_index
                << "\n";
        current = current->next;
    } while (current != head);
}
