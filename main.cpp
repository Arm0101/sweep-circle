#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include "delaunay.h"

using namespace std;

void test(const string &filename, int num_points) {
    vector<Point> input_points;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Error opening file: " << endl;
        return;
    }

    int total_points;
    getline(file, line);
    stringstream(line) >> total_points;

    vector<double> times;
    int points_read = 0;

    while (points_read < total_points && getline(file, line)) {
        stringstream ss(line);
        double x, y;
        char delimiter;
        ss >> x >> delimiter >> y;

        input_points.push_back(Point(x, y));
        points_read++;

        if (input_points.size() == num_points) {
            bool duplicated = false;
            for (int i = 0; i < input_points.size(); i++) {
                for (int j = 0; j < input_points.size(); j++) {
                    if (i!=j && input_points[i].x == input_points[j].x && input_points[i].y == input_points[j].y) {
                        duplicated = true;
                        break;
                    }
                }
            }
            if (duplicated) {
                input_points.clear();
                continue;
            }
            auto start = chrono::high_resolution_clock::now();
            const vector<Triangle> triangles = triangulate(input_points);
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = end - start;

            times.push_back(duration.count());
            input_points.clear();
        }
    }

    file.close();

    double total_time = 0.0;
    double max_time = 0.0;
    for (const double time: times) {
        total_time += time;
        if (time > max_time) {
            max_time = time;
        }
    }

    double average_time = total_time / static_cast<double>(times.size());

    ofstream result_file("test/RESULT_SWEEP_CIRCLE.cpu");
    if (result_file.is_open()) {
        result_file << "Numero de casos " << num_points << endl;
        result_file << "Tiempo total " << total_time << endl;
        result_file << "Tiempo promedio " << average_time << endl;
        result_file << "Tiempo max " << max_time << endl;
        result_file.close();
        cout << "Results saved in test/RESULT_SWEEP_CIRCLE.cpu" << endl;
    } else {
        cerr << "Error opening result file." << endl;
    }
}

int main() {
    test("test/input(10000000).pnt", 5);
    return 0;
    std::vector<Point> input_points = {
        {2.5, 3.45},
        {1.97, 0.47},
        {1.95, 6.54},
        {9.80, 7.16},
        {8.88, 8.90},
    };

    // input_points = {
    //     {3.31764999, 5.64832592},
    //     {6.02465924, 8.89429796},
    //     {7.11434308, 4.38600484},
    //     {3.4438854, 6.41768272},
    //     {9.63765061, 1.70740145}
    // };

    // input_points = {
    //     {7.6644719,  7.89719367},
    //     {4.81609377, 2.76744154},
    //     {1.53724756, 8.02677057},
    //     {4.28467069, 3.22088733},
    //     {1.74825772, 3.57974194}
    // };
    //
    // // wrong
    // input_points = {
    //     {9.4868109, 5.91552423},
    //     {4.17838519, 2.07220474},
    //     {5.67091949, 6.16255558},
    //     {2.79771264, 4.76213489},
    //     {9.76113479, 2.9718843}
    // };

    // input_points = {
    //     {7.2993401, 6.07106814},
    //     {3.60568982, 8.82907912},
    //     {7.54667094, 3.63844434},
    //     {2.9477541, 5.93853308},
    //     {1.51804164, 9.61291437}
    // };
    // fatal
    // input_points = {
    //     {6.2923861,5.60464185},
    //     {2.98782685, 8.23955621},
    //     {4.71923568 ,1.81416428},
    //     {2.68017989 ,9.02705192},
    //     {3.44894972 ,7.79835209}
    // };
    // input_points = {
    //     {4.92247161, 1.05552712},
    //     {7.2551506, 3.35334789},
    //     {8.91306828, 8.37383746},
    //     {7.77254616, 3.26478948},
    //     {9.97901544, 7.92557511}
    // };

    // wrong2
    input_points = {
        {7.58409461, 3.79366236},
        {1.7519613, 8.95195416},
        {9.0356713, 3.41837245},
        {7.40950852, 6.31824286},
        {8.03672304, 7.83267412}
    };

    input_points = {
        {2.59142351,5.73375908},
        {9.88665537, 5.75392384},
        {1.9496729,9.08471611},
        {6.27784709, 5.11270312},
        {3.06695222,3.81823609}
    };

    input_points = {
        {8.01127095,9.1860608},
        {9.85125044,8.70437245},
        {3.29020713,8.14990811},
        {1.06586289,7.9885883},
        {9.42337692,9.86050714}
    };

    input_points = {
        {4.97312246,3.65678663},
        {5.86600176,7.67599412},
        {5.32179284,2.96184111},
        {6.67458054,9.13518965},
        {4.89357819,1.27115774}
    };
    input_points = {
        {7.49341328,3.39617861},
        {9.12331007,2.88290379},
        {4.63522836,6.77278332},
        {3.14142263,9.2894635},
        {4.5758336,7.49264293}
    };

 input_points = {
        {8.75351395474891, 2.034964147716257},
        {5.521559562836532, 5.021118233850411},
        {9.085960501445257, 2.613348945040678},
        {3.911571107866353, 4.945148251342728},
        {2.1582050766222554, 6.903733130116025}
    };
    input_points = {
        {9.0637836, 7.25551942},
        {2.40041717, 1.03582779},
        {3.91197039, 2.49222533},
        {3.3938918, 7.29552673},
        {1.73337936, 8.16798663}
    };


    const std::vector<Triangle> triangles = triangulate(input_points);

    cout << "Triangles: " << endl;
    for (const auto &t: triangles) {
        cout << t.v1 << " " << t.v2 << " " << t.v3 << endl;
    }

    return 0;
}
