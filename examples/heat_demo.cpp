// Portable demonstration of heat spreading on a one-dimensional insulated grid.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<double> diffuse(const std::vector<double>& temperature, double factor) {
    if (temperature.size() < 2 || factor < 0 || factor > 0.5)
        throw std::invalid_argument("need at least two cells and stable factor in [0, 0.5]");
    std::vector<double> next(temperature.size());
    for (std::size_t i = 0; i < temperature.size(); ++i) {
        const double left = i ? temperature[i - 1] : temperature[i];
        const double right = i + 1 < temperature.size() ? temperature[i + 1] : temperature[i];
        next[i] = temperature[i] + factor * (left - 2 * temperature[i] + right);
    }
    return next;
}

int main(int argc, char** argv) {
    try {
        const int cells = argc > 1 ? std::stoi(argv[1]) : 21;
        const int steps = argc > 2 ? std::stoi(argv[2]) : 100;
        const std::string path = argc > 3 ? argv[3] : "heat.csv";
        if (cells < 2 || cells > 100000 || steps < 0 || steps > 100000)
            throw std::invalid_argument("invalid cells or steps");
        std::vector<double> temperature(cells, 0.0);
        temperature[cells / 2] = 100.0;
        const double initial_sum = std::accumulate(temperature.begin(), temperature.end(), 0.0);
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot open output file");
        output << "step,cell,temperature\n";
        for (int step = 0; step <= steps; ++step) {
            for (int cell = 0; cell < cells; ++cell)
                output << step << ',' << cell << ',' << temperature[cell] << '\n';
            if (step < steps) temperature = diffuse(temperature, 0.25);
        }
        const double final_sum = std::accumulate(temperature.begin(), temperature.end(), 0.0);
        const double error = std::abs(initial_sum - final_sum);
        std::cout << "heat conservation error: " << error << '\n';
        return error < 1e-8 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
