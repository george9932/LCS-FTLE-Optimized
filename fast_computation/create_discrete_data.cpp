#include <iostream>
#include <filesystem>
#include "../src/field.hpp"
#include "../src/flow.hpp"
#include "../src/velocity_function.hpp"
#include "../src/io.hpp"
#include "../src/json.hpp"
// Command to compile the C++ code from the lcs directory:
//  g++ -g -Wall -O2 -std=c++17 -fopenmp fast_computation/create_discrete_data.cpp -o bin/CreateDiscreteData

using namespace LCS;

int GetPrecision(double num);

int main() {
    const std::string proj_dir = "fast_computation";
    const std::string json_file = "sim_params.json";

    // Open the JSON file
    std::string json_path = (std::filesystem::path(proj_dir) / json_file).string();
    std::ifstream file(json_path);
    if (!file.is_open()) {
        std::cerr << "Error: unable to open file '" << json_file << "'" << std::endl;
        return 1;
    }
    // Parse the JSON file
    nlohmann::json json;
    file >> json;

    // Read the settings
    const double x_min = json["x_min"];
    const double x_max = json["x_max"];
    const double y_min = json["y_min"];
    const double y_max = json["y_max"];
    const double data_nx = json["data_nx"];
    const double data_ny = json["data_ny"];
    const double t_min = json["t_min"];
    const double t_max = json["t_max"];
    const double data_delta_t = json["data_delta_t"];
    const std::string file_prefix = json["file_prefix"];

    int precision = GetPrecision(data_delta_t);

    std::string data_path = (std::filesystem::path(proj_dir) / "data").string();

    Position<double,2> pos(data_nx,data_ny);
    pos.SetAll(x_min,x_max,y_min,y_max);

    ContinuousVelocity<double, VelocityFunction::DoubleGyreModel<double>,2> double_gyre_vel(data_nx,data_ny,pos);

    double current_time = t_min;
    while (t_max - current_time >= -1e-12) {
        double_gyre_vel.UpdateTime(current_time);
        std::stringstream ss;
        ss << file_prefix << std::fixed << std::setprecision(precision) << double_gyre_vel.GetTime() << ".txt";
        std::string s = (std::filesystem::path(data_path) / ss.str()).string();
        double_gyre_vel.SetAll();
        double_gyre_vel.WriteToFile(s);
        current_time += data_delta_t;
    }
    std::cout << "Discrete " << file_prefix << " data written to files" << std::endl;

    return 0;
}

int GetPrecision(double num) {
    int precision = 0;
    while (std::abs(num - static_cast<long long>(num)) > 0) {
        num *= 10;
        precision++;
    }
    return precision;
}