// /**
//  * @file   discrete_fast_lcs_computation.cpp
//  * @brief  LCS-FTLE-Optimized is a library for performing fast Lagrangian Coherent Structure (LCS) analysis of flow fields. 
//  *         It is a fork of the lcs library (https://github.com/stevenliuyi/lcs) developed by Yi Liu. 
//  *         OpenMP (http://www.openmp.org/) is supported for parallelization. 
//  *         It implements particle advection calculations for discrete flow velocity data, both forward and backward in time, and then uses the Finite-Time Lyapunov Exponents (FTLE) to obtain LCSs.  
//  *         A faster method that the classic FTLE calculations is used in this library. 
//  *         The unidirectional method described in the paper "Fast computation of finite-time Lyapunov exponent fields for unsteady flows" (https://doi.org/10.1063/1.3270044) by Steven Brunton and Clarence Rowley is implemented here. 
//  *         The method ensures that unnecessary particle integrations are not repeated, leading to significant computational time savings.
//  * @author Georgios Mitsos 
//  * @date   21/08/2024
//  */

#include <iostream>
#include <filesystem>
#include "../src/ftle.hpp"
#include "../src/io.hpp"
#include "../src/json.hpp"

using namespace LCS;

void PrintSettings(const nlohmann::json& json);
std::string CreateDataFileName(std::string path, std::string file_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi);
std::string CreateFileName(std::string path, std::string file_prefix, std::string sign_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi);
std::string CreateFTLEFileName(std::string path, std::string file_prefix, std::string sign_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi);
Tensor<Vector<double,2>, 2> CalculateInterpolatedPosition(int nx, int ny, Position<double,2>& UniformGrid, double signed_calc_delta_t, std::string step_flow_maps_path, std::string file_prefix, std::string sign_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi);
int GetPrecision(double num);

int main() {  
    Clock clock;
    clock.Begin();

    const std::string proj_dir = "fast_computation";
    const std::string json_file = "sim_params.json";

    std::cout << "Number of threads: " << omp_get_max_threads() << std::endl; 

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
    const double nx = json["nx"];
    const double ny = json["ny"];
    const double data_nx = json["data_nx"];
    const double data_ny = json["data_ny"];
    const double t_min = json["t_min"];
    const double t_max = json["t_max"];
    const double data_delta_t = json["data_delta_t"];
    const int steps = json["steps"];
    const std::string file_prefix = json["file_prefix"];
    const std::string direction = json["direction"];

    const std::string positive_prefix = "positive_";
    const std::string negative_prefix = "negative_";

    const double calc_delta_t = (t_max - t_min) / steps;
    int precision = GetPrecision(data_delta_t);
    
    // Print settings
    PrintSettings(json);

    // Create paths
    std::string data_path = (std::filesystem::path(proj_dir) / "data").string();
    std::string step_flow_maps_path = (std::filesystem::path(proj_dir) / "step_flow_maps").string();
    std::string results_ftle_path = (std::filesystem::path(proj_dir) / "results" / "ftle").string();

    // Initialize the discrete flow field
    DiscreteFlowField<double,2> DiFlFi(nx,ny,data_nx,data_ny);
    DiFlFi.DataPosition().SetAll(x_min,x_max,y_min,y_max);
    std::string prefix = (std::filesystem::path(data_path) / file_prefix).string();
    DiFlFi.SetVelocityFileNamePrefix(prefix);
    DiFlFi.SetDataTimeRange(t_min,t_max);
    DiFlFi.SetDataDelta(data_delta_t);
    DiFlFi.SetDelta(calc_delta_t);
    DiFlFi.SetPrecision(precision);

    // Initialize variables whose values depend on the direction
    double t_initial;
    double t_final;
    double signed_calc_delta_t;
    std::string sign_prefix; 
    if (direction == "forward") {
        DiFlFi.SetDirection(Direction::Forward);
        t_initial = t_min;
        t_final = t_max;
        signed_calc_delta_t = calc_delta_t;
        sign_prefix = positive_prefix;
        std::cout << "*** FORWARD FTLE CALCULATION BEGINS ***" << std::endl << std::endl;

    } else if (direction == "backward") {
        DiFlFi.SetDirection(Direction::Backward);
        t_initial = t_max;
        t_final = t_min;
        signed_calc_delta_t = -calc_delta_t;
        sign_prefix = negative_prefix;
        std::cout << "*** BACKWARD FTLE CALCULATION BEGINS ***" << std::endl << std::endl;
    } else {
        std::cout << "Direction must be set either to 'forward' or 'backward'" << std::endl;
        return 1;
    }

    Clock clock_step_flowmaps;
    clock_step_flowmaps.Begin();

    // Initialize the Discrete Flow Field (DiFlFi)
    DiFlFi.SetInitialTime(t_initial);
    DiFlFi.InitialPosition().SetAll(x_min,x_max,y_min,y_max);
    std::string s_flow_maps = CreateFileName(step_flow_maps_path,file_prefix,sign_prefix,precision,DiFlFi);
    DiFlFi.InitialPosition().WriteToMemoryMappedFile(s_flow_maps);
    DiFlFi.SetStep(1);
    DiFlFi.CopyInitialPositionToCurrentPosition();
    std::cout << "*** Start step flow maps calculations from t = " << t_initial << " ***" << std::endl << std::endl;
    // Calculate the step flow maps
    for (int i = 0; i < steps; i++) {
        // Run for one timestep at a time, always starting from a uniform grid
        std::cout << "[" << i+1 << "/" << steps << "] Calculate step flow map from t = " << std::fixed << std::setprecision(precision) << DiFlFi.GetTime() << " to t = " << DiFlFi.GetTime()+signed_calc_delta_t << std::endl;
        DiFlFi.Run();
        
        std::string s_flow_maps = CreateFileName(step_flow_maps_path,file_prefix,sign_prefix,precision,DiFlFi);
        DiFlFi.CurrentPosition().WriteToMemoryMappedFile(s_flow_maps);
        DiFlFi.CurrentPosition().SetAll(x_min,x_max,y_min,y_max);
        std::cout << std::endl;
    }
    clock_step_flowmaps.End();
    std::cout << "*** Initial for loop ended successfully (Total calculation time for " << steps << " step flow maps: " << std::setprecision(4) << clock_step_flowmaps.GetTotalElapsedTime() << " s) ***" << std::endl << std::endl;

    Clock clock_fast_calculations;
    clock_fast_calculations.Begin();

    Position<double,2> UniformGrid(nx,ny);
    UniformGrid.SetAll(x_min,x_max,y_min,y_max);
    // Determine the fluid particle trajectories using interpolation
    for (int i = 0; i < steps; i++) {
        DiFlFi.InitialPosition().SetAll(x_min,x_max,y_min,y_max);
        DiFlFi.SetInitialTime(t_final - signed_calc_delta_t * (i+1));
        std::cout << "[" << i+1 << "/" << steps << "] Fast calculation with interpolation from t = " << std::fixed << std::setprecision(precision) << DiFlFi.GetTime() << " to t = " << t_final << std::endl;
        DiFlFi.CopyInitialPositionToCurrentPosition();
        
        for (int ii = 0; ii <= i; ii++) {
            auto new_position = CalculateInterpolatedPosition(nx,ny,UniformGrid,signed_calc_delta_t,step_flow_maps_path,file_prefix,sign_prefix,precision,DiFlFi);
            DiFlFi.CurrentPosition().SetAll(new_position);
            DiFlFi.CurrentPosition().UpdateOutOfBoundTensor(); 
            DiFlFi.UpdateTime();  
        }
        // Compute the FTLE field between the final and the curent positions
        std::cout << "Calculate FTLE field at t = " << DiFlFi.InitialPosition().GetTime() << std::endl;
        FTLE<double,2> ftle(DiFlFi);
        ftle.Calculate();
        std::string s_ftle = CreateFTLEFileName(results_ftle_path,file_prefix,sign_prefix,precision,DiFlFi);
        ftle.WriteToFile(s_ftle);
        std::cout << "Timestep finished successfully!" << std::endl << std::endl;
    }
    
    clock_fast_calculations.End();
    std::cout << "*** Fast calculation ended successfully ***" << std::endl << std::endl;
    std::cout << "Calculation time for " << steps << " step flow maps: " << std::setprecision(4) << clock_step_flowmaps.GetTotalElapsedTime() << " s" << std::endl;
    std::cout << "Calculation time for " << steps << " steps: " << std::setprecision(4) << clock_fast_calculations.GetTotalElapsedTime() << " s" << std::endl;
    clock.End();
    std::cout << "TOTAL CALCULATION TIME: " << std::setprecision(4) << clock.GetTotalElapsedTime() << " s" << std::endl;

    return 0;
}

void PrintSettings(const nlohmann::json& json) {
    std::cout << "*** Settings ***" << std::endl;
    std::cout << "x_min = " << json["x_min"] << std::endl;
    std::cout << "x_max = " << json["x_max"] << std::endl;
    std::cout << "y_min = " << json["y_min"] << std::endl;
    std::cout << "y_max = " << json["y_max"] << std::endl;
    std::cout << "nx = " << json["nx"] << std::endl;
    std::cout << "ny = " << json["ny"] << std::endl;
    std::cout << "data_nx = " << json["data_nx"] << std::endl;
    std::cout << "data_ny = " << json["data_ny"] << std::endl;
    std::cout << "t_min = " << json["t_min"] << std::endl;
    std::cout << "t_max = " << json["t_max"] << std::endl;
    std::cout << "data_delta_t = " << json["data_delta_t"] << std::endl;
    std::cout << "steps = " << json["steps"] << std::endl;
    std::cout << "file_prefix = " << json["file_prefix"] << std::endl;
    std::cout << "direction = " << json["direction"] << std::endl << std::endl;
}

std::string CreateDataFileName(std::string path, std::string file_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi) {
    std::stringstream ss;
    ss << file_prefix << std::fixed << std::setprecision(precision) << DiFlFi.GetTime() << ".txt";
    std::string s = (std::filesystem::path(path) / ss.str()).string();
    return s;
}

std::string CreateFileName(std::string path, std::string file_prefix, std::string sign_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi) {
    std::stringstream ss;
    ss << file_prefix << sign_prefix << std::fixed << std::setprecision(precision) << DiFlFi.GetTime() << ".bin";
    std::string s = (std::filesystem::path(path) / ss.str()).string();
    return s;
}

std::string CreateFTLEFileName(std::string path, std::string file_prefix, std::string sign_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi) {
    std::stringstream ss;
    auto direction = DiFlFi.GetDirection();
    if (direction == Forward) {
        ss << file_prefix << sign_prefix << std::fixed << std::setprecision(precision) << DiFlFi.InitialPosition().GetTime() << "-" << DiFlFi.GetTime() << ".txt";
    } else if (direction == Backward) {
        ss << file_prefix << sign_prefix << std::fixed << std::setprecision(precision) << DiFlFi.GetTime() << "-" << DiFlFi.InitialPosition().GetTime() << ".txt";
    }
    std::string s = (std::filesystem::path(path) / ss.str()).string();
    return s;
}

Tensor<Vector<double,2>, 2> CalculateInterpolatedPosition(int nx, int ny, Position<double,2>& UniformGrid, double signed_calc_delta_t, std::string step_flow_maps_path, std::string file_prefix, std::string sign_prefix, int precision, DiscreteFlowField<double, 2>& DiFlFi) {
    Velocity<double,2> StepPosition(nx,ny,UniformGrid);
    std::string s_flow_maps = CreateFileName(step_flow_maps_path,file_prefix,sign_prefix,precision,DiFlFi);
    StepPosition.ReadFromMemoryMappedFile(s_flow_maps);

    Velocity<double,2> InterpPosition(nx,ny,DiFlFi.CurrentPosition());
    InterpPosition.InterpolateFrom(StepPosition);
    auto new_position = InterpPosition.GetAll();
    return new_position;
}

int GetPrecision(double num) {
    int precision = 0;
    while (std::abs(num - static_cast<long long>(num)) > 0) {
        num *= 10;
        precision++;
    }
    return precision;
}
    
