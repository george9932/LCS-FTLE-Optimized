import numpy as np
import matplotlib.pyplot as plt
import json
import os
import timeit
from multiprocessing import Pool

def ReadJSON(proj_dir, json_file):

    # Construct the relative path to the JSON file
    json_path = os.path.join(proj_dir, json_file)
    # Open and read the JSON file
    with open(json_path, "r") as file:
        data_json = json.load(file)

    # Read settings from the JSON file
    return (
        data_json["x_min"],
        data_json["x_max"],
        data_json["y_min"],
        data_json["y_max"],
        data_json["nx"],
        data_json["ny"],
        data_json["t_min"],
        data_json["t_max"],
        data_json["data_delta_t"],
        data_json["file_prefix"],
        data_json["steps"],
        data_json["direction"]
    )

def CreateFTLEFilePath(direction, proj_dir, file_prefix, sign_prefix, t_final, t_current, precision):

    if direction == "forward":
        filepath = os.path.join(proj_dir, "results", "ftle", file_prefix + sign_prefix + f"{t_current:.{precision}f}" + "-" + f"{t_final:.{precision}f}" + ".txt")
    elif direction == "backward":
        filepath = os.path.join(proj_dir, "results", "ftle", file_prefix + sign_prefix + f"{t_final:.{precision}f}" + "-" + f"{t_current:.{precision}f}" + ".txt")
    return filepath

def GetPrecision(num):

    precision = 0
    while abs(num - int(num)) > 0:
        num *= 10
        precision += 1
    return precision

def PlotContour(args):

    X, Y, ftle, contour_levels, cmap, ftle_path_out, title = args
    plt.contourf(X, Y, ftle, contour_levels, cmap=cmap)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title(title)
    plt.colorbar()
    plt.savefig(ftle_path_out)
    plt.close()

def PlotFTLE(direction, signed_calc_delta_t, steps, proj_dir, file_prefix, sign_prefix, precision, t_final, nx, ny, X, Y, contour_levels, cmap):
    
    for i in range(steps):
        t_current = t_final - signed_calc_delta_t*(i+1)
        title = f"Time = {t_current:.{precision}f}"
        print("[" + str(i+1) + "/" + str(steps) + "] t = " + f"{t_current:.{precision}f}")

        # Construct the relative path to the FTLE files
        ftle_path_in = CreateFTLEFilePath(direction, proj_dir, file_prefix, sign_prefix, t_final, t_current, precision)
        ftle_path_out = ftle_path_in.replace("txt","png")
        
        # Generate the FTLE contour
        ftle = np.genfromtxt(ftle_path_in, skip_header=3).reshape((nx,ny))
        
        args_list = []
        args_list.append((X, Y, ftle, contour_levels, cmap, ftle_path_out,title))

        # Use multiprocessing to parallelize the plotting
        with Pool() as pool:
            pool.map(PlotContour, args_list)

def PlotMain(x_min, x_max, y_min, y_max, nx, ny, precision, t_min, t_max, steps, file_prefix, proj_dir, direction, contour_levels, cmap):

    x = np.linspace(x_min,x_max,nx)
    y = np.linspace(y_min,y_max,ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    calc_delta_t = (t_max-t_min)/steps

    # Plot the FTLE field
    if direction == "forward":
        print("Plotting Forward FTLE Field")
        sign_prefix = "positive_"
        signed_calc_delta_t = calc_delta_t
        t_final = t_max
        PlotFTLE(direction, signed_calc_delta_t, steps, proj_dir, file_prefix, sign_prefix, precision, t_final, nx, ny, X, Y, contour_levels, cmap)
        
    elif direction == "backward":
        print("Plotting Backward FTLE Field")
        sign_prefix = "negative_"
        signed_calc_delta_t = -calc_delta_t
        t_final = t_min
        PlotFTLE(direction, signed_calc_delta_t, steps, proj_dir, file_prefix, sign_prefix, precision, t_final, nx, ny, X, Y, contour_levels, cmap)
    
    else:
        print("Cannot plot if both Forward and Backward are False")

t1 = timeit.default_timer()
if __name__ == "__main__":

    proj_dir = "fast_computation"
    json_file = "sim_params.json"
    x_min, x_max, y_min, y_max, nx, ny, t_min, t_max, data_delta_t, file_prefix, steps, direction = ReadJSON(proj_dir, json_file)
    precision = GetPrecision(data_delta_t)

    cmap = plt.cm.viridis
    contour_levels = 100
    PlotMain(x_min, x_max, y_min, y_max, nx, ny, precision, t_min, t_max, steps, file_prefix, proj_dir, direction, contour_levels, cmap)
t2 = timeit.default_timer()
print(f"Total plotting time: "+f"{t2-t1:.{1}f}"+" s")