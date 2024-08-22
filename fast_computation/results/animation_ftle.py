import matplotlib.pyplot as plt
import matplotlib.animation as animation
import json
import os
import timeit
from PIL import Image
import glob
import sys

def ReadJSON(proj_dir, json_file):

    # Construct the relative path to the JSON file
    json_path = os.path.join(proj_dir, json_file)
    # Open and read the JSON file
    with open(json_path, "r") as file:
        data_json = json.load(file)

    if data_json["direction"] == "forward":
        sign_prefix = "positive_"
    elif data_json["direction"] == "backward":
        sign_prefix = "negative_"
    else:
        print("Error: Invalid direction value")  
        sys.exit(1) 

    # Read settings from the JSON file
    return (
        data_json["t_min"],
        data_json["t_max"],
        data_json["data_delta_t"],
        data_json["file_prefix"],
        sign_prefix,
        data_json["direction"]
    )

def init():
    img.set_data(Image.open(image_files[0]))
    return [img]

def update(frame):
    img.set_data(Image.open(image_files[frame]))
    return [img]

def GetPrecision(num):
    precision = 0
    while abs(num - int(num)) > 0:
        num *= 10
        precision += 1
    return precision


t1 = timeit.default_timer()
if __name__ == "__main__":

    proj_dir = "fast_computation"
    json_file = "sim_params.json"
    t_min, t_max, data_delta_t, file_prefix, sign_prefix, direction = ReadJSON(proj_dir, json_file)

    precision = GetPrecision(data_delta_t)

    # Load and sort all .png files in the current directory
    path_in = os.path.join(proj_dir, "results", "ftle", "*.png")
    image_files = glob.glob(path_in)
    sorted_image_files = sorted(image_files, key=lambda x: float(x.rsplit("_",1)[1].rsplit(".",1)[0].split("-")[1]))
    fps = 25

    fig, ax = plt.subplots()
    ax.axis("off")
    img = ax.imshow(Image.open(sorted_image_files[0]))

    # Create animation
    print("Preparing animation of the FTLE field")
    ani = animation.FuncAnimation(fig, update, frames=len(sorted_image_files), init_func=init, blit=True)
    path_out = path_in.replace("*.png","")
    path_out = path_out + file_prefix + sign_prefix + f"ftle_{t_min:.{precision}f}-{t_max:.{precision}f}" + ".gif"
    ani.save(path_out, fps=fps)
t2 = timeit.default_timer()
print(f"Animation saved  (Total running time: "+f"{t2-t1:.{1}f}"+" s)")