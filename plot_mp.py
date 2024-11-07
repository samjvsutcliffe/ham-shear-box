import matplotlib as mpl
data_name = "plastic_strain"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
import matplotlib.ticker as plticker
import re
import os
import json
import numpy as np
import pandas as pd
import json
import sys
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib import cm
from multiprocessing import Pool




def get_data(filename):
    global data_name
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    vtk_points = data.GetPoints()
    xyz3d = vtk_to_numpy( vtk_points.GetData() )
    xy = xyz3d[:,0:2]
    scalar_names = [reader.GetScalarsNameInFile(i) for i in range(0, reader.GetNumberOfScalarsInFile())]
    scalar_data = data.GetPointData()
    #scalar_names = scalar_data.GetArrayNames()
    def GetScalar(scalar_name):
        return vtk_to_numpy(scalar_data.GetArray(scalar_names.index(scalar_name)))
    lx = GetScalar("size_x")
    ly = GetScalar("size_y")
    damage = GetScalar(data_name)
    #damage = GetScalar("plastic_strain")
    #damage = GetScalar("damage")
    return pd.DataFrame({"coord_x":xy[:,0], "coord_y":xy[:,1],"lx":lx,"ly":ly,
                         "damage":damage,
                         "plastic":GetScalar("plastic_strain"),
                         "sig_xx":GetScalar("sig_xx"),
                         "sig_yy":GetScalar("sig_yy"),
                         "sig_xy":GetScalar("sig_xy")
                         })

def get_data_all(folder,frame_number):
    regex = re.compile(r'sim(_\d+)?_{}'.format(frame_number))
    files = list(filter(regex.search,os.listdir(folder)))
    subframes = [get_data(folder + "/" + f) for f in files]
    df = pd.concat(subframes)
    return df



import subprocess

plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
width = 3.487
height = width / 1.618


output_regex = re.compile("output-\d+")
output_list = list(filter(output_regex.match,os.listdir()))
#output_list.sort(key=lambda x: float(x.split("-")[1]))
print(output_list)
output_dir = "./{}/".format(output_list[int(input())])

# output_dir = "./ham-shear-box/output-8-1.0e+5/"
xlim = [0.03,0.06*3]
ylim = [0,0.10]
#with open(output_dir+"settings.json") as f:
#    json_settings = json.load(f)
#    #print("Water level:{}".format(json_settings["OCEAN-HEIGHT"]))
#    #print("Domain size:{}".format(json_settings["DOMAIN-SIZE"]))
#    #water_height = json_settings["OCEAN-HEIGHT"]
#    water_height = 0
#    xlim = [0,json_settings["DOMAIN-SIZE"][0]]
#    ylim = [0,json_settings["DOMAIN-SIZE"][1]]
#    #ylim[0] = 20


ice_height = 200

plt.close("all")
files = os.listdir(output_dir)
finalcsv = re.compile("sim(_0+)?_\d+.vtk")
files_csvs = list(filter(finalcsv.match,files))

finalpbs = re.compile("sim_pb(_0+)?_\d+.json")
files_pbs = list(filter(finalpbs.match,files))
print(files_pbs)

framenumber_regex = re.compile("\d+")
files_csvs = list(map(lambda x: framenumber_regex.findall(x)[-1],files_csvs))
files_csvs.sort(key=float)
#files_csvs = list(map(lambda x: "sim_{}.vtk".format(x), files_csvs))


print(files_pbs)
files_pbs = list(map(lambda x: framenumber_regex.findall(x)[-1],files_pbs))
files_pbs.sort(key=float)
files_pbs
#files_pbs = list(map(lambda x: "sim_pb_{}.json".format(x), files_pbs))

print("files: {}".format(files_csvs))
dt = 1e4/60
time = []
max_stress = []
damage = []
full_data = []

def get_plot(i,fname,selected):
    plt.clf()
    df = get_data_all(output_dir,fname)
    print("Plot frame {}".format(i),flush=True)
    ax = fig.add_subplot(111,aspect="equal")
    #df = full_data[i]
    patch_list=[]
    #patch = Rectangle(xy=(0,0) ,width=xlim[1], height=water_height,color="blue")
    #patch_sea = [patch]
    #ps = PatchCollection(patch_sea)
    #ax.add_collection(ps)

    for a_x, a_y,lx,ly,damage in zip(df["coord_x"],
                                     df["coord_y"],
                                     df["lx"],
                                     df["ly"],
                                     df["damage"]):
        patch = Rectangle(
            xy=(a_x-lx/2, a_y-ly/2) ,width=lx, height=ly)
        patch_list.append(patch)
    p = PatchCollection(patch_list, cmap=cm.jet, alpha=1)
    df["damage"] = 0
    df["damage"].values[selected] = 1
    p.set_array(df["damage"])
    #p.set_clim([0,1.0])
    ax.add_collection(p)
    fig.colorbar(p,location="bottom",label=data_name)

    # Add the grid
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    #loc = plticker.MultipleLocator(base=h)
    with open(output_dir+"settings.json") as f:
        json_settings = json.load(f)
        h = json_settings["RESOLUTION"]
    loc = plticker.FixedLocator(np.arange(0,0.06*3,h))
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    plt.grid()
    #ax.grid(which='major',linestyle='--')
    plt.title("Frame {}".format(i))
    # plt.savefig("outframes/frame_{:05}.png".format(i))
    # plt.clf()
def plot_pb(i,fname):
    with open(output_dir + "sim_pb_{}.json".format(fname)) as f:
        d = json.load(f)
        for bc in d:
            # print(bc)
            p = np.array(bc["POSITION"][:2])
            normal = np.array(bc["NORMAL"][:2])
            radius = bc["RADIUS"]
            normal[0],normal[1] = -normal[1],normal[0]
            x0 = p - (normal * radius)
            x1 = p + (normal * radius)
            plt.plot([x0[0],x1[0]],[x0[1],x1[1]],c="black")



def plot_graphs(mp_id):
    disp_x = 0
    disp_y = 0
    sig_xx = np.zeros(len(files_csvs))
    sig_yy = np.zeros(len(files_csvs))
    sig_xy = np.zeros(len(files_csvs))
    plastic_strain = np.zeros(len(files_csvs))
    for i,fname in enumerate(files_csvs):
        df = get_data_all(output_dir,fname)
        sig_xx[i] = df["sig_xx"].values[mp_id]
        sig_xy[i] = df["sig_xy"].values[mp_id]
        # plastic_strain[i] = df["damage"].values[mp_id]
        plastic_strain[i] = df["plastic"].values[mp_id]
    plt.plot(-1*sig_xy,label="sig_xy")
    plt.plot(1e7*plastic_strain,label="plastic strain")

current_frame = 0
max_frame = len(files_csvs)-1

data_name = "plastic_strain"
f_0 = files_csvs[0]
df_0 = get_data_all(output_dir,f_0)
with open(output_dir+"settings.json") as f:
    json_settings = json.load(f)
    h = json_settings["RESOLUTION"]
def pick(pick_x,pick_y):
    x_picker = pick_x
    y_picker = pick_y
    diff_x = (df_0["coord_x"].values - x_picker) ** 2
    diff_y = (df_0["coord_y"].values - y_picker) ** 2
    dist = diff_x + diff_y
    return dist.argmin()

def replot(mp_id):
    # plt.figure(1)
    get_plot(current_frame,files_csvs[current_frame],mp_id)
    plot_pb(current_frame,files_pbs[current_frame])
    fig.canvas.draw()
    plt.figure(2)
    plt.clf()
    plot_graphs(mp_id)
    plt.gcf().canvas.draw()
    plt.pause(0.1)




# def on_press(event):
#     global current_frame,data_name
#     # print('press', event.key)
#     sys.stdout.flush()
#     # print(event.key)
#     if event.key == 'p':
#         data_name = "plastic_strain"
#         replot()
#     if event.key == 'd':
#         data_name = "damage"
#         replot()
#     if event.key == 'y':
#         data_name = "damage-ybar"
#         replot()
#     if event.key == 'x':
#         data_name = "disp_x"
#         replot()
#     if event.key == 'c':
#         data_name = "disp_y"
#         replot()
#     if event.key == 'i':
#         #data_name = "sig_xy"
#         data_name = "index"
#         replot()
#     if event.key == 'right':
#         current_frame = min(current_frame + 1,max_frame)
#         replot()
#     if event.key == 'left':
#         current_frame = max(current_frame - 1,0)
#         replot()
#     if event.key == 'up':
#         current_frame = min(current_frame + 10,max_frame)
#         replot()
#     if event.key == 'down':
#         current_frame = max(current_frame - 10,0)
#         replot()


def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
    if event.xdata:
        mp_id = pick( event.xdata, event.ydata)
        print("MP ID:",mp_id)
        replot(mp_id)



print("Plotting {}, {}",current_frame,files_csvs[current_frame])

# replot()
# plt.show()


# print("Mp index: ",mp_id)
# fig = plt.figure()#plt.figure(figsize=(16,9),dpi=200)
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
# get_plot(0,files_csvs[0],mp_id)
# # plt.show()

# disp_x = 0
# disp_y = 0
# sig_xx = np.zeros(len(files_csvs))
# sig_yy = np.zeros(len(files_csvs))
# sig_xy = np.zeros(len(files_csvs))
# plastic_strain = np.zeros(len(files_csvs))
# # mp_id = 0
# for i,fname in enumerate(files_csvs):
#     df = get_data_all(output_dir,fname)
#     sig_xx[i] = df["sig_xx"].values[mp_id]
#     sig_xy[i] = df["sig_xy"].values[mp_id]
#     # plastic_strain[i] = df["damage"].values[mp_id]
#     plastic_strain[i] = df["plastic"].values[mp_id]

# fig2 = plt.figure()
# plt.plot(-1*sig_xy,label="sig_xy")
# plt.plot(1e7*plastic_strain,label="plastic strain")
# plt.show()

fig = plt.figure()#plt.figure(figsize=(16,9),dpi=200)
# fig.canvas.mpl_connect('key_press_event', on_press)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
mp_id = pick(0.06 * 1.5,(0.06 * 0.5) + (2 * h))
replot(mp_id)
plt.show()




#def wrapper(x):
#    get_plot(x[0],x[1])
#if __name__ == '__main__':
#    with Pool(8) as p:
#        p.map(wrapper, enumerate(files_csvs))
#cmd_str = "ffmpeg -y -framerate 60 -pattern_type glob -i 'outframes/*.png' -c:v libx264 -pix_fmt yuv420p out.mp4"
#subprocess.run(cmd_str, shell=True)


