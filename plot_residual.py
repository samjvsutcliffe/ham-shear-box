import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import os
import numpy as np
import re
mpl.use("pdf")

from scipy import integrate

#plt.rc('font', family='serif', serif='Times')
## plt.rc('text', usetex=True)
#plt.rc('xtick', labelsize=8)
#plt.rc('ytick', labelsize=8)
#plt.rc('axes', labelsize=8)
plt.style.use("seaborn-paper")
plt.rc('font', family='serif', serif='Times')
# plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=8)
plt.rc('axes', labelsize=8)
# mpl.rcParams.update({"pgf.preamble": "\n".join([ # plots will use this preamble
#         r"\usepackage[utf8]{inputenc}",
#         r"\usepackage[T1]{fontenc}",
#         r"\usepackage[detect-all,locale=DE]{siunitx}",
#         ])})
# width = 3.487
# height = width / 1.618

width = 3.487 * 1
height = width / 1.5

prop_cycle = plt.rcParams['axes.prop_cycle']
colours = prop_cycle.by_key()['color']

def extract_vals(f):
    output,refine,load = f.split("-")
    #refine = float(refine)
    return refine,float(load)

#top_dir = "./plastic-damage-residual/"
top_dir = "/nobackup/rmvn14/paper-1/plastic-damage-residual/"
regex = re.compile(r'^output-.*')
folders = list(filter(regex.search,os.listdir(top_dir)))

print(folders)
unique_ids = []
for f in folders:
    u = f.split("-")[1]
    if u not in unique_ids:
        unique_ids.append(u)

unique_ids.sort(key=lambda x: float(x.split("_")[-2]))
# unique_ids.sort()
print(unique_ids)
# unique_ids = ["4"]

prop_cycle = plt.rcParams['axes.prop_cycle']
colours = prop_cycle.by_key()['color']
fig = plt.figure(figsize=(width,height),dpi=200)
fig = plt.figure(figsize=(width,height),dpi=200)

load_zeroing = True
# load_zeroing = False
# load_combined = True
load_combined = False
load_clipping = False

def get_load(filename):
    mpm = pd.read_csv(top_dir+filename)
    if load_clipping:
        mpm = mpm[mpm["disp"] >= 0.01e-3]
    
    if len(mpm["load"]) > 0:
        mpm = mpm[mpm["disp"]<3e-3]
        # mpm["load-diff"] = mpm["l-left"] + mpm["l-right"]
        if load_combined:
            mpm["load"] = mpm["load-diff"]
        if load_zeroing:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
        mpm["stress"] = mpm["load"] / (0.06 - mpm["disp"])
    return mpm

E = 1e9
nu = 0.24
G = (E / ( 2 * (1 - nu)))
print("G actual: {}GPa".format(1e-9*G))
for colour,unique_id in zip(colours,unique_ids):
    plt.figure(1)
    unreg = re.compile(r'^output-{}-.*'.format(unique_id))
    folders = list(filter(unreg.search,os.listdir(top_dir)))
    folders.sort(key=lambda x: float(x.split("-")[2]))
    folders_filtered = list(filter(lambda x: x.split("-")[2] == "100000.0",folders))
    for i in folders_filtered:
        print("loading folder: ",i)
        mpm = get_load("./{}/disp.csv".format(i))
        if len(mpm["load"]) > 0:
            #if load_zeroing:
            #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
            l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label="d = {}".format(i.split("_")[-2]),marker=".")
            # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load-diff"].values,label=i,marker="x",c=l[0].get_color())
            print("Shear modulus {}GPa".format(1e-9*mpm["load"].max()/mpm["disp"].values[mpm["load"].argmax()]))
            maxload = (1e-3/0.06)*mpm["load"].max()
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Shear stress (kPa)")
    plt.legend()
    # plt.legend(fontsize="5")
    plt.tight_layout()
    plt.savefig("load-disp.pdf")

    # plt.figure()
    # plt.title(unique_id)
    # for i in folders:
    #     print("loading folder: ",i)
    #     mpm = get_load("./{}/disp.csv".format(i))
    #     if len(mpm["load"]) > 0:
    #         #if load_zeroing:
    #         #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
    #         l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
    #         # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["l-left"].values,label=i,marker=".")
    #         # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["l-right"].values,label=i,marker=".")
    #         maxload = 200
    #         maxp=mpm["plastic"].max()
    #         maxd=mpm["damage"].max()
    #         maxp=0.1e0
    #         maxd=1e3
    #         plt.plot(1e3*mpm["disp"].values,maxload*mpm["plastic"].values/maxp,label="",marker="x",ls="--",c=l[0].get_color())
    #         plt.plot(1e3*mpm["disp"].values,maxload*mpm["damage"].values/maxd,label="",marker="o",ls="--",c=l[0].get_color())
    # plt.xlabel("Displacement (mm)")
    # plt.ylabel("Load (N)")
    # plt.legend()
    # plt.figure()
    # for i in folders:
    #     print("loading folder: ",i)
    #     mpm = pd.read_csv("./{}/disp.csv".format(i))
    #     plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
    # plt.xlabel("Displacement (mm)")
    # plt.ylabel("Load (N)")
    # plt.legend()
    # plt.savefig("load-disp-{}.pdf".format(unique_id))

    surcharge = []
    peak = []
    residual = []
    plt.figure(2)
    for f in folders:
        refine,load = extract_vals(f)
        mpm = get_load("./{}/disp.csv".format(f))
        # mpm["load"] = mpm["l-left"]
        if len(mpm["load"]) > 0:
            #if load_zeroing:
            #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
            width = 0.06
            p = mpm["load"].max()/width
            r = mpm["load"].values[-1]/width
            #residual_window = 0.25
            #residual_back = round(len(mpm["load"].values) * (1 - residual_window))
            #r = mpm["load"].values[residual_back:].mean()/width
            surcharge.append(load)
            peak.append(p)
            residual.append(r)
            # plt.scatter(load,r)
            # plt.scatter(load,p)

    if len(peak) > 0:
        peak = [x for y, x in sorted(zip(surcharge, peak))]
        residual = [x for y, x in sorted(zip(surcharge, residual))]
        surcharge = sorted(surcharge)

        unique_id = "d = {}".format(unique_id.split("_")[-2])
        m,b = np.polyfit(surcharge, peak, 1)
        #unique_id = "D_res = "+unique_id.split("_")[-1]
        #plt.scatter(surcharge,peak,label="Peak - {} - {:.2f}, {:.2f}kN".format(unique_id,np.arctan(m)*180/np.pi,b*1e-3),color=colour)
        #p = plt.plot(surcharge,peak,color=colour)
        #plt.axline((0,b),slope=m,c=p[0].get_color())
        m,b = np.polyfit(surcharge, residual, 1)
        #plt.scatter(surcharge,residual,label="Residual - {} - $\SI{{{:.2f}}}{{\degree}}$, $\SI{{{:.2f}}}{{\kilo\pascal}}$".format(unique_id,np.arctan(m)*180/np.pi,b*1e-3),color=colour,marker="x")
        plt.scatter(surcharge,residual,label="Residual - {} - {:.2f} deg, {:.2f} kPa".format(unique_id,np.arctan(m)*180/np.pi,b*1e-3),color=colour,marker="x") # r = plt.plot(surcharge,residual,color=colour,ls="--")
        plt.axline((0,b),slope=m,c=colour,ls="--")

        plt.axline((0,0),slope=np.tan(30 * np.pi/180),ls="-")
        plt.axline((0,131e3),slope=np.tan(42 * np.pi/180),ls="-")
        ticformat = ticker.FuncFormatter(lambda x,pos: "{0:g}".format(x*1e-3))
        plt.gca().xaxis.set_major_formatter(ticformat)
        plt.gca().yaxis.set_major_formatter(ticformat)
        plt.xlim([0,500e3])
        plt.ylim([0,500e3])
        plt.xlabel("Normal load (kPa)")
        plt.ylabel("Shear stress (kPa)")
        # plt.legend(fontsize="5")
        # plt.legend()
        plt.legend()
        plt.tight_layout()
        plt.savefig("frictional.pdf")
plt.show()
