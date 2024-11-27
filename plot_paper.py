import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re

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
width = 3.487
height = width / 1.618

def extract_vals(f):
    output,refine,load = f.split("-")
    #refine = float(refine)
    return refine,float(load)

regex = re.compile(r'^output-.*')
folders = list(filter(regex.search,os.listdir("./")))
folders = [
    "output-4.0_4_1.0_100.0-100000.0",
    "output-4.0_4_1.0_100.0-200000.0",
    "output-4.0_4_1.0_100.0-300000.0"
    # "output-SE_4.0_2_0.5_1.0_10000.0-50000.0" ,
    # "output-SE_4.0_2_0.5_1.0_10000.0-100000.0",
    # "output-SE_4.0_2_0.5_1.0_10000.0-150000.0",
    # "output-SE_4.0_2_0.5_1.0_10000.0-200000.0",
    # "output-SE_4.0_2_0.5_1.0_10000.0-250000.0",
    # "output-SE_4.0_2_0.5_1.0_10000.0-300000.0"
]
names = [
    "Load: 100kPa",
    "Load: 200kPa",
    "Load: 300kPa"
]
print(folders)

prop_cycle = plt.rcParams['axes.prop_cycle']
colours = prop_cycle.by_key()['color']
plt.figure(1)
plt.figure(2)

#load_zeroing = True
load_zeroing = False
load_clipping = False

def get_load(filename):
    mpm = pd.read_csv(filename)
    if load_clipping:
        mpm = mpm[mpm["disp"] >= 0.01e-3]
    if len(mpm["load"]) > 0:
        if load_zeroing:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
    return mpm

plt.figure(1)
for i in folders:
    print("loading folder: ",i)
    mpm = get_load("./{}/disp.csv".format(i))
    if len(mpm["load"]) > 0:
        #if load_zeroing:
        #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
        l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
        maxload = (1e-3/0.06)*mpm["load"].max()
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend(names)

plt.figure()
for i in folders:
    print("loading folder: ",i)
    mpm = get_load("./{}/disp.csv".format(i))
    if len(mpm["load"]) > 0:
        #if load_zeroing:
        #    mpm["load"] = mpm["load"] - mpm["load"].values[0]
        l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
        # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["l-left"].values,label=i,marker=".")
        # plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["l-right"].values,label=i,marker=".")
        maxload = 200
        maxp=mpm["plastic"].max()
        maxd=mpm["damage"].max()
        maxp=0.1e0
        maxd=1e3
        #plt.plot(1e3*mpm["disp"].values,maxload*mpm["plastic"].values/maxp,label="",marker="x",ls="--",c=l[0].get_color())
        #plt.plot(1e3*mpm["disp"].values,maxload*mpm["damage"].values/maxd,label="",marker="o",ls="--",c=l[0].get_color())
plt.xlabel("Displacement (mm)")
plt.ylabel("Load (N)")
plt.legend(names)
plt.savefig("load-disp.pdf")
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

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colour=colors[0]
    m,b = np.polyfit(surcharge, peak, 1)
    #unique_id = "D_res = "+unique_id.split("_")[-1]
    plt.scatter(surcharge,peak,label="Peak - {:.2f}, {:.2f}kN".format(np.arctan(m)*180/np.pi,b*1e-3),color=colour)
    #p = plt.plot(surcharge,peak,color=colour)
    plt.axline((0,b),slope=m,c=colour,ls="--")

    m,b = np.polyfit(surcharge, residual, 1)
    plt.scatter(surcharge,residual,label="Residual - {:.2f}, {:.2f}kN".format(np.arctan(m)*180/np.pi,b*1e-3),color=colour,marker="x")
    #r = plt.plot(surcharge,residual,color=colour,ls="--")
    plt.axline((0,b),slope=m,c=colour,ls="--")

    plt.axline((0,0),slope=np.tan(30 * np.pi/180),ls="-")
    plt.axline((0,131e3),slope=np.tan(42 * np.pi/180),ls="-")
    plt.xlim([0,500e3])
    plt.ylim([0,500e3])
    plt.xlabel("Normal load (Pa)")
    plt.ylabel("Shear stress (Pa)")
    plt.legend()
    plt.savefig("frictional.pdf")
plt.show()
