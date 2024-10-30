import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re

from scipy import integrate

def extract_vals(f):
    output,refine,load = f.split("-")
    #refine = float(refine)
    return refine,float(load)

regex = re.compile(r'^output-.*')
folders = list(filter(regex.search,os.listdir("./")))

print(folders)
unique_ids = []
for f in folders:
    u = f.split("-")[1]
    if u not in unique_ids:
        unique_ids.append(u)

#unique_ids.sort(key=float)
print(unique_ids)
# unique_ids = ["4"]

prop_cycle = plt.rcParams['axes.prop_cycle']
colours = prop_cycle.by_key()['color']
plt.figure(1)
plt.figure(2)

for colour,unique_id in zip(colours,unique_ids):
    plt.figure(1)
    unreg = re.compile(r'^output-{}-.*'.format(unique_id))
    folders = list(filter(unreg.search,os.listdir("./")))
    folders.sort(key=lambda x: float(x.split("-")[2]))
    for i in folders:
        print("loading folder: ",i)
        mpm = pd.read_csv("./{}/disp.csv".format(i))
        if len(mpm["load"]) > 0:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
            l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
            maxload = (1e-3/0.06)*mpm["load"].max()
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Load (N)")
    plt.legend()

    plt.figure()
    plt.title(unique_id)
    for i in folders:
        print("loading folder: ",i)
        mpm = pd.read_csv("./{}/disp.csv".format(i))
        if len(mpm["load"]) > 0:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
            l=plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
            maxload = 200
            maxp=mpm["plastic"].max()
            maxd=mpm["damage"].max()
            maxp=0.1e0
            maxd=1e3
            plt.plot(1e3*mpm["disp"].values,maxload*mpm["plastic"].values/maxp,label="",marker="x",ls="--",c=l[0].get_color())
            plt.plot(1e3*mpm["disp"].values,maxload*mpm["damage"].values/maxd,label="",marker="o",ls="--",c=l[0].get_color())
    plt.xlabel("Displacement (mm)")
    plt.ylabel("Load (N)")
    plt.legend()
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
        mpm = pd.read_csv("./{}/disp.csv".format(f))
        if len(mpm["load"]) > 0:
            mpm["load"] = mpm["load"] - mpm["load"].values[0]
            width = 0.06
            p = mpm["load"].max()/width
            residual_window = 0.5
            residual_back = round(len(mpm["load"].values) * (1 - residual_window))
            r = mpm["load"].values[residual_back:].mean()/width
            surcharge.append(load)
            peak.append(p)
            residual.append(r)
            # plt.scatter(load,r)
            # plt.scatter(load,p)

    if len(peak) > 0:
        peak = [x for y, x in sorted(zip(surcharge, peak))]
        residual = [x for y, x in sorted(zip(surcharge, residual))]
        surcharge = sorted(surcharge)

        m,b = np.polyfit(surcharge, peak, 1)
        #unique_id = "D_res = "+unique_id.split("_")[-1]
        plt.scatter(surcharge,peak,label="Peak - {} - {:.2f}, {:.2f}kN".format(unique_id,np.arctan(m)*180/np.pi,b*1e-3),color=colour)
        p = plt.plot(surcharge,peak,color=colour)
        plt.axline((0,b),slope=m,c=p[0].get_color())
        m,b = np.polyfit(surcharge, residual, 1)
        plt.scatter(surcharge,residual,label="Residual - {} - {:.2f}, {:.2f}kN".format(unique_id,np.arctan(m)*180/np.pi,b*1e-3),color=colour,marker="x")
        r = plt.plot(surcharge,residual,color=colour)
        plt.axline((0,b),slope=m,c=r[0].get_color())

        plt.axline((0,0),slope=np.tan(30 * np.pi/180),ls="--")
        plt.axline((0,131e3),slope=np.tan(42 * np.pi/180),ls="--")
        plt.xlim([0,500e3])
        plt.ylim([0,500e3])
        plt.xlabel("Normal load")
        plt.ylabel("Shear stress")
        plt.legend()
        plt.savefig("frictional.pdf")
plt.show()
