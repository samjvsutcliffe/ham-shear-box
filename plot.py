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

print(unique_ids)


for unique_id in unique_ids:
    plt.figure(1)
    unreg = re.compile(r'^output-{}-.*'.format(unique_id))
    folders = list(filter(unreg.search,os.listdir("./")))
    for i in folders:
        print("loading folder: ",i)
        mpm = pd.read_csv("./{}/disp.csv".format(i))
        plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
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
            width = 0.06
            p = mpm["load"].max()/width
            r = mpm["load"].values[-1]/width
            surcharge.append(load)
            peak.append(p)
            residual.append(r)
            # plt.scatter(load,r)
            # plt.scatter(load,p)

    if len(peak) > 0:
        peak = [x for y, x in sorted(zip(surcharge, peak))]
        residual = [x for y, x in sorted(zip(surcharge, residual))]
        surcharge = sorted(surcharge)
        
        plt.scatter(surcharge,peak,label="Peak - {}".format(unique_id))
        plt.scatter(surcharge,residual,label="Residual - {}".format(unique_id))

        p = plt.plot(surcharge,peak)
        r = plt.plot(surcharge,residual)

        m,b = np.polyfit(surcharge, peak, 1)
        print(m)
        plt.axline((0,b),slope=m,label="Peak, {:.2f}, {:.2f}kN".format(np.arctan(m)*180/np.pi,b*1e-3),c=p[0].get_color())
        m,b = np.polyfit(surcharge, residual, 1)
        plt.axline((0,b),slope=m,label="Residual, {:.2f}, {:.2f}kN".format(np.arctan(m)*180/np.pi,b*1e-3),c=r[0].get_color())
        
        
        
        plt.axline((0,0),slope=np.tan(30 * np.pi/180),ls="--")
        plt.axline((0,131e3),slope=np.tan(42 * np.pi/180),ls="--")
        plt.xlim([0,500e3])
        plt.ylim([0,500e3])
        plt.legend()
        plt.savefig("frictional.pdf")
plt.show()
