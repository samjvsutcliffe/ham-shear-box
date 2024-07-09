import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import re

from scipy import integrate

def extract_vals(f):
    output,refine,load = f.split("-")
    return float(refine),float(load)

regex = re.compile(r'^output-.*')
folders = list(filter(regex.search,os.listdir("./")))


plt.figure()
for i in folders:
    print("loading folder: ",i)
    mpm = pd.read_csv("./{}/disp.csv".format(i))
    plt.plot(1e3*mpm["disp"].values,(1e-3/0.06)*mpm["load"].values,label=i,marker=".")
plt.xlabel("Displacement (mm)")
plt.ylabel("Shear stress (kN/m^2)")
plt.legend()

surcharge = []
peak = []
residual = []
plt.figure()
for f in folders:
    refine,load = extract_vals(f)
    mpm = pd.read_csv("./{}/disp.csv".format(f))
    width = 0.06
    p = mpm["load"].max()/width
    r = mpm["load"].values[-1]/width
    surcharge.append(load)
    peak.append(p)
    residual.append(r)
    # plt.scatter(load,r)
    # plt.scatter(load,p)
peak = [x for y, x in sorted(zip(surcharge, peak))]
residual = [x for y, x in sorted(zip(surcharge, residual))]
surcharge = sorted(surcharge)
peak = np.array(peak)
residual = np.array(residual)
surcharge = np.array(surcharge)

plt.scatter(surcharge,peak,label="Peak")
plt.scatter(surcharge,residual,label="Residual")
m,b = np.polyfit(surcharge, peak, 1)
print(m)
plt.axline((0,b),slope=m,label="Peak, {:.2f}, {:.2f}kN".format(np.arctan(m)*180/np.pi,b*1e-3))
m,b = np.polyfit(surcharge, residual, 1)
plt.axline((0,b),slope=m,label="Residual, {:.2f}, {:.2f}kN".format(np.arctan(m)*180/np.pi,b*1e-3))



plt.plot(surcharge,peak)
plt.plot(surcharge,residual)
plt.axline((0,0),slope=np.sin(30 * np.pi/180),ls="--")
plt.axline((0,131e3),slope=np.sin(42 * np.pi/180),ls="--")
plt.xlim([0,500e3])
plt.ylim([0,500e3])
plt.xlabel("Surcharge kPa")
plt.xlabel("Shear stress kPa")
plt.legend()
plt.show()
