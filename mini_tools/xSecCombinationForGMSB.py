import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

xSecC1C1 = pickle.load(open("../plotter/data/xSec_SMS_C1C1_13TeV.pkl", "rb"))
xSecC1N2 = pickle.load(open("../plotter/data/xSec_SMS_N2C1_13TeV.pkl", "rb"))
xSecComb = pickle.load(open("../plotter/data/xSec_SMS_TChiNG_13TeV.pkl", "rb"))
xSecGMSB = pickle.load(open("../plotter/data/xSec_GMSB.pkl", "rb"))

cx = []
cy = []
nx = []
ny = []
cox = []
coy = []
gx = []
gy = []

coyUp = []
coyDn = []


for key in xSecC1C1:
    # print key, xSecC1C1[key]
    cx.append(key)
    cy.append(xSecC1C1[key])
for key in xSecC1N2:
    # print key, xSecC1N2[key]
    nx.append(key)
    ny.append(xSecC1N2[key])
for key in xSecComb:
    # print key, xSecComb[key]
    # print key, xSecComb[key][1] * xSecComb[key][0]
    coyUp.append(xSecComb[key][0] + xSecComb[key][1] * xSecComb[key][0])
    coyDn.append(xSecComb[key][0] - xSecComb[key][1] * xSecComb[key][0])
    cox.append(key)
    coy.append(xSecComb[key])

arGMSB = {}

for key in xSecGMSB:
    temp = []
    # arGMSB[key] = {}
    for key2 in xSecGMSB[key]:
        temp.append(xSecGMSB[key][key2][0])
        # print key, key2, xSecGMSB[key][key2]
    # print np.mean(temp)
    arGMSB[key] = (np.mean(temp), 0.)

for key in arGMSB:
    # print key, arGMSB[key]
    gx.append(key)
    gy.append(arGMSB[key])

# cxSort = [_ for _, x in sorted(zip(cy, cx))]
cxSort = sorted(cx)
cySort = [x[0] for _, x in sorted(zip(cx, cy))]
nxSort = sorted(nx)
nySort = [x[0] for _, x in sorted(zip(nx, ny))]
coxSort = sorted(cox)
coySort = [x[0] for _, x in sorted(zip(cox, coy))]
coyUpSort = [x for _, x in sorted(zip(cox, coyUp))]
coyDnSort = [x for _, x in sorted(zip(cox, coyDn))]
# print len(gx), len(gy)

gxSort = sorted(gx)
gySort = [x[0] for _, x in sorted(zip(gx, gy))]
# print len(gxSort), len(gySort)
# print gxSort
# print gySort


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


coxSort = np.array(coxSort)
coySort = np.array(coySort)
coyDnSort = np.array(coyDnSort)
coyUpSort = np.array(coyUpSort)
gxSort = np.array(gxSort)
gySort = np.array(gySort)
p0 = [500., 0.01, 0.]

# popt, pcov = curve_fit(func, coxSort, coySort, p0)
# print popt
# print len(coxSort)

# z = np.polyfit(coxSort, coySort, 15)
# p = np.poly1d(z)

f = interp1d(coxSort, coySort, kind='cubic')
fUp = interp1d(coxSort, coyUpSort, kind='cubic')
fDn = interp1d(coxSort, coyDnSort, kind='cubic')

plt.figure()
# plt.plot(cxSort, cySort, label="C1C1")
# plt.plot(nxSort, nySort, label="C1N2")
plt.plot(coxSort, coySort, "-o", label="C1C1+C1N2")
plt.plot(coxSort, coyUpSort, "o-", label="C1C1+C1N2 +")
plt.plot(coxSort, coyDnSort, "o-", label="C1C1+C1N2 +")
plt.plot(gxSort, gySort, "o", label="GMSB")
# plt.plot(gxSort, func(gxSort, *popt), label="fit")
# plt.plot(gxSort, p(gxSort), label="fit")
plt.plot(gxSort, f(gxSort), label="fit")
plt.plot(gxSort, fUp(gxSort), label="fit+")
plt.plot(gxSort, fDn(gxSort), label="fit-")
plt.yscale('log')
plt.grid()
plt.legend(loc="best")
# plt.show()

erUp = []
erDn = []
for x in gxSort:
    erUp.append((fUp(x).item() - f(x).item()) / f(x).item() * 100.)
    erDn.append((fDn(x).item() - f(x).item()) / f(x).item() * 100.)
    # erDn.append(fDn(x).item())
# print np.round(erUp, 3)
# print len(erUp), len(gxSort)
# print erDn
toSave = {}
for i in range(len(erUp)):
    toSave[gxSort[i]] = erUp[i]

# print toSave
output = open("../plotter/data/xSec_GMSB_uncert.pkl", 'wb')
pickle.dump(toSave, output)
