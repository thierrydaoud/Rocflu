import glob, json, os, re, sys
import numpy as np
sys.path.insert(0, "/blue/bala1s/t.daoud/code_versions/2026/Rocflu_2d/tests")
from rocpicl_tests import outputs
from concurrent.futures import ProcessPoolExecutor
BIN = 0.65 / 1499.0
RUNS = [("AUSM+ 3D 1x","flux4_d3",1),("AUSM+ 2D 1x","flux4_d2",1),("AUSM+ 2D 2x","ref2x_flux4",2),("AUSM+ 2D 4x","ref4x_flux4",4),
        ("AUSM+UP 3D 1x","d3",1),("AUSM+UP 2D 1x","d2",1),("AUSM+UP 2D 2x","ref2x_flux5",2),("AUSM+UP 2D 4x","ref4x_flux5",4)]
def times(rd):
    out = []
    for f in glob.glob(rd + "/post/p0/*.vtu"):
        with open(f, "rb") as h:
            head = h.read(4000).decode(errors="ignore")
        out.append((float(re.search(r"Name=.TIME.[^>]*>\s*([-+.\dEe]+)", head).group(1)), int(os.path.basename(f)[:-4])))
    return out
def one(a):
    name, d, k, T = a
    t, s = min(times(d), key=lambda p: abs(p[0] - T))
    C, P, V, R = [], [], [], []
    for pdir in glob.glob(d + "/post/p[0-9]*"):
        f = "%s/%08d.vtu" % (pdir, s)
        if not os.path.exists(f): continue
        _, _, _, arr = outputs.read_post_vtu(f)
        C.append(np.array(outputs.read_post_vtu(f, centroids=True))); P.append(np.array(arr["Pressure"])); V.append(np.array(arr["Velocity"])); R.append(np.array(arr["Gas Density"]))
    C, P, V, R = map(np.concatenate, (C, P, V, R))
    r = np.hypot(C[:,0], C[:,1]); th = np.arctan2(C[:,1], C[:,0])
    w = BIN / k
    b = ((r - 0.1) / w).astype(int)
    n = b.max() + 1
    cnt = np.bincount(b, minlength=n); ok = cnt > 0
    pm = np.bincount(b, P, n)[ok] / cnt[ok]; rm = np.bincount(b, r, n)[ok] / cnt[ok]; dm = np.bincount(b, R, n)[ok] / cnt[ok]
    f2 = rm[np.where(pm > 2 * 101325.0)[0].max()]
    near = rm > f2 - 0.03
    pk = pm[near].max()
    idx = np.where(near & (pm > 0.5 * pk))[0].max()
    fh = rm[idx] + (pm[idx] - 0.5 * pk) / (pm[idx] - pm[idx + 1]) * (rm[idx + 1] - rm[idx])
    ut = -V[:,0] * np.sin(th) + V[:,1] * np.cos(th)
    m = (r > f2 - 0.05) & (r < f2 - 0.01)
    return dict(name=name, T=T*1e6, t=t, f2=f2, fh=fh, pk=pk, ut=float(np.sqrt((ut[m]**2).mean())), w=float(np.sqrt((V[m,2]**2).mean())), ps=float(P[m].std()/P[m].mean()), r=rm.tolist(), p=pm.tolist(), rho=dm.tolist())
jobs = [(n, d, k, T) for n, d, k in RUNS for T in (5e-6, 10e-6, 20e-6)]
with ProcessPoolExecutor(12) as ex:
    res = list(ex.map(one, jobs))
json.dump(res, open("ref_compare.json", "w"))
for d in res:
    print("%-15s %2.0f us t %.4e  front 2amb %.4f  half-peak %.4f  peak p %.3e | ut rms %6.1f  w rms %5.1f  p scatter %5.2f%%" % (d["name"], d["T"], d["t"], d["f2"], d["fh"], d["pk"], d["ut"], d["w"], 100*d["ps"]))
