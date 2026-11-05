import glob, json, os, re, sys
import numpy as np
sys.path.insert(0, "/blue/bala1s/t.daoud/code_versions/2026/Rocflu_2d/tests")
from rocpicl_tests import outputs
from concurrent.futures import ProcessPoolExecutor
RUNS = ["aligned_d2", "aligned_d3", "sheared_d2", "sheared_d3"]
TIMES = [5e-6, 10e-6, 20e-6]
BIN = 0.65 / 1499.0
W, X0 = 14207.862, 0.248
P1, P2, R1, R2, U2 = 101325.0, 2.0265e8, 1.2047767, 7.2076398, 11832.979
def times(rd):
    out = []
    for f in glob.glob(rd + "/post/p0/*.vtu"):
        with open(f, "rb") as h:
            head = h.read(4000).decode(errors="ignore")
        out.append((float(re.search(r"Name=.TIME.[^>]*>\s*([-+.\dEe]+)", head).group(1)), int(os.path.basename(f)[:-4])))
    return out
def one(a):
    rd, T = a
    t, s = min(times(rd), key=lambda p: abs(p[0] - T))
    C, R, P, V = [], [], [], []
    for pdir in glob.glob(rd + "/post/p[0-9]*"):
        f = "%s/%08d.vtu" % (pdir, s)
        if not os.path.exists(f): continue
        _, _, _, arr = outputs.read_post_vtu(f)
        C.append(np.array(outputs.read_post_vtu(f, centroids=True)))
        R.append(np.array(arr["Gas Density"])); P.append(np.array(arr["Pressure"])); V.append(np.array(arr["Velocity"]))
    C, R, P, V = map(np.concatenate, (C, R, P, V))
    b = ((C[:,0] - 0.1) / BIN).astype(int)
    u = np.unique(b)
    xm = np.array([C[b==k,0].mean() for k in u]); pm = np.array([P[b==k].mean() for k in u]); rm = np.array([R[b==k].mean() for k in u]); um = np.array([V[b==k,0].mean() for k in u])
    xe = X0 + W * t
    # front at half the pressure jump, and the outermost bin above twice ambient
    k = np.where(pm > 0.5 * (P1 + P2))[0].max()
    xhalf = xm[k] + (pm[k] - 0.5*(P1+P2)) / (pm[k] - pm[k+1]) * (xm[k+1] - xm[k])
    x2 = xm[np.where(pm > 2 * P1)[0].max()]
    # plateau behind the shock
    m = (C[:,0] > xe - 0.05) & (C[:,0] < xe - 0.01)
    # mass-based front: exact integral of rho in x for a single shock
    return dict(run=os.path.basename(rd), T=T*1e6, t=t, step=s, xe=xe, xhalf=xhalf, x2=x2, rho=R[m].mean(), p=P[m].mean(), u=V[m,0].mean(), vrms=float(np.sqrt((V[m,1]**2).mean())), wrms=float(np.sqrt((V[m,2]**2).mean())), rhostd=R[m].std()/R[m].mean(), pstd=P[m].std()/P[m].mean(), pmax=pm.max()/P2, x=xm.tolist(), pm=pm.tolist(), rm=rm.tolist(), um=um.tolist())
with ProcessPoolExecutor(6) as ex:
    res = list(ex.map(one, [(os.path.join(os.getcwd(), r), T) for r in RUNS for T in TIMES]))
json.dump(res, open("planar.json", "w"))
print("exact: rho %.4f p %.4e u %.1f" % (R2, P2, U2))
for d in res:
    print("%-11s %2.0f us step %5d  err half %+6.2f mm  err 2amb %+6.2f mm | rho %+6.2f%% p %+6.2f%% u %+6.2f%% | scatter rho %.2f%% p %.2f%% | vrms %6.1f wrms %6.1f | pmax/p2 %.3f" % (d["run"], d["T"], d["step"], 1e3*(d["xhalf"]-d["xe"]), 1e3*(d["x2"]-d["xe"]), 100*(d["rho"]/R2-1), 100*(d["p"]/P2-1), 100*(d["u"]/U2-1), 100*d["rhostd"], 100*d["pstd"], d["vrms"], d["wrms"], d["pmax"]))
