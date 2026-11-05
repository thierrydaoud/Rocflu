import glob, os, re, sys
import numpy as np
sys.path.insert(0, "/blue/bala1s/t.daoud/code_versions/2026/Rocflu_2d/tests")
from rocpicl_tests import outputs
from concurrent.futures import ProcessPoolExecutor
def times(rd):
    out = []
    for f in glob.glob(rd + "/post/p0/*.vtu"):
        with open(f, "rb") as h:
            head = h.read(4000).decode(errors="ignore")
        out.append((float(re.search(r"Name=.TIME.[^>]*>\s*([-+.\dEe]+)", head).group(1)), int(os.path.basename(f)[:-4])))
    return out
def load(rd, T):
    t, s = min(times(rd), key=lambda p: abs(p[0] - T))
    C, P, V, R = [], [], [], []
    for pdir in glob.glob(rd + "/post/p[0-9]*"):
        f = "%s/%08d.vtu" % (pdir, s)
        if not os.path.exists(f): continue
        _, _, _, arr = outputs.read_post_vtu(f)
        C.append(np.array(outputs.read_post_vtu(f, centroids=True))); P.append(np.array(arr["Pressure"])); V.append(np.array(arr["Velocity"])); R.append(np.array(arr["Gas Density"]))
    return t, s, *map(np.concatenate, (C, P, V, R))
def planar(a):
    rd, T = a
    t, s, C, P, V, R = load(rd, T)
    xe = 0.248 + 14207.862 * t
    # 24 rows in y, 5 layers in z where present
    yb = np.minimum(((C[:,1] + 0.00325) / 0.0065 * 24).astype(int), 23)
    fr = []
    for j in range(24):
        m = yb == j
        fr.append(C[m,0][P[m] > 2 * 101325.0].max())
    fr = np.array(fr)
    z = "zspread n/a"
    if np.ptp(C[:,2]) > 1e-5:
        zb = np.round(C[:,2] / 2e-4 - 0.5).astype(int)
        z = "front by layer mm " + " ".join("%+.2f" % (1e3*(C[zb==k,0][P[zb==k] > 2*101325.].max() - xe)) for k in range(5))
    return "%-11s %2.0f us front by y row: mean %+6.2f min %+6.2f max %+6.2f std %5.2f mm | %s" % (os.path.basename(rd), T*1e6, 1e3*(fr.mean()-xe), 1e3*(fr.min()-xe), 1e3*(fr.max()-xe), 1e3*fr.std(), z)
def wedge(a):
    name, rd, T = a
    t, s, C, P, V, R = load(rd, T)
    r = np.hypot(C[:,0], C[:,1])
    th = np.arctan2(C[:,1], C[:,0])
    tb = np.minimum(((th + np.radians(0.5)) / np.radians(1.0) * 24).astype(int), 23)
    fr = np.array([r[tb==j][P[tb==j] > 2*101325.].max() for j in range(24)])
    ut = -V[:,0]*np.sin(th) + V[:,1]*np.cos(th)
    # expected from the IC: u along x, so u_t about -u sin(theta), under 0.5 degree
    m = (r > fr.mean() - 0.05) & (r < fr.mean() - 0.01)
    return "%-9s %2.0f us front by theta row: mean %.4f min %.4f max %.4f std %5.2f mm | behind front: ut rms %7.1f w rms %6.1f p scatter %5.2f%%" % (name, T*1e6, fr.mean(), fr.min(), fr.max(), 1e3*fr.std(), np.sqrt((ut[m]**2).mean()), np.sqrt((V[m,2]**2).mean()), 100*P[m].std()/P[m].mean())
with ProcessPoolExecutor(8) as ex:
    for l in ex.map(planar, [(os.path.join(os.getcwd(), r), T) for r in ("aligned_d2","aligned_d3","sheared_d2","sheared_d3") for T in (5e-6,10e-6,20e-6)]): print(l)
    print("--- wedge gas only")
    for l in ex.map(wedge, [(n, "/blue/bala1s/t.daoud/scratch_claude/wedge_gasonly/" + n, T) for n in ("d2","d3","o1_d2","o1_d3") for T in (5e-6,10e-6,20e-6)]): print(l)
