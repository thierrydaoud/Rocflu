import glob, json, math, os, re, sys
sys.path.insert(0, "/blue/bala1s/t.daoud/code_versions/2026/Rocflu_2d/tests")
from rocpicl_tests import outputs
from concurrent.futures import ProcessPoolExecutor
W = "/blue/bala1s/t.daoud/JoshWedge/FluidParticles2D/"
RUNS = {"3D": "/blue/bala1s/t.daoud/scratch_claude/bin_wedge", "2D100s1": W + "dz100um", "2D100s2": W + "dz100um_s2", "2D125s1": W + "dz125um"}
TIMES = [5e-6, 10e-6, 20e-6]
BIN = 0.65 / 1499.0
def times(rd):
    out = []
    for f in glob.glob(rd + "/post/p0/*.vtu"):
        with open(f, "rb") as h:
            head = h.read(4000).decode(errors="ignore")
        out.append((float(re.search(r'Name="TIME"[^>]*>\s*([-+.\dEe]+)', head).group(1)), int(os.path.basename(f)[:-4])))
    return out
def one(a):
    name, rd, T = a
    t, s = min(times(rd), key=lambda p: abs(p[0] - T))
    B = {}
    for pdir in glob.glob(rd + "/post/p[0-9]*"):
        f = "%s/%08d.vtu" % (pdir, s)
        if not os.path.exists(f): continue
        _, _, _, arr = outputs.read_post_vtu(f)
        cen = outputs.read_post_vtu(f, centroids=True)
        for k, c in enumerate(cen):
            r = math.hypot(c[0], c[1]); b = int((r - 0.1) / BIN)
            e = B.setdefault(b, [0, 0.0, 0.0, 0.0, 0.0, 0.0])
            v = arr["Velocity"][k] if "Velocity" in arr else None
            e[0] += 1; e[1] += r; e[2] += arr["Gas Density"][k]; e[3] += arr["Pressure"][k]; e[4] += (v[0] * c[0] + v[1] * c[1]) / r
            e[5] += arr["Particle Volume Fraction"][k]
    d = outputs.particle_dumps(rd)
    p = outputs.read_particle_dump(d[s])
    pr = [math.hypot(y[0], y[1]) for y in p["y"]]
    idx = sorted(B)
    return name, "%.0f" % (T * 1e6), {"t": t, "step": s, "r": [B[b][1] / B[b][0] for b in idx], "rho": [B[b][2] / B[b][0] for b in idx], "p": [B[b][3] / B[b][0] for b in idx], "ur": [B[b][4] / B[b][0] for b in idx], "phi": [B[b][5] / B[b][0] for b in idx], "pr": pr}
jobs = [(n, rd, T) for n, rd in RUNS.items() for T in TIMES]
out = {}
with ProcessPoolExecutor(6) as ex:
    for n, T, d in ex.map(one, jobs):
        out.setdefault(n, {})[T] = d
json.dump(out, open("/blue/bala1s/t.daoud/scratch_claude/ndim2_baseline/profiles.json", "w"))
print("DONE")
