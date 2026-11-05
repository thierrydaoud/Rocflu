# Remap the wedge grid to a planar tube: r -> x, theta -> y, optional shear of x with y.
import math, sys
src, dst, shear_deg = sys.argv[1], sys.argv[2], float(sys.argv[3])
H = 6.5e-3
L = open(src).read().split("\n")
n = int(L[1])
vals, i = [], 2
while len(vals) < 3 * n:
    l = L[i]; vals += [l[k:k+16] for k in range(0, len(l), 16) if l[k:k+16].strip()]; i += 1
assert len(vals) == 3 * n
X = list(map(float, vals[:n])); Y = list(map(float, vals[n:2*n])); Z = vals[2*n:]
th = [math.atan2(y, x) for x, y in zip(X, Y)]
tmax = max(abs(t) for t in th)
s = math.tan(math.radians(shear_deg))
ny = [t / tmax * H / 2 for t in th]
nx = [math.hypot(x, y) + s * v for x, y, v in zip(X, Y, ny)]
out = [L[0], L[1]]
for block in (nx, ny, list(map(float, Z))):
    for k in range(0, n, 5):
        out.append("".join("%16.9E" % v for v in block[k:k+5]))
out += L[i:]
open(dst, "w").write("\n".join(out))
print(dst, "n", n, "theta max deg %.4f" % math.degrees(tmax), "x %.5f %.5f" % (min(nx), max(nx)), "y %.5f %.5f" % (min(ny), max(ny)), "z", min(map(float,Z)), max(map(float,Z)))
