# Exact Riemann solution for the wedge states, gamma 1.4.
import math
g=1.4
rl,ul,pl=7.2076398498353,11832.9791498099548,202650000.0
rr,ur,pr=1.2047767301700,0.0,101325.0
cl=math.sqrt(g*pl/rl); cr=math.sqrt(g*pr/rr)
def f(p,rk,pk,ck):
    if p>pk:
        A=2/((g+1)*rk);B=(g-1)/(g+1)*pk
        return (p-pk)*math.sqrt(A/(p+B))
    return 2*ck/(g-1)*((p/pk)**((g-1)/(2*g))-1)
lo,hi=1.0,1e10
for _ in range(300):
    m=math.sqrt(lo*hi)
    if f(m,rl,pl,cl)+f(m,rr,pr,cr)+(ur-ul)>0: hi=m
    else: lo=m
p=m; u=0.5*(ul+ur)+0.5*(f(p,rr,pr,cr)-f(p,rl,pl,cl))
W=ur+cr*math.sqrt((g+1)/(2*g)*p/pr+(g-1)/(2*g))
rs=rr*((p/pr+(g-1)/(g+1))/((g-1)/(g+1)*p/pr+1))
rlstar=rl*(p/pl)**(1/g) if p<pl else rl*((p/pl+(g-1)/(g+1))/((g-1)/(g+1)*p/pl+1))
print("p* %.6e u* %.3f shock speed %.3f rho*R %.5f rho*L %.5f  (p*-pL)/pL %.2e (u*-uL)/uL %.2e"%(p,u,W,rs,rlstar,(p-pl)/pl,(u-ul)/ul))
print("left state c %.1f, u+c %.1f; rarefaction head reaches shock from x=0.1 at t=%.1f us"%(cl,ul+cl,(0.248-0.1)/(ul+cl-W)*1e6))
for t in (5,10,20): print(t,"us shock x = %.5f"%(0.248+W*t*1e-6))
