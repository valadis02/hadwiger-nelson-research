import math, sys, time, numpy as np
from upper import voronoi_cell, sublattices
from upper2 import feasible_x
ks=[int(x) for x in sys.argv[1].split(",")]; t0=time.time()
best=None
for k in ks:
    hits=0
    for a in np.linspace(0,0.5,9):
        for b in np.linspace(math.sqrt(1-a*a),2.0,10):
            B=np.array([[1.0,0.0],[a,b]])
            cell=voronoi_cell(B)
            for M in sublattices(k):
                g=feasible_x(B,M,cell,XMAX=6.0)
                if g:
                    hits+=1
                    if hits<=3: print(f"k={k} a={a:.3f} b={b:.3f} M={M.tolist()} x in {[(round(p,4),round(q,4)) for p,q in g]}",flush=True)
    print(f"k={k}: {hits} hits ({time.time()-t0:.0f}s)",flush=True)
