"""Detail of the free-network outcome: pivot sizes after the chol fold, Cxx scale,
dof, and the test_fit fixture (testDatumIsADeclaredPrior) itself."""
import numpy as np, scipy.linalg as sla
from m1_qr_roundtrip import QRSystem
def detail(A, l):
    nc = A.shape[1]
    loc = QRSystem(nc); loc.add(A, l)
    R1raw = loc.R[:nc,:nc]
    _, _, N, _, n, ltl, m = loc.blocks(np.arange(nc), np.arange(nc, nc))
    g = QRSystem(nc); path = g.merge(N, n, ltl, m)
    R1 = g.R[:nc,:nc]; dR = np.abs(np.diag(R1)); tol = nc*np.spacing(dR.max())
    x, st = g.solve()
    out = dict(path=path, st=st, rawPivotRel=np.abs(np.diag(R1raw)).min()/np.abs(np.diag(R1raw)).max(),
               foldPivotRel=dR.min()/dR.max(), tolRel=tol/dR.max())
    if x is not None:
        Ri = np.linalg.inv(R1); Q = Ri@Ri.T
        out['maxSigmaX_cofactor'] = np.sqrt(np.diag(Q)).max()
        out['minNorm_x_dev'] = np.linalg.norm(x - np.linalg.pinv(A)@l)
    return out
# test_fit fixture: 4-point levelling network, the pinned refusal
pairs=[(0,1),(1,2),(2,3),(3,0),(0,2)]; A=np.zeros((5,4))
for k,(i,j) in enumerate(pairs): A[k,i]=-1; A[k,j]=1
l = A@np.array([0,1.5,-0.7,2.2]) + 0.01*np.random.RandomState(5).randn(5)
print('test_fit fixture:', detail(A,l))
# a few random rings with chords
rs=np.random.RandomState(0)
for u in (5,6,8,10,20,50):
    A=np.zeros((2*u,u))
    for i in range(u): A[i,i]=-1; A[i,(i+1)%u]=1
    for k in range(u):
        i,j=rs.choice(u,2,replace=False); A[u+k,i]=-1; A[u+k,j]=1
    l=A@rs.randn(u)+0.01*rs.randn(2*u)
    d=detail(A,l); print('u=%3d'%u, {k:(('%.2e'%v) if isinstance(v,float) else v) for k,v in d.items()})
