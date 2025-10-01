"""
topopt_rect.py
SIMP topology optimization 2D (plane stress), fixed base, load at middle top.
Run: python topopt_rect.py
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# --------------------
# Problem parameters
nelx = 80        # numero elementi in x (modifica per risoluzione)
nely = 40        # numero elementi in y
volfrac = 0.4    # frazione di volume permessa
penal = 3.0      # penalizzazione SIMP
rmin = 1.5       # filtro di sensibilità (radius)
E0 = 1.0         # modulo elastico materiale
Emin = 1e-9      # piccolo valore per materiale "vuoto"
nu = 0.3         # Poisson ratio

# Solver / OC params
maxiter = 120
tol = 0.01

# --------------------
# Preparazione mesh / elementi
ndof = 2 * (nelx + 1) * (nely + 1)

# element stiffness matrix (4-node quad, plane stress)
def lk():
    E = 1.0
    nu = 0.3
    k = np.array([
        [  0.5 - nu/6,    0.125 + nu/8,  -0.25 - nu/12,  -0.125 + 3*nu/8,
           -0.25 + nu/12, -0.125 - nu/8,   nu/6,         0.125 - 3*nu/8],
        [ 0.125 + nu/8,   0.5 - nu/6,    0.125 - 3*nu/8,  nu/6,
          -0.125 - nu/8, -0.25 + nu/12, -0.125 + 3*nu/8, -0.25 - nu/12],
        [-0.25 - nu/12,  0.125 - 3*nu/8, 0.5 - nu/6,   -0.125 - nu/8,
           nu/6,         -0.125 + 3*nu/8, -0.25 + nu/12, 0.125 + nu/8],
        [-0.125 + 3*nu/8, nu/6,         -0.125 - nu/8,  0.5 - nu/6,
           0.125 + nu/8, -0.25 - nu/12,  0.125 - 3*nu/8, -0.25 + nu/12],
        [-0.25 + nu/12, -0.125 - nu/8,   nu/6,          0.125 + nu/8,
            0.5 - nu/6,  0.125 - 3*nu/8, -0.25 - nu/12, -0.125 + 3*nu/8],
        [-0.125 - nu/8, -0.25 + nu/12,  -0.125 + 3*nu/8, -0.25 - nu/12,
           0.125 - 3*nu/8, 0.5 - nu/6,   0.125 + nu/8,   nu/6],
        [ nu/6,         -0.125 + 3*nu/8, -0.25 + nu/12,  0.125 - 3*nu/8,
          -0.25 - nu/12, 0.125 + nu/8,    0.5 - nu/6,   -0.125 - nu/8],
        [ 0.125 - 3*nu/8, -0.25 - nu/12,  0.125 + nu/8,  -0.25 + nu/12,
          -0.125 + 3*nu/8, nu/6,         -0.125 - nu/8,  0.5 - nu/6]
    ])
    return k * E / (1 - nu**2)

KE = lk()

# prepare index vectors for assembling global stiffness
nodenrs = np.arange((nelx+1)*(nely+1)).reshape((nely+1, nelx+1))
edofVec = (2 * nodenrs[:-1,:-1] + 2).reshape(-1)
# build edof matrix (8 dofs per element)
edofMat = np.zeros((nelx*nely, 8), dtype=int)
count = 0
for ely in range(nely):
    for elx in range(nelx):
        n1 = (nodenrs[ely, elx])
        n2 = (nodenrs[ely, elx+1])
        n3 = (nodenrs[ely+1, elx+1])
        n4 = (nodenrs[ely+1, elx])
        edofs = np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n3, 2*n3+1, 2*n4, 2*n4+1])
        edofMat[count, :] = edofs
        count += 1

# force vector: apply downward point load at middle top node
F = np.zeros((ndof, 1))
mid = (nelx // 2)
top_node = (nely) * (nelx + 1) + mid
F[2*top_node + 1, 0] = -1.0  # downward

# fixed degrees: fix entire bottom edge (y=0)
fixed = []
for i in range(nelx+1):
    node = i
    fixed += [2*node, 2*node+1]
fixed = np.array(fixed)

free = np.setdiff1d(np.arange(ndof), fixed)

# filter: for each element compute neighbors within rmin
n_el = nelx * nely
iH = []
jH = []
sH = []
for ely in range(nely):
    for elx in range(nelx):
        e1 = ely*nelx + elx
        neigh_x_min = max(elx - int(np.floor(rmin)), 0)
        neigh_x_max = min(elx + int(np.floor(rmin)), nelx-1)
        neigh_y_min = max(ely - int(np.floor(rmin)), 0)
        neigh_y_max = min(ely + int(np.floor(rmin)), nely-1)
        for k in range(neigh_y_min, neigh_y_max+1):
            for l in range(neigh_x_min, neigh_x_max+1):
                e2 = k*nelx + l
                fac = rmin - np.sqrt((elx - l)**2 + (ely - k)**2)
                if fac > 0:
                    iH.append(e1)
                    jH.append(e2)
                    sH.append(fac)
H = sp.coo_matrix((sH, (iH, jH)), shape=(n_el, n_el)).tocsr()
Hs = np.array(H.sum(axis=1)).reshape(-1)

# initial design variables
x = volfrac * np.ones(n_el, dtype=float)
xold = x.copy()
dc = np.zeros(n_el, dtype=float)
ce = np.zeros(n_el, dtype=float)

# Start iterations
for loop in range(1, maxiter+1):
    # Assemble global stiffness
    sK = ((Emin + x**penal * (E0 - Emin))).repeat(8*8)  # placeholder
    K = sp.lil_matrix((ndof, ndof))
    # Efficient assembly: accumulate KE scaled per element
    for iel in range(n_el):
        edofs = edofMat[iel, :]
        ke = (Emin + x[iel]**penal * (E0 - Emin)) * KE
        for i in range(8):
            for j in range(8):
                K[edofs[i], edofs[j]] += ke[i, j]
    K = K.tocsr()

    # Solve system
    # split into free dofs
    # convert to csr then use spla.spsolve on submatrix
    Kff = K[free, :][:, free]
    Ff = F[free, 0]
    try:
        Uf = spla.spsolve(Kff, Ff)
    except Exception as e:
        print("Linear solver failed:", e)
        break
    U = np.zeros((ndof, ))
    U[free] = Uf

    # element-wise compliance and sensitivity
    for iel in range(n_el):
        edofs = edofMat[iel, :]
        ue = U[edofs]
        ce[iel] = (ue @ (KE @ ue))
        dc[iel] = -penal * (x[iel]**(penal-1)) * (E0 - Emin) * ce[iel]

    # filter sensitivities
    dcn = np.asarray((H @ (x * dc)) / Hs)  # filtered numerator
    xH = np.asarray((H @ x) / Hs)
    dc = dcn / (np.maximum(1e-9, xH))

    # Optimality Criteria update of design variables
    l1 = 0; l2 = 1e9; move = 0.2
    while (l2 - l1) / (l1 + l2 + 1e-9) > 1e-3:
        lmid = 0.5 * (l2 + l1)
        B = np.maximum(0.0, np.minimum(1.0, x * np.sqrt(-dc / lmid)))
        # apply move limits
        xnew = np.maximum(x - move, np.minimum(x + move, B))
        if xnew.mean() - volfrac > 0:
            l1 = lmid
        else:
            l2 = lmid
    x = xnew

    # compute change
    change = np.max(np.abs(x - xold))
    xold = x.copy()

    # print progress
    c_total = (ce * (Emin + x**penal * (E0 - Emin))).sum()
    print(f"Iter: {loop:3d}, Compliance: {c_total:.4f}, Vol: {x.mean():.3f}, ch: {change:.4f}")

    # stopping criteria
    if change < tol:
        print("Converged.")
        break

# plot result
plt.figure(figsize=(9,4))
x_plot = x.reshape((nely, nelx))
plt.imshow(1 - x_plot[::-1,:], cmap='gray', interpolation='none', extent=[0, nelx, 0, nely])
plt.title('Topologia ottimizzata (nero = materiale)')
plt.xlabel('nelx')
plt.ylabel('nely')
plt.colorbar(label='material density')
plt.tight_layout()
plt.savefig('topopt_result.png', dpi=300)
plt.show()

