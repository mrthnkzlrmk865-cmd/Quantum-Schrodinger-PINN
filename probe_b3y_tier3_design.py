# DESIGN-VERIFICATION PROBE for the Tier 3 concept. Synthetic single-channel chain (Ny = 1); NOT a B3Y comparison:
# no B3Y point, no RGF, no B3X reference is used. Question: does "extended-precision residual of the UNMERGED operator"
# separate (i) operator representation error, (ii) linear-solver error, and does the refinement converge?
import os; os.environ["OPENBLAS_NUM_THREADS"]="1"
import numpy as np, mpmath as mp, scipy.sparse as sp, scipy.sparse.linalg as spla, math
mp.mp.dps = 60
LD, CLD = np.longdouble, np.clongdouble
print("longdouble: eps=%.3e, mantissa bits=%d"%(np.finfo(LD).eps, np.finfo(LD).nmant+1))
def bump(x, a, b, d, A): return A*0.5*(np.tanh((x-a)/d)-np.tanh((x-b)/d))
Eperp = 0.5*(4*math.pi)**2          # channel-4 threshold (worst-case-like near-threshold channel)
E = 88.5134459632352                # same E as the worst B3X point; here only a label of the energy scale
for inv in (512, 1024, 2048):
    dx = 1.0/inv; nb = 103*inv//512; x0 = 1.9-nb*dx; N = 2*nb+inv+1
    xs = x0+dx*np.arange(N); V = bump(xs,2.35,2.45,0.01,5.0)+bump(xs,2.60,2.75,0.02,3.0)   # doubles = the exact inputs
    t = 1.0/(2*dx*dx)                                                                       # exact power of two
    # ---------- exact lattice answer: mpmath transfer recursion from the exact inputs ----------
    x_ex = (mp.mpf(E)-mp.mpf(Eperp))*mp.mpf(dx)**2
    th = 2*mp.asin(mp.sqrt(x_ex/2)); lam_ex = mp.expjpi(th/mp.pi)
    psi_next, psi = lam_ex, mp.mpf(1)
    for j in range(N, 0, -1):
        prev = ((2*mp.mpf(t)+mp.mpf(Eperp)+mp.mpf(float(V[j-1]))-mp.mpf(E))*psi - mp.mpf(t)*psi_next)/mp.mpf(t)
        psi_next, psi = psi, prev
    psi1, psi0 = psi_next, psi
    a_ex = (psi1*lam_ex - psi0)/(lam_ex-1/lam_ex)*(1) if False else None
    A = mp.matrix([[1,1],[1/lam_ex, lam_ex]]); sol = mp.lu_solve(A, mp.matrix([psi1, psi0])); a_, b_ = sol[0], sol[1]
    t_exact, r_exact = 1/a_, b_/a_
    # ---------- double-precision method (like the sparse references): rounded matrix + LU ----------
    c = 1.0-(E-Eperp)*dx*dx; s = math.sqrt(1.0-c*c); lam_d = complex(c, s); gam_d = 2*t*s
    diag = ((E-2*t)-Eperp-V).astype(complex); diag[0] += t*lam_d; diag[-1] += t*lam_d
    M = sp.diags([np.full(N-1,t,complex), diag, np.full(N-1,t,complex)], [-1,0,1], format="csc")
    lu = spla.splu(M); src = -t*(1/lam_d-lam_d)/math.sqrt(gam_d)
    b = np.zeros(N, complex); b[0] = src
    xd = lu.solve(b)
    def S_from(x1, xN, lam, gam):                      # wave-function extraction, left incidence, one channel
        a = gam**-0.5; return np.sqrt(gam)*(x1-a), np.sqrt(gam)*xN        # (r, t)
    r_d, t_d = S_from(xd[0], xd[-1], lam_d, gam_d)
    # ---------- extended precision pieces ----------
    Eext, Epext, tL = LD(E), LD(Eperp), LD(t); VL = V.astype(LD)
    xL = (Eext-Epext)*LD(dx)**2; thL = 2*np.arcsin(np.sqrt(xL/2)); lamL = np.exp(1j*thL).astype(CLD); gamL = 2*tL*np.sin(thL)
    P = (Eext-Epext-VL).astype(CLD)                                          # the UNMERGED physical part
    def Aop_exact(c_):                                                      # P c + t (c_{j-1} - 2 c_j + c_{j+1}), ghosts c_0 = lam c_1, c_{N+1} = lam c_N
        cm = np.empty_like(c_); cp = np.empty_like(c_)
        cm[1:] = c_[:-1]; cm[0] = lamL*c_[0]; cp[:-1] = c_[1:]; cp[-1] = lamL*c_[-1]
        return P*c_ + tL*(cm-2*c_+cp)
    srcL = -tL*(1/lamL-lamL)*gamL**LD(-0.5)
    bL = np.zeros(N, CLD); bL[0] = srcL
    def refine(apply_op, x_start, label, iters=8):
        x = x_start.astype(CLD); hist = []
        for k in range(iters):
            r = bL-apply_op(x); rn = float(np.max(np.abs(r))); dxk = lu.solve(r.astype(complex)); x = x+dxk.astype(CLD)
            hist.append((rn, float(np.max(np.abs(dxk)))))
        return x, hist
    # exact-operator refinement (representation error removed)
    x_ex_op, h_ex = refine(Aop_exact, xd, "exact operator")
    # rounded-matrix refinement: residual against the ROUNDED double matrix entries (exactly), in extended precision
    Md = M.astype(complex)
    def Aop_rounded(c_):
        c64 = c_; out = np.empty_like(c64)
        dg = diag.astype(CLD); out = dg*c64
        out[1:] += LD(t)*c64[:-1]; out[:-1] += LD(t)*c64[1:]
        return out
    x_rd_op, h_rd = refine(Aop_rounded, xd, "rounded matrix")
    r_x, t_x = S_from(x_ex_op[0], x_ex_op[-1], lamL, gamL)      # exact-lattice result (extended precision data)
    r_m, t_m = S_from(x_rd_op[0], x_rd_op[-1], lam_d, gam_d)    # solution of the rounded matrix with double lead data
    def to_mp(z): 
        f=lambda v: mp.mpf(np.format_float_positional(v, precision=24, unique=False))
        return mp.mpc(f(z.real), f(z.imag))
    # operator B: ROUNDED bulk diagonal, but EXACT (extended) lead terms and exact extraction -> isolates the bulk rounding
    dgB = diag.astype(CLD).copy(); tl_d = CLD(t*lam_d)
    dgB[0] += -tl_d + tL*lamL; dgB[-1] += -tl_d + tL*lamL
    def Aop_B(c_):
        out = dgB*c_; out[1:] += tL*c_[:-1]; out[:-1] += tL*c_[1:]; return out
    x_B, h_B = refine(Aop_B, xd, "bulk-rounded, exact lead")
    r_B, t_B = S_from(x_B[0], x_B[-1], lamL, gamL)
    e_tot = abs(complex(t_d)-complex(t_exact)); e_xp = float(abs(to_mp(np.complex128(0)+t_x)-t_exact)) if False else float(abs(to_mp(t_x)-t_exact))
    e_bulk = float(abs(to_mp(t_B)-t_exact))
    e_rep = abs(complex(t_m)-complex(t_exact)); e_sol = abs(complex(t_d)-complex(t_m))
    print("\ndx=1/%d N=%d  |t_exact|^2=%.6f  sin(k dx)=%.2e"%(inv,N,float(abs(t_exact)**2),float(mp.sin(th))))
    print("  error of the double-precision solve  |t_dbl - t_exact|          = %.2e"%e_tot)
    print("  extended-refined exact-operator      |t_XP  - t_exact|          = %.2e   (measured in mpmath; the Tier 3 candidate)"%e_xp)
    print("  bulk-rounded diagonal + exact lead   |t_B   - t_exact|          = %.2e   <- bulk operator rounding alone"%e_bulk)
    print("  rounded-matrix solution (extended)   |t_fl(A) - t_exact|        = %.2e   <- operator representation error (+lead-data error)"%e_rep)
    print("  double LU vs rounded matrix solution |t_dbl - t_fl(A)|          = %.2e   <- linear-solver error"%e_sol)
    print("  refinement history (exact op): residual, update:",["%.1e/%.1e"%h for h in h_ex[:5]])
    print("  refinement history (rounded) : residual, update:",["%.1e/%.1e"%h for h in h_rd[:5]])
