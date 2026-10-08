"""Independent re-check of the scalar inequalities of Proposition 3.1 (mpmath interval arithmetic).

The verifier (cap/verify_step.cpp, CAPD) writes, for each case, the enclosures it computed:
the growth factor G, the return time, the end time s_end of the hull, the bound q_KA on the
last window shell, and the enclosure of the renormalized image of the box.  This script reads
those numbers, widens each of them outward by a relative 1e-14 (so that the decimal printing
of the doubles cannot help), and re-checks in 60-digit interval arithmetic, with formulas typed
in afresh here:

  (deep)   1 - W_d R_{-KB-1} s_end > 0 and W_d / (1 - W_d R_{-KB-1} s_end) / G_lo <= W_d;
  (ahead)  R_{KA+1} s_end q_KA^2 <= theta,  r/2 < 1,  c1 = 4 R_{KA+2} s_end theta <= 1/2,
           1/2 + c1 <= G_lo;
  (eps)    Lambda^{2 alpha} / (r G_lo) < 1 for the largest alpha of the case;
  (box)    the image enclosure lies in the open box  |Y_k - w*_k| < rho_k, the centre w* and
           radii rho read from cap/cases/ (not from the verifier's output);
  (start)  h = Y_2 - c Y_1 < 0 on the closed box at the start, i.e. w*_2 + rho_2 < c^2 = 1/4;
  (growth) G_lo > 1, and the exponent gamma = 1/3 - log G_lo / log Lambda.

It also prints the constants used in Theorem 1.1 (Table 1 of the paper).  It does not repeat
the ODE enclosures, which only CAPD computes.
"""
import re
import sys
from pathlib import Path

from mpmath import iv, mp, mpf

iv.dps = 60
mp.dps = 30
ROOT = Path(__file__).resolve().parent.parent
CASES = [  # name, KB, KA, Lambda, alpha_hi (as fractions)
    ("L6_5_K80_7", 80, 7, (6, 5), (377, 1000)),
    ("L13_10_K60_6", 60, 6, (13, 10), (189, 500)),
    ("L3_2_K50_8", 50, 8, (3, 2), (73, 200)),
    ("L17_10_K40_8", 40, 8, (17, 10), (87, 250)),
    ("L9_5_K40_8", 40, 8, (9, 5), (339, 1000)),
]
WD = iv.mpf(1)
THETA = iv.mpf("1e-31")
EPS0 = iv.mpf("1e-12")
WIDEN = mpf("1e-14")


def lo(x):
    """A lower bound for a logged lower bound."""
    x = mpf(x)
    return iv.mpf(x - abs(x) * WIDEN - mpf("1e-300"))


def hi(x):
    x = mpf(x)
    return iv.mpf(x + abs(x) * WIDEN + mpf("1e-300"))


def ok(cond, what, failures):
    if not cond:
        failures.append(what)
    return cond


def check(name, KB, KA, lam, ahi):
    failures = []
    log = (ROOT / "cap/logs" / f"{name}.log").read_text()
    res = (ROOT / "cap/logs" / f"{name}_result.txt").read_text().splitlines()
    if "PASS" not in log.split("\n"):
        failures.append("verifier did not print PASS")
    Lam = iv.mpf(lam[0]) / lam[1]
    alpha = iv.mpf(ahi[0]) / ahi[1]
    r = Lam ** (iv.mpf(2) / 3)
    Rk = lambda k: r ** k
    g = re.search(r"G true in \[([^,]+), ([^\]]+)\]", log)
    Glo, Ghi = lo(g.group(1)), hi(g.group(2))
    t = [l for l in res if l.startswith("returnTime")][0].split()
    send = hi(t[4])
    q = re.search(r"ahead tail: qKA=\[[^,]+, ([^\]]+)\]", log)
    qKA = hi(q.group(1))
    # deep tail
    den = 1 - WD * Rk(-KB - 1) * send
    ok(den.a > 0, "deep Riccati denominator", failures)
    ok((WD / den / Glo).b <= WD.a, "deep tail closure", failures)
    # ahead tail
    ok((Rk(KA + 1) * send * qKA ** 2).b <= THETA.a, "ahead base", failures)
    ok((r / 2).b < 1, "r/2 < 1", failures)
    c1 = 4 * Rk(KA + 2) * send * THETA
    ok(c1.b <= 0.5, "ahead cond1", failures)
    ok((iv.mpf(0.5) + c1).b <= Glo.a, "ahead cond2", failures)
    # viscosity contraction
    epsfac = Lam ** (2 * alpha) / (r * Glo)
    ok(epsfac.b < 1, "eps contraction", failures)
    # box inclusion against the stored centre and radii
    cen = {}
    for line in (ROOT / "cap/cases" / f"center_{name}.txt").read_text().split("\n"):
        if line.strip():
            k, v = line.split()
            cen[int(k)] = mpf(v)
    cen[0], cen[1] = mpf(1), mpf("0.5")
    rad = {}
    for line in (ROOT / "cap/cases" / f"radii_{name}.txt").read_text().split("\n"):
        if line.strip():
            k, v = line.split()
            rad[int(k)] = mpf(v)
    # the section function is negative on the closed start box (Y_1 = c = 1/2)
    ok((iv.mpf(cen[2]) + iv.mpf(rad[2]) - iv.mpf("0.25")).b < 0, "h < 0 at the start", failures)
    worst = mpf(0)
    nk = 0
    for line in res:
        if line.startswith("#") or not line[:1] in "-0123456789":
            continue
        f = line.split()
        k = int(f[0])
        plo, phi = lo(f[3]), hi(f[4])
        c = iv.mpf(cen[k])
        rr = iv.mpf(rad[k])
        dev = max(abs(mpf((phi - c).b)), abs(mpf((plo - c).a)))
        worst = max(worst, dev / mpf(rr.a))
        ok((phi - c).b < rr.a and (c - plo).b < rr.a, f"box inclusion k={k}", failures)
        nk += 1
    ok(nk == KB + KA - 1, "number of window coordinates", failures)
    ok(Glo.a > 1, "growth G > 1", failures)
    gamma = iv.mpf(1) / 3 - iv.log(Glo) / iv.log(Lam)
    athr = iv.mpf(1) / 3 + iv.log(Glo) / (2 * iv.log(Lam))
    # constants of Theorem 1.1: T_max <= C_T / A, A >= C_A * nu with C_A at alpha = alpha_hi
    CT = Lam ** (iv.mpf(1) / 3) * send / (r ** KB * (1 - 1 / (r * Glo)))
    CA = Lam ** (iv.mpf(1) / 3) * Lam ** ((2 * alpha - iv.mpf(2) / 3) * KB) / EPS0
    return failures, dict(
        Glo=mpf(Glo.a), Ghi=mpf(Ghi.b), send=mpf(send.b), gamma=mpf(gamma.b), athr=mpf(athr.a),
        epsfac=mpf(epsfac.b), worst=worst, CT=mpf(CT.b), CA=mpf(CA.b),
    )


def main():
    allok = True
    print(f"{'case':14s} {'G_lo':>20s} {'G_hi':>20s} {'s_end':>12s} {'gamma<=':>9s} "
          f"{'a_thr>=':>9s} {'epsfac<=':>10s} {'worst':>7s} {'C_T':>10s} {'C_A':>10s}  result")
    for name, KB, KA, lam, ahi in CASES:
        fails, d = check(name, KB, KA, lam, ahi)
        allok &= not fails
        print(f"{name:14s} {mp.nstr(d['Glo'], 17):>20s} {mp.nstr(d['Ghi'], 17):>20s} "
              f"{mp.nstr(d['send'], 8):>12s} {mp.nstr(d['gamma'], 5):>9s} {mp.nstr(d['athr'], 5):>9s} "
              f"{mp.nstr(d['epsfac'], 6):>10s} {mp.nstr(d['worst'], 4):>7s} {mp.nstr(d['CT'], 4):>10s} "
              f"{mp.nstr(d['CA'], 4):>10s}  {'ok' if not fails else 'FAILED: ' + ', '.join(fails)}")
    print("all scalar inequalities re-checked" if allok else "RE-CHECK FAILED")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
