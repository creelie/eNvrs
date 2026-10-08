// Computer-assisted verification of one renormalized front step for the Katz-Pavlovic dyadic model
//   du_n/dt = -nu N_n^{2 alpha} u_n + N_{n-1} u_{n-1}^2 - N_n u_n u_{n+1},  N_n = Lambda^n.
//
// Front coordinates at front index n (window k = -KB..KA):
//   Y_k = y_{n+k}/y_n(start),  y_m = Lambda^{m/3} u_m,  ds = r^n y_n(start) dtau,  r = Lambda^{2/3},
//   dY_k/ds = R_k (Y_{k-1}^2 - Y_k Y_{k+1}) - E_k Y_k,   R_k = r^k,  E_k = eps Lambda^{2 alpha k},
// eps in [0, eps0] and alpha in [alphaLo, alphaHi] are constant parameters during a step.
// Start: Y_0 = 1, Y_1 = c.  End: first crossing of h = Y_2 - c Y_1 from - to +.
// Step map: Y'_k = Y_{k+1}(s*)/G, G = Y_1(s*), eps' = eps Lambda^{2 alpha}/(r G).
//
// The truncated field sets Y_{-KB-1} = Y_{KA+1} = 0.  The true window equations differ by
//   e_{-KB} in [0, R_{-KB} Wd^2],  e_{KA} in [-R_{KA} q_KA q_{KA+1}, 0],
// whose effect is bounded by a cooperative comparison argument (the paper, Lemma 4.3).
// All arithmetic is interval arithmetic (CAPD, filib); the program prints PASS only if every
// inequality (V1)-(V9) of Section 5 of the paper is verified.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "capd/capdlib.h"
#define MYHOE_EPS 1.e-30
#include "MyHOE.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/TimeMap.hpp"
#include "capd/poincare/PoincareMap.hpp"
using namespace std;
using namespace capd;
using capd::autodiff::Node;
// CAPD solver whose high-order enclosure uses a trial-remainder floor of 1e-30 instead of 1e-300.
// The floor only enlarges the trial set; the inclusion test that validates the enclosure is unchanged.
typedef capd::dynsys::OdeSolver<IMap, capd::dynsys::ILastTermsStepControl, capd::dynsys::MyHOE> MSolver;
typedef capd::poincare::PoincareMap<MSolver> MPoincareMap;
typedef capd::poincare::TimeMap<MSolver> MTimeMap;

static int KB = 60, KA = 6;
inline int nn() { return KB + KA + 1; }
inline int idx(int k) { return k + KB; }

void vfield(Node, Node in[], int, Node out[], int, Node params[], int) {
  for (int k = -KB; k <= KA; ++k) {
    Node Y = in[idx(k)];
    Node R = params[idx(k)], E = params[nn() + idx(k)];
    Node tr;
    if (k == -KB) tr = -(R * (Y * in[idx(k + 1)]));
    else if (k == KA) tr = R * (in[idx(k - 1)] ^ 2);
    else tr = R * ((in[idx(k - 1)] ^ 2) - Y * in[idx(k + 1)]);
    out[idx(k)] = tr - E * Y;
  }
}

static interval LAM, RR, ALPHA, C, EPS;
interval Rk(int k) { return exp(interval(2 * k) / interval(3) * log(LAM)); }

void setParams(IMap& f) {
  for (int k = -KB; k <= KA; ++k) {
    f.setParameter(idx(k), Rk(k));
    f.setParameter(nn() + idx(k), EPS * exp(interval(2.0) * ALPHA * interval(k) * log(LAM)));
  }
}

// Phi: P (state at end section) -> next state at start section (truncated tail: Y'_KA = 0)
IVector Phi(const IVector& P) {
  IVector y(nn());
  interval G = P[idx(1)];
  for (int k = -KB; k <= KA; ++k) {
    if (k == 0) { y[idx(k)] = interval(1.0); continue; }
    if (k == 1) { y[idx(k)] = C; continue; }
    if (k + 1 > KA) { y[idx(k)] = interval(0.0); continue; }
    y[idx(k)] = P[idx(k + 1)] / G;
  }
  return y;
}
IMatrix DPhi(const IVector& P) {
  IMatrix J(nn(), nn());
  interval G = P[idx(1)];
  for (int k = -KB; k <= KA; ++k) {
    if (k == 0 || k == 1 || k + 1 > KA) continue;
    J[idx(k)][idx(k + 1)] += interval(1.0) / G;
    J[idx(k)][idx(1)] += -P[idx(k + 1)] / sqr(G);
  }
  return J;
}

int main(int argc, char** argv) {
  cout.precision(17);
  if (argc < 15) {
    cerr << "usage: step2 KB KA LamNum LamDen aLoNum aLoDen aHiNum aHiDen eps0 center radii out Wd theta\n";
    return 1;
  }
  KB = atoi(argv[1]); KA = atoi(argv[2]);
  LAM = interval(atof(argv[3])) / interval(atof(argv[4]));
  argv += 2;
  interval aLo = interval(atof(argv[3])) / interval(atof(argv[4]));
  interval aHi = interval(atof(argv[5])) / interval(atof(argv[6]));
  ALPHA = interval(aLo.leftBound(), aHi.rightBound());
  double eps0 = atof(argv[7]);
  EPS = interval(0.0, eps0);
  string centerFile = argv[8], radiiFile = argv[9], out = argv[10];
  interval Wd = interval(atof(argv[11]));     // deep tail bound (normalized), hypothesis (H_d)
  interval theta = interval(atof(argv[12]));  // ahead tail: Y_k <= theta * (1/2)^{k-KA-1}, k > KA
  RR = Rk(1);
  C = interval(0.5);
  int d = nn();
  ofstream logf(out + ".log");
  logf.precision(17);
  logf << "Lambda=" << LAM << " KB=" << KB << " KA=" << KA << " alpha in " << ALPHA << " eps in " << EPS << " c=" << C
      << " Wd=" << Wd << " theta=" << theta << "\n";

  IMap vf(vfield, d, d, 2 * d);
  setParams(vf);
  MSolver solver(vf, 20);
  IVector origin(d), normal(d);
  normal[idx(2)] = interval(1.0);
  normal[idx(1)] = -C;
  IAffineSection section(origin, normal);
  MPoincareMap pm(solver, section, poincare::MinusPlus);

  IVector x0(d), r0(d);
  {
    ifstream in(centerFile);
    int k; double v;
    while (in >> k >> v) if (k >= -KB && k <= KA) x0[idx(k)] = interval(v);
  }
  x0[idx(0)] = interval(1.0);
  x0[idx(1)] = C;
  bool pointMode = (radiiFile == "none");
  if (!pointMode) {
    ifstream in(radiiFile);
    int k; double v;
    while (in >> k >> v) {
      if (k == 0 || k == 1 || k < -KB || k > KA) continue;
      r0[idx(k)] = interval(-v, v);
    }
  }

  // 1) center
  C1Rect2Set sc(x0);
  IMatrix monoC(d, d);
  interval rtC;
  IVector Pc = pm(sc, monoC, rtC);
  IMatrix DPc = pm.computeDP(Pc, monoC, rtC);
  IVector PhiC = Phi(Pc);
  logf << "center: return time " << rtC << "  G " << Pc[idx(1)] << "\n";
  {
    IMatrix JC = DPhi(Pc) * DPc;
    ofstream jf(out + "_Jcenter.txt");
    jf.precision(17);
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j) jf << mid(JC[i][j]).leftBound() << (j + 1 < d ? " " : "\n");
    ofstream pf(out + "_PhiCenter.txt");
    pf.precision(17);
    for (int k = -KB; k <= KA; ++k)
      pf << k << " " << x0[idx(k)].leftBound() << " " << PhiC[idx(k)].leftBound() << " " << PhiC[idx(k)].rightBound() << "\n";
  }
  if (pointMode) { logf << "point mode done\n"; return 0; }

  bool ok = true;
  auto fail = [&](const string& s) { ok = false; logf << "FAIL: " << s << "\n"; };

  // 2) C1 Poincare map of the box
  C1Rect2Set sB(x0, r0);
  IMatrix monoB(d, d);
  interval rtB;
  IVector PB = pm(sB, monoB, rtB);
  IMatrix DPB = pm.computeDP(PB, monoB, rtB);
  IVector PBmv = Pc + DPB * r0;
  IVector PBi(d);
  for (int i = 0; i < d; ++i)
    if (!intersection(PB[i], PBmv[i], PBi[i])) fail("empty P intersection");
  IMatrix DPhiB = DPhi(PBi) * DPB;
  IVector PhiMV = PhiC + DPhiB * r0;
  IVector PhiD = Phi(PBi);
  IVector PhiB(d);
  for (int i = 0; i < d; ++i)
    if (!intersection(PhiMV[i], PhiD[i], PhiB[i])) fail("empty Phi intersection");
  logf << "box: return time " << rtB << "  G " << PBi[idx(1)] << "\n";

  // 3) hull of truncated trajectories over [0, sEnd]
  // optional settings (environment): SEXTRA = length of the hull beyond the return time (default 0.02),
  // HMAX = maximal time step of the hull integration (default: none)
  double sExtra = getenv("SEXTRA") ? atof(getenv("SEXTRA")) : 0.02;
  double sEnd = rtB.rightBound() + sExtra;
  MSolver solver2(vf, 20);
  if (getenv("HMAX")) solver2.setMaxStep(atof(getenv("HMAX")));
  MTimeMap tm(solver2);
  C0Rect2Set hs(x0, r0);
  tm.stopAfterStep(true);
  vector<IVector> encl;
  vector<interval> times;
  IVector hull = x0 + r0;
  interval prev(0.0);
  do {
    tm(sEnd, hs);
    interval st = solver2.getStep();
    const MSolver::SolutionCurve& curve = solver2.getCurve();
    IVector v = curve(interval(0, 1) * st);
    encl.push_back(v);
    times.push_back(prev + interval(0, 1) * st);
    hull = intervalHull(hull, v);
    prev = tm.getCurrentTime();
  } while (!tm.completed());
  logf << "hull: " << encl.size() << " steps up to s=" << prev << " (sExtra " << sExtra << ")\n";

  // 4) tails
  interval Rdeep = Rk(-KB - 1);
  interval sMax = interval(sEnd);
  // deep: max_{k<-KB} Y_k(s) <= Wd/(1 - Wd Rdeep s)
  interval den = interval(1.0) - Wd * Rdeep * sMax;
  if (!(den.leftBound() > 0)) fail("deep tail Riccati");
  interval WdHat = Wd / den;
  // a priori perturbation bound D (bootstrap), iterate twice
  IVector D(d);
  for (int i = 0; i < d; ++i) D[i] = interval(0.0);
  IVector z(d);
  interval qKA, qKA1;
  bool bootOk = false;
  for (int it = 0; it < 6; ++it) {
    IVector hullD = hull;
    for (int i = 0; i < d; ++i) hullD[i] = hull[i] + interval(-1, 1) * D[i];
    qKA = interval(0.0, std::max(abs(hullD[idx(KA)].leftBound()), abs(hullD[idx(KA)].rightBound())));
    // ahead tail: q_{KA+1} = theta_{KA+1} + R_{KA+1} sEnd q_KA^2
    qKA1 = theta + Rk(KA + 1) * sMax * sqr(qKA);
    IVector ebar(d);
    ebar[idx(-KB)] = Rk(-KB) * sqr(WdHat);
    ebar[idx(KA)] = Rk(KA) * qKA * qKA1;
    // Metzler majorant of the Jacobian over hullD
    IMatrix J = vf[hullD];
    vector<vector<double>> M(d, vector<double>(d, 0.0));
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j) {
        if (i == j) M[i][j] = J[i][j].rightBound();
        else M[i][j] = std::max(abs(J[i][j].leftBound()), abs(J[i][j].rightBound()));
      }
    // zeta = (lambda0 I - M)^{-1} (ebar + floor), solved by Gauss-Seidel style iteration (nonrigorous)
    double lambda0 = 1.0;
    for (int i = 0; i < d; ++i) {
      double rs = M[i][i];
      for (int j = 0; j < d; ++j) if (j != i) rs += M[i][j];
      lambda0 = std::max(lambda0, rs + 1.0);
    }
    vector<double> zeta(d, 0.0), rhs(d);
    for (int i = 0; i < d; ++i) rhs[i] = ebar[i].rightBound() + 1e-250;
    for (int sweep = 0; sweep < 4000; ++sweep) {
      double ch = 0;
      for (int i = 0; i < d; ++i) {
        double s = rhs[i];
        for (int j = 0; j < d; ++j) if (j != i) s += M[i][j] * zeta[j];
        double nz = s / (lambda0 - M[i][i]);
        ch = std::max(ch, abs(nz - zeta[i]) / std::max(nz, 1e-300));
        zeta[i] = nz;
      }
      if (ch < 1e-14) break;
    }
    // rigorous lambda and eta
    interval lam(0.0), eta(0.0);
    for (int i = 0; i < d; ++i) {
      interval s(0.0);
      for (int j = 0; j < d; ++j) s += interval(M[i][j]) * interval(zeta[j]);
      lam = interval(std::max(lam.rightBound(), (s / interval(zeta[i])).rightBound()));
      eta = interval(std::max(eta.rightBound(), (ebar[i] / interval(zeta[i])).rightBound()));
    }
    interval fac = (exp(lam * sMax) - interval(1.0)) / lam;
    bool fits = true;
    for (int i = 0; i < d; ++i) {
      z[i] = eta * interval(zeta[i]) * fac;
      if (!(z[i].rightBound() < D[i].rightBound())) fits = false;  // strict, for the continuity argument
    }
    logf << "bootstrap it " << it << ": lambda=" << lam << " eta=" << eta << " z_front(-1,1,2)=" << z[idx(-1)].rightBound()
        << "," << z[idx(1)].rightBound() << "," << z[idx(2)].rightBound() << " z_deep=" << z[idx(-KB)].rightBound()
        << " z_KA=" << z[idx(KA)].rightBound() << "\n";
    if (fits) { bootOk = true; break; }
    for (int i = 0; i < d; ++i) D[i] = interval(2.0) * z[i] + interval(1e-300);
  }
  if (!bootOk) fail("perturbation bootstrap");

  // 5) section crossing for the true trajectories
  {
    // at s = 0 the true and truncated states coincide and lie in the closed box, where h < 0
    interval h0 = (x0[idx(2)] + r0[idx(2)]) - C * (x0[idx(1)] + r0[idx(1)]);
    if (!(h0.rightBound() < 0)) fail("section: h not negative at the start");
    logf << "h at start <= " << h0.rightBound() << "\n";
  }
  interval m(1e300), xi = z[idx(2)] + C * z[idx(1)];
  IVector F(d);
  for (int i = 0; i < d; ++i) F[i] = interval(0.0);
  int nCross = 0;
  for (size_t j = 0; j < encl.size(); ++j) {
    IVector v = encl[j];
    for (int i = 0; i < d; ++i) v[i] += interval(-1, 1) * z[i];
    interval h = v[idx(2)] - C * v[idx(1)];
    IVector fv = vf(v);
    interval dh = fv[idx(2)] - C * fv[idx(1)];
    if (h.rightBound() < 0) continue;
    if (dh.leftBound() > 0) {
      ++nCross;
      m = interval(std::min(m.leftBound(), dh.leftBound()));
      for (int i = 0; i < d; ++i) F[i] = interval(0.0, std::max(F[i].rightBound(), std::max(abs(fv[i].leftBound()), abs(fv[i].rightBound()))));
    } else {
      fail("section: neither negative nor monotone on a step");
    }
  }
  {
    // at s = sEnd both the truncated and the true trajectories are past the section: h > 0
    IVector v = IVector(hs);
    interval h = v[idx(2)] - C * v[idx(1)];
    if (!(h.leftBound() - xi.rightBound() > 0)) fail("section: not crossed by sEnd");
  }
  F[idx(-KB)] += interval(0.0, (Rk(-KB) * sqr(WdHat)).rightBound());
  F[idx(KA)] += interval(0.0, (Rk(KA) * qKA * qKA1).rightBound());
  logf << "crossing steps " << nCross << " m=" << m << " xi=" << xi << "\n";
  IVector Infl(d);
  for (int i = 0; i < d; ++i) Infl[i] = interval(-1, 1) * (z[i].rightBound() + (F[i] * xi / m).rightBound());

  // 6) true image
  IVector PTrue = PBi + Infl;
  // ahead tail chain: for k >= KA+2, theta_{k-1} >= q_k/G with theta_k = theta 2^{-(k-KA-1)}
  interval Glo = interval(PTrue[idx(1)].leftBound());
  {
    interval beta = interval(0.5);
    interval th1 = theta;  // theta_{KA+1}
    // base: R_{KA+1} sEnd qKA^2 <= theta
    if (!((Rk(KA + 1) * sMax * sqr(qKA)).rightBound() <= theta.leftBound())) fail("ahead base");
    // beta r < 1 and the k = KA+2 conditions (monotone in k afterwards)
    if (!((beta * RR).rightBound() < 1.0)) fail("beta r < 1");
    interval cond1 = interval(4.0) * Rk(KA + 2) * sMax * th1;  // <= beta
    if (!(cond1.rightBound() <= beta.leftBound())) fail("ahead cond1");
    if (!((beta + cond1).rightBound() <= Glo.leftBound())) fail("ahead cond2");
    logf << "ahead tail: qKA=" << qKA << " qKA1=" << qKA1 << "\n";
  }

  IMatrix DPhiT = DPhi(PTrue);
  IVector PhiTrue = PhiB + DPhiT * Infl;
  // ahead tail contribution to the last window coordinate: Y'_KA in [0, qKA1/G]
  PhiTrue[idx(KA)] = PhiTrue[idx(KA)] + interval(0.0, (qKA1 / Glo).rightBound());
  interval GTrue = PTrue[idx(1)];
  logf << "G true in " << GTrue << "\n";

  // 7) checks
  double worst = 0;
  int worstK = 0;
  ofstream rf(out + "_result.txt");
  rf.precision(17);
  rf << "# k center radius Phi_lo Phi_hi ratio\n";
  for (int k = -KB; k <= KA; ++k) {
    if (k == 0 || k == 1) continue;
    int i = idx(k);
    double rad = r0[i].rightBound();
    double dev = std::max(abs((PhiTrue[i] - x0[i]).leftBound()), abs((PhiTrue[i] - x0[i]).rightBound()));
    double ratio = dev / rad;
    rf << k << " " << x0[i].leftBound() << " " << rad << " " << PhiTrue[i].leftBound() << " " << PhiTrue[i].rightBound() << " " << ratio << "\n";
    if (!(ratio < 1.0)) fail("box inclusion at k=" + to_string(k));
    if (ratio > worst) { worst = ratio; worstK = k; }
  }
  logf << "worst inclusion ratio " << worst << " at k=" << worstK << "\n";
  // deep tail closure: Y_{-KB}(s*)/G <= Wd and WdHat/G <= Wd
  if (!((PTrue[idx(-KB)] / GTrue).rightBound() <= Wd.leftBound())) fail("deep exit value");
  if (!((WdHat / GTrue).rightBound() <= Wd.leftBound())) fail("deep tail closure");
  // eps contraction and growth versus viscosity
  interval epsFac = exp(interval(2.0) * interval(aHi.rightBound()) * log(LAM)) / (RR * GTrue);
  logf << "eps factor Lambda^{2 alphaHi}/(r G) <= " << epsFac.rightBound() << "\n";
  if (!(epsFac.rightBound() < 1.0)) fail("eps contraction / growth");
  // positivity of the start section values used for normalisation
  if (!(GTrue.leftBound() > 0)) fail("G positive");
  rf << "G " << GTrue.leftBound() << " " << GTrue.rightBound() << "\n";
  rf << "returnTime " << rtB.leftBound() << " " << rtB.rightBound() << " sEnd " << sEnd << "\n";
  rf << "epsFactor " << epsFac.rightBound() << "\n";
  logf << (ok ? "PASS" : "NOT VERIFIED") << "\n";
  cout << (ok ? "PASS" : "NOT VERIFIED") << "  worst ratio " << worst << " at k=" << worstK << "  G in " << GTrue << "\n";
  return ok ? 0 : 2;
}
