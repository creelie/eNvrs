/* Independent non-rigorous check, written separately from the CAPD verifier and the Python code.
 * Integrates the Katz-Pavlovic model in front coordinates on the full lattice (every deep shell,
 * AH ahead shells, viscous term included) with a fixed-step classical Runge-Kutta method in long
 * double, detects the section Y_2 = c Y_1 by a secant iteration on the last step, renormalizes,
 * and repeats.  It reports the growth factors G_n, the viscosity parameter eps_n and the largest
 * distance of the renormalized window from the centre, in units of the box radii.
 *
 * usage: front_steps Lnum Lden alpha eps0 nsteps KB KA center radii seed [h]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef long double R;
#define AH 14
#define MAXS 4096

static int N;              /* number of shells in the state: 0..N-1 */
static R Rk[MAXS], Ek[MAXS];

static void field(const R *y, R *dy) {
  for (int m = 0; m < N; ++m) {
    R in = (m > 0) ? y[m - 1] * y[m - 1] : 0;
    R out = (m + 1 < N) ? y[m] * y[m + 1] : 0;
    dy[m] = Rk[m] * (in - out) - Ek[m] * y[m];
  }
}

static void rk4(R *y, R h) {
  static R k1[MAXS], k2[MAXS], k3[MAXS], k4[MAXS], t[MAXS];
  field(y, k1);
  for (int m = 0; m < N; ++m) t[m] = y[m] + h / 2 * k1[m];
  field(t, k2);
  for (int m = 0; m < N; ++m) t[m] = y[m] + h / 2 * k2[m];
  field(t, k3);
  for (int m = 0; m < N; ++m) t[m] = y[m] + h * k3[m];
  field(t, k4);
  for (int m = 0; m < N; ++m) y[m] += h / 6 * (k1[m] + 2 * k2[m] + 2 * k3[m] + k4[m]);
}

int main(int argc, char **argv) {
  if (argc < 11) { fprintf(stderr, "usage: front_steps Lnum Lden alpha eps0 nsteps KB KA center radii seed [h]\n"); return 1; }
  R Lam = strtold(argv[1], 0) / strtold(argv[2], 0);
  R alpha = strtold(argv[3], 0), eps = strtold(argv[4], 0);
  int nsteps = atoi(argv[5]), KB = atoi(argv[6]), KA = atoi(argv[7]);
  const char *cf = argv[8], *rf = argv[9];
  unsigned seed = (unsigned)atoi(argv[10]);
  R h = argc > 11 ? strtold(argv[11], 0) : 1e-4L;
  R r = powl(Lam, 2.0L / 3), c = 0.5L;
  static R cen[MAXS], rad[MAXS], y[MAXS], prev[MAXS];
  int W = KB + KA + 1;
  for (int i = 0; i < W; ++i) { cen[i] = 0; rad[i] = 0; }
  { FILE *f = fopen(cf, "r"); int k; double v; while (fscanf(f, "%d %lf", &k, &v) == 2) if (k >= -KB && k <= KA) cen[k + KB] = v; fclose(f); }
  { FILE *f = fopen(rf, "r"); int k; double v; while (fscanf(f, "%d %lf", &k, &v) == 2) if (k >= -KB && k <= KA) rad[k + KB] = v; fclose(f); }
  cen[KB] = 1; cen[KB + 1] = c;
  /* initial state: front index n = KB, shells 0..KB+KA+AH; window = centre + 0.99 * (random sign) * radius */
  int n = KB; N = n + KA + AH + 1;
  srand(seed);
  for (int m = 0; m < N; ++m) y[m] = 0;
  for (int k = -KB; k <= KA; ++k) {
    R d = 0;
    if (k != 0 && k != 1) d = (seed == 0 ? 1 : ((rand() & 1) ? 1 : -1)) * 0.99L * rad[k + KB];
    if (k == KA) d = 0.99L * rad[k + KB];
    y[n + k] = cen[k + KB] + d;
    if (y[n + k] < 0) y[n + k] = 0;
  }
  R worst = 0, Gmin = 1e9, Gmax = 0, tsum = 0, logY = 0; int worstk = 0;
  for (int step = 0; step < nsteps; ++step) {
    for (int m = 0; m < N; ++m) { int k = m - n; Rk[m] = powl(r, (R)k); Ek[m] = eps * powl(Lam, 2 * alpha * k); }
    R s = 0, hprev = y[n + 2] - c * y[n + 1];
    if (!(hprev < 0)) { printf("section condition fails at start\n"); return 2; }
    for (;;) {
      memcpy(prev, y, sizeof(R) * N);
      rk4(y, h);
      R hv = y[n + 2] - c * y[n + 1];
      if (hv >= 0) {
        /* secant iteration on the substep length tau in [0, h] */
        R a = 0, b = h, fa = hprev, fb = hv;
        for (int it = 0; it < 60; ++it) {
          R tau = a - fa * (b - a) / (fb - fa);
          memcpy(y, prev, sizeof(R) * N); rk4(y, tau);
          R ft = y[n + 2] - c * y[n + 1];
          if (fabsl(ft) < 1e-30L) { a = b = tau; break; }
          if (ft < 0) { a = tau; fa = ft; } else { b = tau; fb = ft; }
          if (b - a < 1e-28L) break;
        }
        s += a;
        break;
      }
      hprev = hv; s += h;
      if (s > 5) { printf("no crossing\n"); return 3; }
    }
    R G = y[n + 1];
    if (G < Gmin) Gmin = G;
    if (G > Gmax) Gmax = G;
    /* physical time increment in units of 1/(r^n Y_n): accumulate sum_j s_j / (r^j Y_j) relative to step 0 */
    tsum += s * expl(-(R)step * logl(r) - logY);
    logY += logl(G);
    for (int m = 0; m < N; ++m) y[m] /= G;
    eps *= powl(Lam, 2 * alpha) / (r * G);
    ++n;
    if (n + KA + AH + 1 > N) { y[N] = 0; ++N; }
    for (int k = -KB; k <= KA; ++k) {
      if (k == 0 || k == 1) continue;
      R dev = fabsl(y[n + k] - cen[k + KB]) / rad[k + KB];
      if (dev > worst) { worst = dev; worstk = k; }
    }
    if (step % 25 == 0 || step == nsteps - 1)
      printf("step %3d  G=%.15Lf  eps=%.4Le  s=%.10Lf  worst box ratio %.4Lf (k=%d)\n", step, G, eps, s, worst, worstk);
  }
  printf("min G %.15Lf  max G %.15Lf  worst ratio %.4Lf  time sum (units of 1/(r^n0 Y_n0)) %.10Lf\n", Gmin, Gmax, worst, tsum);
  return 0;
}
