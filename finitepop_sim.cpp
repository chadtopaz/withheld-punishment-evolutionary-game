#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Unordered pairs k = 0, 1, 2 are (i, j) = (D, P), (D, A), (P, A) with i < j.
static const int PI_[3] = {0, 0, 1};
static const int PJ_[3] = {1, 2, 2};
// Position of the ordered transitions i -> j and j -> i in the R table `trans`
// (rows 1->2, 1->3, 2->1, 2->3, 3->1, 3->2; 0-based here).
static const int FWD[3] = {0, 1, 3};
static const int BWD[3] = {2, 4, 5};

struct Rates {
  double dM[3][3];   // dM[k][m] = M(j, m) - M(i, m), so that pi_j - pi_i = dM[k] . x
  double om, mu3, sel;
  Rates(const NumericMatrix& M, double mu, double w, double dpi_max) {
    for (int k = 0; k < 3; k++) for (int m = 0; m < 3; m++) dM[k][m] = M(PJ_[k], m) - M(PI_[k], m);
    om = 1.0 - mu; mu3 = mu / 3.0; sel = 0.5 * w / dpi_max;
  }
  // probability per step of a change within pair k
  inline double pair_weight(int k, const double* x) const {
    return om * x[PI_[k]] * x[PJ_[k]] + mu3 * (x[PI_[k]] + x[PJ_[k]]);
  }
  // probability per step that a type-i individual becomes type j (pair k, i < j)
  inline double forward_rate(int k, const double* x) const {
    double dpi = dM[k][0] * x[0] + dM[k][1] * x[1] + dM[k][2] * x[2];
    double p = 0.5 + sel * dpi;
    return x[PI_[k]] * (om * x[PJ_[k]] * p + mu3);
  }
};

// [[Rcpp::export]]
NumericVector event_probs_cpp(NumericMatrix M, int N, double mu, double w, double dpi_max, IntegerVector n) {
  Rates rt(M, mu, w, dpi_max);
  double x[3]; for (int i = 0; i < 3; i++) x[i] = (double)n[i] / (double)N;
  NumericVector out(6);
  for (int k = 0; k < 3; k++) {
    double R = rt.pair_weight(k, x), f = rt.forward_rate(k, x);
    out[FWD[k]] = f; out[BWD[k]] = R - f;
  }
  return out;
}

// Event-driven local update process with reintroduction for a three-strategy game.
// [[Rcpp::export]]
List simulate_chain_cpp(NumericMatrix M, int N, double mu, double w, double dpi_max,
                        double T_gen, IntegerVector n0, bool stop_at_vertex, bool track_episodes,
                        double lo, double hi, double record_dt, double max_events) {
  Rates rt(M, mu, w, dpi_max);
  int n[3] = {(int)n0[0], (int)n0[1], (int)n0[2]};
  const double Nd = (double)N;
  const double max_steps = T_gen * Nd;         // horizon in elementary steps
  double t_steps = 0.0, n_events = 0.0;
  double t_ext = NA_REAL; int lost_first = NA_INTEGER;
  bool armed = false; int n_ep = 0; double time_hi = 0.0;
  const bool recording = record_dt > 0;
  const double rec_step = record_dt * Nd;
  double next_rec = 0.0;
  std::vector<double> rec_t, rec_D, rec_P, rec_A;
  double x[3], R[3];
  while (true) {
    for (int i = 0; i < 3; i++) x[i] = (double)n[i] / Nd;
    for (int k = 0; k < 3; k++) R[k] = rt.pair_weight(k, x);
    double q = R[0] + R[1] + R[2];
    double hold;
    if (q <= 0.0) hold = R_PosInf;
    else { double u = R::runif(0.0, 1.0); hold = std::floor(std::log(u) / std::log1p(-q)) + 1.0; }
    double t_next = t_steps + hold;            // step at which the next change occurs
    if (recording) {                           // the state is constant on [t_steps, t_next)
      while (next_rec < t_next && next_rec <= max_steps) {
        rec_t.push_back(next_rec / Nd); rec_D.push_back(x[0]); rec_P.push_back(x[1]); rec_A.push_back(x[2]);
        next_rec += rec_step;
      }
    }
    if (track_episodes && x[2] > hi) time_hi += std::min(hold, max_steps - t_steps);
    if (t_next > max_steps) { t_steps = max_steps; break; }   // no further change within the horizon
    t_steps = t_next; n_events += 1.0;
    // draw the pair with probability R_k / q, then the direction; v is uniform on [0, q)
    double v = R::runif(0.0, 1.0) * q;
    int k; double cum;
    if (v < R[0]) { k = 0; cum = 0.0; }
    else if (v < R[0] + R[1]) { k = 1; cum = R[0]; }
    else { k = 2; cum = R[0] + R[1]; }
    double v_res = v - cum;
    if (R[k] <= 0.0) {                         // reachable only through rounding at the top of [0, q)
      k = 0; if (R[1] > R[k]) k = 1; if (R[2] > R[k]) k = 2;
      v_res = 0.0;
    }
    if (v_res < 0.0) v_res = 0.0;
    const int i = PI_[k], j = PJ_[k];
    bool forward = v_res < rt.forward_rate(k, x);
    if (forward && n[i] == 0) forward = false;         // guards against rounding at the boundary
    else if (!forward && n[j] == 0) forward = true;
    if (forward) { n[i] -= 1; n[j] += 1; } else { n[j] -= 1; n[i] += 1; }
    if (n[0] < 0 || n[1] < 0 || n[2] < 0 || n[0] > N || n[1] > N || n[2] > N) stop("count out of range");
    if (lost_first == NA_INTEGER) {
      for (int m = 0; m < 3; m++) if (n[m] == 0) { t_ext = t_steps / Nd; lost_first = m + 1; break; }
    }
    if (track_episodes) {
      double xA = (double)n[2] / Nd;
      if (!armed && xA < lo) armed = true;
      else if (armed && xA > hi) { armed = false; n_ep += 1; }
    }
    if (stop_at_vertex && (n[0] == N || n[1] == N || n[2] == N)) break;
    if (max_events > 0 && n_events >= max_events) break;
  }
  const bool at_vertex = (n[0] == N || n[1] == N || n[2] == N);
  IntegerVector fin = IntegerVector::create(n[0], n[1], n[2]);
  RObject rec = R_NilValue;
  if (recording) {
    int m = rec_t.size(); NumericMatrix Rm(m, 4);
    for (int r = 0; r < m; r++) { Rm(r, 0) = rec_t[r]; Rm(r, 1) = rec_D[r]; Rm(r, 2) = rec_P[r]; Rm(r, 3) = rec_A[r]; }
    colnames(Rm) = CharacterVector::create("gen", "xD", "xP", "xA");
    rec = Rm;
  }
  return List::create(Named("t_ext") = t_ext, Named("lost_first") = lost_first,
                      Named("at_vertex") = at_vertex, Named("t_vertex") = at_vertex ? t_steps / Nd : NA_REAL,
                      Named("final") = fin, Named("n_ep") = n_ep, Named("time_hi") = time_hi,
                      Named("t_run") = t_steps / Nd, Named("n_events") = n_events, Named("rec") = rec);
}
