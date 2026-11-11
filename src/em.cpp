#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

static const double LOG2PI = 1.8378770664093453;

// EM iterations for the diagonal Gaussian mixture with known groups, equality
// constraints (parameters sharing a label id are pooled) and fixed parameters
// (id 0), when profile membership has no covariates and the membership
// logits are free (closed-form proportion updates). Arrays mu, s2 and
// var_floor are p x K x G in column-major order (variable fastest); pi is
// K x G. Returns the updated parameters, the log-likelihood at the returned
// parameters, the number of iterations and the convergence flag, mirroring
// the pure-R implementation in R/em.R.
// [[Rcpp::export]]
List em_core_cpp(const NumericMatrix& Y, const NumericMatrix& O, const IntegerVector& g,
                 int G, NumericVector mu, NumericVector s2, NumericMatrix pi,
                 bool equal_prop, const IntegerVector& mean_id, const IntegerVector& var_id,
                 int n_mean, int n_var, const NumericVector& var_floor, int maxiter,
                 double tol, bool check) {
  const int N = Y.nrow(), p = Y.ncol(), K = pi.nrow();
  const int PK = p * K, PKG = PK * G;
  NumericVector mu_c = clone(mu), s2_c = clone(s2);
  NumericMatrix pi_c = clone(pi);
  const double* y = Y.begin();     // column-major: y[i + N * j]
  const double* o = O.begin();
  const int* gp = g.begin();
  double* muv = mu_c.begin();
  double* s2v = s2_c.begin();
  const double* floorv = var_floor.begin();
  const int* mid = mean_id.begin();
  const int* vid = var_id.begin();
  std::vector<double> tau(static_cast<size_t>(N) * K);   // tau[i + N * k]
  std::vector<double> lconst(PKG), inv(PKG), lp(K), logpi(K * G);
  std::vector<double> A(PKG), B(PKG), C(PKG);
  const int n_pool = std::max(std::max(n_mean, n_var), 1);
  std::vector<double> num(n_pool), den(n_pool);
  std::vector<int> ng(G, 0);
  for (int i = 0; i < N; ++i) ng[gp[i] - 1]++;

  auto estep = [&]() -> double {
    for (int idx = 0; idx < PKG; ++idx) {
      inv[idx] = 1.0 / s2v[idx];
      lconst[idx] = -0.5 * (LOG2PI + std::log(s2v[idx]));
    }
    for (int gg = 0; gg < G; ++gg)
      for (int k = 0; k < K; ++k) logpi[k + K * gg] = std::log(pi_c(k, gg));
    double total = 0.0;
    for (int i = 0; i < N; ++i) {
      const int gi = gp[i] - 1;
      double m = R_NegInf;
      for (int k = 0; k < K; ++k) {
        double v = logpi[k + K * gi];
        const int base = p * k + PK * gi;
        for (int j = 0; j < p; ++j) {
          const int ij = i + N * j;
          if (o[ij] != 0.0) {
            const double d = y[ij] - muv[base + j];
            v += lconst[base + j] - 0.5 * d * d * inv[base + j];
          }
        }
        lp[k] = v;
        if (v > m) m = v;
      }
      double s = 0.0;
      for (int k = 0; k < K; ++k) s += std::exp(lp[k] - m);
      const double lse = m + std::log(s);
      total += lse;
      for (int k = 0; k < K; ++k) tau[i + N * k] = std::exp(lp[k] - lse);
    }
    return total;
  };

  auto mstep = [&]() {
    std::fill(A.begin(), A.end(), 0.0);
    std::fill(B.begin(), B.end(), 0.0);
    std::fill(C.begin(), C.end(), 0.0);
    for (int i = 0; i < N; ++i) {
      const int gi = gp[i] - 1;
      for (int k = 0; k < K; ++k) {
        const double t = tau[i + N * k];
        const int base = p * k + PK * gi;
        for (int j = 0; j < p; ++j) {
          const int ij = i + N * j;
          if (o[ij] != 0.0) {
            const double yy = y[ij];
            A[base + j] += t * yy;
            B[base + j] += t;
            C[base + j] += t * yy * yy;
          }
        }
      }
    }
    if (n_mean > 0) {
      std::fill(num.begin(), num.begin() + n_mean, 0.0);
      std::fill(den.begin(), den.begin() + n_mean, 0.0);
      for (int idx = 0; idx < PKG; ++idx) {
        const int id = mid[idx];
        if (id > 0) {
          const double w = 1.0 / s2v[idx];
          num[id - 1] += A[idx] * w;
          den[id - 1] += B[idx] * w;
        }
      }
      for (int idx = 0; idx < PKG; ++idx) {
        const int id = mid[idx];
        if (id > 0 && den[id - 1] > 0.0) muv[idx] = num[id - 1] / den[id - 1];
      }
    }
    if (n_var > 0) {
      std::fill(num.begin(), num.begin() + n_var, 0.0);
      std::fill(den.begin(), den.begin() + n_var, 0.0);
      for (int idx = 0; idx < PKG; ++idx) {
        const int id = vid[idx];
        if (id > 0) {
          const double m = muv[idx];
          num[id - 1] += C[idx] - 2.0 * m * A[idx] + m * m * B[idx];
          den[id - 1] += B[idx];
        }
      }
      for (int idx = 0; idx < PKG; ++idx) {
        const int id = vid[idx];
        if (id > 0 && den[id - 1] > 0.0) {
          double v = num[id - 1] / den[id - 1];
          if (v < floorv[idx]) v = floorv[idx];
          s2v[idx] = v;
        }
      }
    }
    if (K > 1) {
      if (equal_prop) {
        std::vector<double> pk(K, 0.0);
        for (int k = 0; k < K; ++k)
          for (int i = 0; i < N; ++i) pk[k] += tau[i + N * k];
        double tot = 0.0;
        for (int k = 0; k < K; ++k) { pk[k] = std::max(pk[k] / N, 1e-12); tot += pk[k]; }
        for (int gg = 0; gg < G; ++gg)
          for (int k = 0; k < K; ++k) pi_c(k, gg) = pk[k] / tot;
      } else {
        std::vector<double> pk(K * G, 0.0);
        for (int k = 0; k < K; ++k)
          for (int i = 0; i < N; ++i) pk[k + K * (gp[i] - 1)] += tau[i + N * k];
        for (int gg = 0; gg < G; ++gg) {
          double tot = 0.0;
          for (int k = 0; k < K; ++k) {
            pk[k + K * gg] = std::max(pk[k + K * gg] / ng[gg], 1e-12);
            tot += pk[k + K * gg];
          }
          for (int k = 0; k < K; ++k) pi_c(k, gg) = pk[k + K * gg] / tot;
        }
      }
    }
  };

  double ll = R_NegInf, ll_old = R_NegInf;
  bool converged = false;
  int iter = 0;
  for (int it = 1; it <= maxiter; ++it) {
    ll = estep();
    iter = it;
    if (!R_finite(ll)) break;
    if (check && std::fabs(ll - ll_old) < tol) { converged = true; break; }
    ll_old = ll;
    mstep();
  }
  if (!converged) ll = estep();
  return List::create(_["mu"] = mu_c, _["s2"] = s2_c, _["pi"] = pi_c, _["ll"] = ll,
                      _["iter"] = iter, _["converged"] = converged);
}
