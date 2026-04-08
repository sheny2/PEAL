// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <map>
#include <vector>

using namespace Rcpp;
using namespace arma;

// Pattern encoding for MAR missingness
std::string make_key(const Row<double>& y_row) {
  std::string key = "";
  for(unsigned int i=0; i<y_row.n_elem; ++i) {
    key += (std::isfinite(y_row[i])) ? "1" : "0";
  }
  return key;
}

// [[Rcpp::export]]
List get_summaries_rs_cpp(const mat& X, const mat& Y, const vec& X_slope,
                          const CharacterVector& sites, const CharacterVector& patients) {
  int n = X.n_rows; int px = X.n_cols; int R = Y.n_cols;
  std::map<std::string, std::vector<int>> site_map;
  for(int i=0; i<n; ++i) site_map[as<std::string>(sites[i])].push_back(i);
  
  List Sh; 
  for(auto const& [site_name, row_indices] : site_map) {
    std::map<std::string, int> pat_to_local_id;
    int p_counter = 0;
    std::vector<int> valid_rows; std::vector<int> valid_local_pids; std::vector<std::string> row_keys;
    
    for(int idx : row_indices) {
      std::string key = make_key(Y.row(idx));
      if(key.find('1') == std::string::npos) continue; 
      
      std::string pid = as<std::string>(patients[idx]);
      if(pat_to_local_id.find(pid) == pat_to_local_id.end()) pat_to_local_id[pid] = p_counter++;
      valid_rows.push_back(idx);
      valid_local_pids.push_back(pat_to_local_id[pid]);
      row_keys.push_back(key);
    }
    if(valid_rows.empty()) continue;
    
    int pz = 1 + 2 * p_counter; // Site Int + (Pat Int + Pat Slope)
    List site_list;
    std::map<std::string, std::vector<int>> pattern_map;
    for(size_t i=0; i<row_keys.size(); ++i) pattern_map[row_keys[i]].push_back(i);
    
    for(auto const& [key, rel_indices] : pattern_map) {
      int Nh = rel_indices.size();
      mat Xs(Nh, px), Ys_full(Nh, R); uvec pids(Nh); vec slope_vals(Nh);
      std::vector<uword> temp_idx;
      for(size_t c=0; c<key.length(); ++c) if(key[c] == '1') temp_idx.push_back(c);
      uvec o = conv_to<uvec>::from(temp_idx);
      
      for(int i=0; i<Nh; ++i) {
        int o_idx = valid_rows[rel_indices[i]];
        Xs.row(i) = X.row(o_idx); Ys_full.row(i) = Y.row(o_idx);
        pids[i] = valid_local_pids[rel_indices[i]]; slope_vals[i] = X_slope[o_idx];
      }
      mat Ys = Ys_full.cols(o);
      mat ShX = Xs.t() * Xs; mat ShXZ = zeros(px, pz); mat ShZ = zeros(pz, pz); mat ShZY = zeros(pz, o.n_elem);
      
      ShXZ.col(0) = sum(Xs, 0).t(); ShZ(0,0) = Nh; ShZY.row(0) = sum(Ys, 0);
      
      for(int i=0; i<Nh; ++i) {
        int pi = 1 + pids[i]; 
        int ps = 1 + p_counter + pids[i]; // Stacked indexing
        double v = slope_vals[i];
        ShXZ.col(pi) += Xs.row(i).t(); ShXZ.col(ps) += Xs.row(i).t() * v;
        ShZ(0,pi) += 1; ShZ(pi,0) += 1; ShZ(0,ps) += v; ShZ(ps,0) += v;
        ShZ(pi,pi) += 1; ShZ(pi,ps) += v; ShZ(ps,pi) += v; ShZ(ps,ps) += v*v;
        ShZY.row(pi) += Ys.row(i); ShZY.row(ps) += Ys.row(i) * v;
      }
      site_list[key] = List::create(Named("key")=key, Named("s")=o.n_elem, Named("idx_outcomes")=o+1,
                                    Named("ShX")=ShX, Named("ShXZ")=ShXZ, Named("ShXY")=Xs.t()*Ys,
                                          Named("ShZ")=ShZ, Named("ShZY")=ShZY, Named("ShYY")=Ys.t()*Ys,
                                                Named("Nh")=Nh, Named("pzh")=pz, Named("mh")=p_counter);
    }
    Sh[site_name] = site_list;
  }
  return Sh;
}

// [[Rcpp::export]]
List compute_site_stats_rs_cpp(List ShPat, vec par, int R, int px, mat Corr_f) {
  double lp1 = 0, lp2 = 0; mat b1 = zeros(R*px, R*px); vec b2 = zeros(R*px);
  double Ntot = 0; int K = ShPat.size();
  
  for (int h = 0; h < K; ++h) {
    List S_h = ShPat[h]; if (S_h.size() == 0) continue;
    double tu = par[0], tvi = par[1 + h], tvs = par[1 + K + h]; // Offset indexing
    int pzh = as<List>(S_h[0])["pzh"]; int mh = as<List>(S_h[0])["mh"];
    
    vec v0 = zeros(pzh); v0[0] = tu; 
    for(int i=1; i<=mh; ++i) v0[i] = tvi; 
    for(int i=mh+1; i<pzh; ++i) v0[i] = tvs; 
    
    mat Th_inv = kron(eye(R, R), diagmat(1.0 / v0));
    mat Sxx = zeros(R*px, R*px), Sxz = zeros(R*px, R*pzh), Szz = zeros(R*pzh, R*pzh);
    vec sxy = zeros(R*px), szy = zeros(R*pzh); double yty = 0;
    
    for (int k = 0; k < S_h.size(); ++k) {
      List S = S_h[k]; uvec o = as<uvec>(S["idx_outcomes"]) - 1; 
      mat Ri = inv_sympd(Corr_f.submat(o, o) + eye(o.n_elem,o.n_elem)*1e-9);
      lp1 += as<int>(S["Nh"]) * log_det_sympd(Corr_f.submat(o, o) + eye(o.n_elem,o.n_elem)*1e-9);
      
      mat Ri_emb = zeros(R, R); Ri_emb.submat(o, o) = Ri;
      
      // OPTIMIZATION: Zero-copy pointer mapping for matrices
      NumericMatrix ShX_R = as<NumericMatrix>(S["ShX"]);
      mat ShX_arma(ShX_R.begin(), ShX_R.nrow(), ShX_R.ncol(), false);
      
      NumericMatrix ShXZ_R = as<NumericMatrix>(S["ShXZ"]);
      mat ShXZ_arma(ShXZ_R.begin(), ShXZ_R.nrow(), ShXZ_R.ncol(), false);
      
      NumericMatrix ShZ_R = as<NumericMatrix>(S["ShZ"]);
      mat ShZ_arma(ShZ_R.begin(), ShZ_R.nrow(), ShZ_R.ncol(), false);
      
      NumericMatrix ShXY_R = as<NumericMatrix>(S["ShXY"]);
      mat ShXY_arma(ShXY_R.begin(), ShXY_R.nrow(), ShXY_R.ncol(), false);
      
      NumericMatrix ShZY_R = as<NumericMatrix>(S["ShZY"]);
      mat ShZY_arma(ShZY_R.begin(), ShZY_R.nrow(), ShZY_R.ncol(), false);
      
      NumericMatrix ShYY_R = as<NumericMatrix>(S["ShYY"]);
      mat ShYY_arma(ShYY_R.begin(), ShYY_R.nrow(), ShYY_R.ncol(), false);
      
      Sxx += kron(Ri_emb, ShX_arma); 
      Sxz += kron(Ri_emb, ShXZ_arma); 
      Szz += kron(Ri_emb, ShZ_arma);
      
      mat ShXY_f = zeros(px, R); ShXY_f.cols(o) = ShXY_arma * Ri;
      mat ShZY_f = zeros(pzh, R); ShZY_f.cols(o) = ShZY_arma * Ri;
      sxy += vectorise(ShXY_f); szy += vectorise(ShZY_f); yty += trace(Ri * ShYY_arma);
      Ntot += (double)as<int>(S["Nh"]) * o.n_elem;
    }
    mat Ah = Th_inv + Szz; mat Ah_i = inv_sympd(Ah);
    lp1 += R * sum(log(v0)) + log_det_sympd(Ah);
    mat proj = Sxz * Ah_i; b1 += Sxx - proj * Sxz.t(); b2 += sxy - proj * szy;
    lp2 += yty - as_scalar(szy.t() * Ah_i * szy);
  }
  return List::create(Named("lp1")=lp1, Named("lp2")=lp2, Named("b1")=b1, Named("b2")=b2, Named("Ntot")=Ntot);
}

// [[Rcpp::export]]
vec compute_gradient_rs_cpp(List ShPat, vec par, int R, int px, mat Corr_f, vec beta, double s2) {
  int K = ShPat.size(); vec grad = zeros(par.size());
  for (int h = 0; h < K; ++h) {
    List S_h = ShPat[h]; if (S_h.size() == 0) continue;
    double tu = par[0], tvi = par[1 + h], tvs = par[1 + K + h];
    int pzh = as<List>(S_h[0])["pzh"]; int mh = as<List>(S_h[0])["mh"];
    vec v0 = zeros(pzh); v0[0] = tu; for(int i=1; i<=mh; ++i) v0[i] = tvi; for(int i=mh+1; i<pzh; ++i) v0[i] = tvs;
    
    mat Th_inv = kron(eye(R, R), diagmat(1.0 / v0));
    mat Sxz = zeros(R * px, R * pzh), Szz = zeros(R * pzh, R * pzh); vec szy = zeros(R * pzh);
    for (int k = 0; k < S_h.size(); ++k) {
      List S = S_h[k]; uvec o = as<uvec>(S["idx_outcomes"]) - 1; 
      mat Ri = inv_sympd(Corr_f.submat(o, o) + eye(o.n_elem,o.n_elem)*1e-9); 
      mat Ri_emb = zeros(R, R); Ri_emb.submat(o, o) = Ri;
      
      // OPTIMIZATION: Zero-copy pointer mapping for matrices
      NumericMatrix ShXZ_R = as<NumericMatrix>(S["ShXZ"]);
      mat ShXZ_arma(ShXZ_R.begin(), ShXZ_R.nrow(), ShXZ_R.ncol(), false);
      
      NumericMatrix ShZ_R = as<NumericMatrix>(S["ShZ"]);
      mat ShZ_arma(ShZ_R.begin(), ShZ_R.nrow(), ShZ_R.ncol(), false);
      
      NumericMatrix ShZY_R = as<NumericMatrix>(S["ShZY"]);
      mat ShZY_arma(ShZY_R.begin(), ShZY_R.nrow(), ShZY_R.ncol(), false);
      
      Sxz += kron(Ri_emb, ShXZ_arma); 
      Szz += kron(Ri_emb, ShZ_arma);
      mat ShZY_f = zeros(pzh, R); ShZY_f.cols(o) = ShZY_arma * Ri; 
      szy += vectorise(ShZY_f);
    }
    mat Ah_i = inv_sympd(Th_inv + Szz); vec w = Ah_i * (szy - Sxz.t() * beta); vec d_Ah = Ah_i.diag();
    
    auto calc_g = [&](int start, int end, double theta, int g_idx) {
      double tr = 0, dat = 0;
      for(int r=0; r<R; ++r) { for(int i=start; i<=end; ++i) { int idx = r*pzh + i; tr += (1.0/theta) - (1.0/(theta*theta))*d_Ah[idx]; dat += w[idx]*w[idx]; } }
      grad[g_idx] += 0.5 * ( (1.0/(s2*theta*theta))*dat - tr );
    };
    calc_g(0, 0, tu, 0); calc_g(1, mh, tvi, 1+h); calc_g(mh+1, pzh-1, tvs, 1+K+h);
  }
  return grad;
}