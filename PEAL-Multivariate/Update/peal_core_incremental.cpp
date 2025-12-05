// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <map>
#include <vector>

using namespace Rcpp;
using namespace arma;

std::string make_key_v2(const Row<double>& y_row) {
  std::string key = "";
  for(unsigned int i=0; i<y_row.n_elem; ++i) key += (std::isfinite(y_row[i])) ? "1" : "0";
  return key;
}

// ----------------------------------------------------------------
// 1. CONTROLLED SUMMARY GENERATOR
//    Uses a PRE-DEFINED patient mapping (registry) to ensure 
//    indices match the old model.
// ----------------------------------------------------------------
// [[Rcpp::export]]
List get_summaries_controlled_cpp(const mat& X, const mat& Y, 
                                  const CharacterVector& sites, 
                                  const CharacterVector& patients,
                                  List patient_registry) { 
  // patient_registry: List where name is site, value is CharacterVector of PatientIDs
  
  int px = X.n_cols;
  int R = Y.n_cols;
  
  // Convert Registry to fast C++ Maps
  std::map<std::string, std::map<std::string, int>> global_map;
  std::map<std::string, int> site_pzh_map;
  
  CharacterVector reg_sites = patient_registry.names();
  for(int i=0; i<reg_sites.size(); ++i) {
    std::string s = as<std::string>(reg_sites[i]);
    CharacterVector pids = patient_registry[i];
    
    // pzh = 1 (site intercept) + number of patients
    site_pzh_map[s] = 1 + pids.size();
    
    for(int j=0; j<pids.size(); ++j) {
      std::string pid = as<std::string>(pids[j]);
      // Patient j is at index j+1 (0 is site intercept)
      global_map[s][pid] = j + 1;
    }
  }
  
  // Group Data by Site
  std::map<std::string, std::vector<int> > site_row_map;
  for(int i=0; i<X.n_rows; ++i) {
    site_row_map[as<std::string>(sites[i])].push_back(i);
  }
  
  List Sh; 
  
  for(auto const& [site_name, row_indices] : site_row_map) {
    
    // If this site isn't in registry (Brand new site in new data),
    // We skip or handle separately. Ideally registry is updated in R before calling this.
    if(global_map.find(site_name) == global_map.end()) continue;
    
    auto& local_pat_map = global_map[site_name];
    int pz = site_pzh_map[site_name];
    
    List site_list;
    std::map<std::string, std::vector<int>> pattern_map;
    std::vector<std::string> row_keys;
    std::vector<int> valid_rows;
    
    for(int idx : row_indices) {
      std::string key = make_key_v2(Y.row(idx));
      bool is_empty = true;
      for(char c : key) if(c=='1') is_empty = false;
      if(is_empty) continue;
      
      row_keys.push_back(key);
      valid_rows.push_back(idx);
    }
    
    for(size_t i=0; i<row_keys.size(); ++i) {
      pattern_map[row_keys[i]].push_back(i);
    }
    
    for(auto const& [key, rel_indices] : pattern_map) {
      int Nh = rel_indices.size();
      mat Xs(Nh, px);
      mat Ys_full(Nh, R); 
      uvec pids(Nh); 
      
      std::vector<uword> temp_idx;
      for(size_t c=0; c<key.length(); ++c) if(key[c] == '1') temp_idx.push_back(c);
      uvec idx_outcomes = conv_to<uvec>::from(temp_idx);
      int s = idx_outcomes.n_elem;
      
      for(int i=0; i<Nh; ++i) {
        int original_idx = valid_rows[rel_indices[i]];
        Xs.row(i) = X.row(original_idx);
        Ys_full.row(i) = Y.row(original_idx);
        
        std::string pid = as<std::string>(patients[original_idx]);
        // LOOKUP ID from Global Map
        pids[i] = local_pat_map[pid]; 
      }
      mat Ys = Ys_full.cols(idx_outcomes);
      
      // --- Compute Cross Products ---
      mat ShX = Xs.t() * Xs;
      mat ShXY = Xs.t() * Ys;
      mat ShYY = Ys.t() * Ys;
      
      // Z components depend on pz (Full size from registry)
      mat ShXZ = zeros(px, pz);
      ShXZ.col(0) = sum(Xs, 0).t(); 
      for(int i=0; i<Nh; ++i) ShXZ.col(pids[i]) += Xs.row(i).t();
      
      mat ShZ = zeros(pz, pz);
      ShZ(0,0) = Nh; 
      for(int i=0; i<Nh; ++i) {
        ShZ(0, pids[i]) += 1;
        ShZ(pids[i], 0) += 1;
        ShZ(pids[i], pids[i]) += 1;
      }
      
      mat ShZY = zeros(pz, s);
      ShZY.row(0) = sum(Ys, 0); 
      for(int i=0; i<Nh; ++i) ShZY.row(pids[i]) += Ys.row(i);
      
      uvec idx_outcomes_r = idx_outcomes + 1; 
      
      List pat_list = List::create(
        Named("key") = key, Named("s") = s, Named("idx_outcomes") = idx_outcomes_r,
        Named("ShX") = ShX, Named("ShXZ") = ShXZ, Named("ShXY") = ShXY,
              Named("ShZ") = ShZ, Named("ShZY") = ShZY, Named("ShYY") = ShYY,
                    Named("Nh") = Nh, Named("pzh") = pz
      );
      site_list[key] = pat_list;
    }
    Sh[site_name] = site_list;
  }
  return Sh;
}

// ----------------------------------------------------------------
// 2. PAD SUMMARIES
//    Resizes "Old" matrices to match "New" dimensions 
//    (fills new rows/cols with zeros)
// ----------------------------------------------------------------
// [[Rcpp::export]]
List pad_summaries_cpp(List ShPat, List target_sizes_list) {
  // target_sizes_list: Name=Site, Value=Integer (New pzh)
  
  CharacterVector sites = ShPat.names();
  List Sh_New = clone(ShPat);
  
  // Map for fast lookup
  std::map<std::string, int> size_map;
  CharacterVector tn = target_sizes_list.names();
  for(int i=0; i<tn.size(); ++i) size_map[as<std::string>(tn[i])] = as<int>(target_sizes_list[i]);
  
  for(int h=0; h<sites.size(); ++h) {
    std::string site = as<std::string>(sites[h]);
    if(size_map.find(site) == size_map.end()) continue;
    
    int target_pz = size_map[site];
    List site_list = Sh_New[h];
    
    for(int k=0; k<site_list.size(); ++k) {
      List S = site_list[k];
      int old_pz = S["pzh"];
      
      if(old_pz >= target_pz) continue; // No padding needed
      
      // PAD ShZ (pz x pz)
      mat ShZ = as<mat>(S["ShZ"]);
      mat ShZ_new = zeros(target_pz, target_pz);
      ShZ_new.submat(0, 0, old_pz-1, old_pz-1) = ShZ;
      
      // PAD ShXZ (px x pz)
      mat ShXZ = as<mat>(S["ShXZ"]);
      mat ShXZ_new = zeros(ShXZ.n_rows, target_pz);
      ShXZ_new.cols(0, old_pz-1) = ShXZ;
      
      // PAD ShZY (pz x s)
      mat ShZY = as<mat>(S["ShZY"]);
      mat ShZY_new = zeros(target_pz, ShZY.n_cols);
      ShZY_new.rows(0, old_pz-1) = ShZY;
      
      // Update List
      S["ShZ"] = ShZ_new;
      S["ShXZ"] = ShXZ_new;
      S["ShZY"] = ShZY_new;
      S["pzh"] = target_pz;
      
      site_list[k] = S;
    }
    Sh_New[h] = site_list;
  }
  return Sh_New;
}

// ----------------------------------------------------------------
// 3. PROFILE LIKELIHOOD (Optimized Woodbury)
// ----------------------------------------------------------------
// [[Rcpp::export]]
List compute_site_stats_cpp(List ShPat, vec par, int R, int px, std::string corstr, double rho_fixed, mat Corr_full, bool estimate_rho) {
  double lpterm1 = 0, lpterm2 = 0;
  mat bterm1 = zeros(R * px, R * px);
  vec bterm2 = zeros(R * px);
  double Nobs_total = 0, Nsum_rows = 0;
  CharacterVector site_names = ShPat.names();
  int K = site_names.size();
  double theta_u = par[0];
  
  for (int h = 0; h < K; ++h) {
    std::string site = as<std::string>(site_names[h]);
    List S_h = ShPat[site];
    if (S_h.size() == 0) continue;
    double theta_vh = par[1 + h];
    
    List first_pat = S_h[0];
    int pzh = first_pat["pzh"];
    vec v0_diag = zeros(pzh); v0_diag[0] = theta_u;
    for(int i=1; i<pzh; ++i) v0_diag[i] = theta_vh;
    double logdet_Theta_h = R * sum(log(v0_diag));
    mat V0_inv = diagmat(1.0 / v0_diag);
    mat Theta_h_inv = kron(eye(R, R), V0_inv);
    
    mat SxxR_sum = zeros(R * px, R * px);
    mat Sxz_sum  = zeros(R * px, R * pzh);
    mat Szz_sum  = zeros(R * pzh, R * pzh);
    vec sxy_sum  = zeros(R * px);
    vec szy_sum  = zeros(R * pzh);
    double ytildeY_sum = 0, logdet_Rcorr_sum = 0, Nh_site = 0;
    
    for (int k = 0; k < S_h.size(); ++k) {
      List S = S_h[k];
      uvec o = as<uvec>(S["idx_outcomes"]) - 1; 
      int s = S["s"], Nh = S["Nh"];
      Nh_site += Nh; Nobs_total += Nh * s;
      
      mat ShX = as<mat>(S["ShX"]), ShXZ = as<mat>(S["ShXZ"]), ShXY = as<mat>(S["ShXY"]);
      mat ShZ = as<mat>(S["ShZ"]), ShZY = as<mat>(S["ShZY"]), ShYY = as<mat>(S["ShYY"]);
      
      mat Rinve_s; double ld_s = 0;
      if (corstr == "independence") { Rinve_s = eye(s, s); } 
      else if (corstr == "exchangeable") {
        double rho_curr = (estimate_rho) ? par[par.size()-1] : rho_fixed;
        double a = 1.0/(1.0-rho_curr); double b = rho_curr/(1.0-rho_curr+rho_curr*s);
        Rinve_s = a * (eye(s, s) - b * ones(s, s));
        ld_s = (s - 1) * log(1.0 - rho_curr) + log(1.0 - rho_curr + rho_curr * s);
      } else {
        mat Coo = Corr_full.submat(o, o);
        mat Uoo = chol(Coo + eye(s,s)*1e-9);
        Rinve_s = inv(trimatu(Uoo)); Rinve_s = Rinve_s * Rinve_s.t();
        ld_s = 2.0 * sum(log(Uoo.diag()));
      }
      logdet_Rcorr_sum += Nh * ld_s;
      mat Rinve_embed = zeros(R, R); Rinve_embed.submat(o, o) = Rinve_s;
      SxxR_sum += kron(Rinve_embed, ShX); Sxz_sum  += kron(Rinve_embed, ShXZ); Szz_sum  += kron(Rinve_embed, ShZ);
      mat ShXY_full = zeros(ShXY.n_rows, R); ShXY_full.cols(o) = ShXY * Rinve_s; sxy_sum += vectorise(ShXY_full);
      mat ShZY_full = zeros(ShZY.n_rows, R); ShZY_full.cols(o) = ShZY * Rinve_s; szy_sum += vectorise(ShZY_full);
      ytildeY_sum += trace(Rinve_s * ShYY);
    }
    
    mat A_h = Theta_h_inv + Szz_sum;
    mat A_h_inv; double logdet_Ah; double sign;
    bool success = log_det(logdet_Ah, sign, A_h);
    if(success) { A_h_inv = inv_sympd(A_h); } else { A_h_inv = inv(A_h + eye(size(A_h))*1e-6); log_det(logdet_Ah, sign, A_h); }
    
    lpterm1 += logdet_Rcorr_sum + logdet_Theta_h + logdet_Ah;
    mat term_xz = Sxz_sum * A_h_inv;
    bterm1 += SxxR_sum - term_xz * Sxz_sum.t();
    bterm2 += sxy_sum - term_xz * szy_sum;
    lpterm2 += ytildeY_sum - as_scalar(szy_sum.t() * A_h_inv * szy_sum);
    Nsum_rows += Nh_site;
  }
  return List::create(Named("lpterm1")=lpterm1, Named("lpterm2")=lpterm2, Named("bterm1")=bterm1, Named("bterm2")=bterm2, Named("Nobs_total")=Nobs_total, Named("Nsum_rows")=Nsum_rows);
}

// ----------------------------------------------------------------
// 4. ANALYTIC GRADIENT (Hybrid)
// ----------------------------------------------------------------
// [[Rcpp::export]]
vec compute_gradient_theta_cpp(List ShPat, vec par, int R, int px, std::string corstr, double rho_fixed, mat Corr_full, bool estimate_rho, vec beta, double sigma2) { 
  CharacterVector site_names = ShPat.names();
  int K = site_names.size();
  vec grad = zeros(1 + K);
  double theta_u = par[0];
  
  for (int h = 0; h < K; ++h) {
    std::string site = as<std::string>(site_names[h]);
    List S_h = ShPat[site];
    if (S_h.size() == 0) continue;
    double theta_vh = par[1 + h];
    List first_pat = S_h[0];
    int pzh = first_pat["pzh"];
    vec v0_diag = zeros(pzh); v0_diag[0] = theta_u;
    for(int i=1; i<pzh; ++i) v0_diag[i] = theta_vh;
    mat V0_inv = diagmat(1.0 / v0_diag);
    mat Theta_h_inv = kron(eye(R, R), V0_inv);
    mat Szz_sum = zeros(R * pzh, R * pzh); vec szy_sum = zeros(R * pzh); mat Sxz_sum = zeros(R * px, R * pzh); 
    
    for (int k = 0; k < S_h.size(); ++k) {
      List S = S_h[k];
      uvec o = as<uvec>(S["idx_outcomes"]) - 1; int s = S["s"];
      mat Rinve_s;
      if (corstr == "exchangeable") {
        double rho_curr = (estimate_rho) ? par[par.size()-1] : rho_fixed;
        double a = 1.0/(1.0-rho_curr); double b = rho_curr/(1.0-rho_curr+rho_curr*s);
        Rinve_s = a * (eye(s,s) - b*ones(s,s));
      } else if (corstr == "independence") { Rinve_s = eye(s, s); }
      else { mat Coo = Corr_full.submat(o, o); Rinve_s = inv(symmatu(Coo + eye(s,s)*1e-9)); }
      
      mat Rinve_embed = zeros(R, R); Rinve_embed.submat(o, o) = Rinve_s;
      mat ShZ = as<mat>(S["ShZ"]); mat ShZY = as<mat>(S["ShZY"]); mat ShXZ = as<mat>(S["ShXZ"]);
      Szz_sum += kron(Rinve_embed, ShZ); Sxz_sum += kron(Rinve_embed, ShXZ);
      mat ShZY_full = zeros(ShZY.n_rows, R); ShZY_full.cols(o) = ShZY * Rinve_s; szy_sum += vectorise(ShZY_full);
    }
    
    mat A_h = Theta_h_inv + Szz_sum; mat A_h_inv;
    try { A_h_inv = inv_sympd(A_h); } catch(...) { A_h_inv = inv(A_h + eye(size(A_h))*1e-8); }
    vec v_ZR = szy_sum - Sxz_sum.t() * beta; vec w = A_h_inv * v_ZR; vec diag_A_inv = A_h_inv.diag();
    
    double tr_term_u = 0, data_term_u = 0;
    for(int r=0; r<R; ++r) { int idx = r * pzh; tr_term_u += (1.0/theta_u) - (1.0/(theta_u*theta_u)) * diag_A_inv[idx]; data_term_u += w[idx] * w[idx]; }
    data_term_u = (1.0 / (sigma2 * theta_u * theta_u)) * data_term_u; grad[0] += 0.5 * (data_term_u - tr_term_u);
    
    double tr_term_v = 0, data_term_v = 0;
    for(int r=0; r<R; ++r) { for(int p=1; p<pzh; ++p) { int idx = r * pzh + p; tr_term_v += (1.0/theta_vh) - (1.0/(theta_vh*theta_vh)) * diag_A_inv[idx]; data_term_v += w[idx] * w[idx]; }}
    data_term_v = (1.0 / (sigma2 * theta_vh * theta_vh)) * data_term_v; grad[1 + h] = 0.5 * (data_term_v - tr_term_v);
  }
  return grad;
}