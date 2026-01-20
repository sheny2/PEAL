// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <map>
#include <vector>

using namespace Rcpp;
using namespace arma;

// ----------------------------------------------------------------
// HELPER: Generate Pattern Key "101" from a row
// ----------------------------------------------------------------
std::string make_key(const Row<double>& y_row) {
  std::string key = "";
  for(unsigned int i=0; i<y_row.n_elem; ++i) {
    key += (std::isfinite(y_row[i])) ? "1" : "0";
  }
  return key;
}

// ----------------------------------------------------------------
// 1. FAST SUMMARY GENERATOR
//    Replaces R-side model.matrix logic. 
//    Builds X'X, X'Z, Z'Z, etc., without ever creating the full Z matrix.
// ----------------------------------------------------------------
// [[Rcpp::export]]
List get_summaries_fast_cpp(const mat& X, const mat& Y, 
                            const CharacterVector& sites, 
                            const CharacterVector& patients) {
  
  int n = X.n_rows;
  int px = X.n_cols;
  int R = Y.n_cols;
  
  // 1. Hierarchical Map: Site -> Row Indices
  std::map<std::string, std::vector<int> > site_map;
  for(int i=0; i<n; ++i) {
    std::string s = as<std::string>(sites[i]);
    site_map[s].push_back(i);
  }
  
  List Sh; 
  
  // 2. Iterate over Sites
  for(auto const& [site_name, row_indices] : site_map) {
    
    // Map Patient ID -> Local Integer ID (0 to N_pats-1)
    std::map<std::string, int> pat_to_local_id;
    int p_counter = 0;
    
    std::vector<int> valid_rows; 
    std::vector<int> valid_local_pids;
    std::vector<std::string> row_keys;
    
    // Filter valid rows and build local patient mapping
    for(int idx : row_indices) {
      std::string key = make_key(Y.row(idx));
      
      // Check if row is completely empty (all NA)
      bool is_empty = true;
      for(char c : key) if(c=='1') is_empty = false;
      if(is_empty) continue; 
      
      std::string pid = as<std::string>(patients[idx]);
      if(pat_to_local_id.find(pid) == pat_to_local_id.end()) {
        pat_to_local_id[pid] = p_counter++;
      }
      
      valid_rows.push_back(idx);
      valid_local_pids.push_back(pat_to_local_id[pid]);
      row_keys.push_back(key);
    }
    
    if(valid_rows.empty()) continue;
    
    // Dimensions for Z: [Intercept_Site, Pat_1 ... Pat_k]
    int n_pats = p_counter; 
    int pz = 1 + n_pats; 
    
    List site_list;
    
    // Group by Pattern within Site
    std::map<std::string, std::vector<int>> pattern_map;
    for(size_t i=0; i<row_keys.size(); ++i) {
      pattern_map[row_keys[i]].push_back(i); // Store index relative to 'valid_rows'
    }
    
    for(auto const& [key, rel_indices] : pattern_map) {
      
      int Nh = rel_indices.size();
      
      mat Xs(Nh, px);
      mat Ys_full(Nh, R); 
      uvec pids(Nh); 
      
      // Identify observed columns
      std::vector<uword> temp_idx;
      for(size_t c=0; c<key.length(); ++c) {
        if(key[c] == '1') temp_idx.push_back(c);
      }
      uvec idx_outcomes = conv_to<uvec>::from(temp_idx);
      int s = idx_outcomes.n_elem;
      
      // Extract Data
      for(int i=0; i<Nh; ++i) {
        int original_idx = valid_rows[rel_indices[i]];
        Xs.row(i) = X.row(original_idx);
        Ys_full.row(i) = Y.row(original_idx);
        pids[i] = valid_local_pids[rel_indices[i]];
      }
      mat Ys = Ys_full.cols(idx_outcomes);
      
      // --- Compute Cross Products ---
      
      // 1. ShX = X'X
      mat ShX = Xs.t() * Xs;
      
      // 2. ShXZ = X'Z
      // Z col 0 is site intercept (ones). Z col k+1 is patient k indicator.
      mat ShXZ = zeros(px, pz);
      ShXZ.col(0) = sum(Xs, 0).t(); // Sum columns for Site Intercept
      
      for(int i=0; i<Nh; ++i) {
        int pid = pids[i];
        ShXZ.col(pid + 1) += Xs.row(i).t();
      }
      
      // 3. ShXY = X'Y
      mat ShXY = Xs.t() * Ys;
      
      // 4. ShZ = Z'Z (Count Matrix)
      mat ShZ = zeros(pz, pz);
      ShZ(0,0) = Nh; // Site-Site
      
      for(int i=0; i<Nh; ++i) {
        int pid = pids[i];
        // Site-Patient overlap
        ShZ(0, pid+1) += 1;
        ShZ(pid+1, 0) += 1;
        // Patient-Patient diagonal
        ShZ(pid+1, pid+1) += 1;
      }
      
      // 5. ShZY = Z'Y
      mat ShZY = zeros(pz, s);
      ShZY.row(0) = sum(Ys, 0); // Site sum
      
      for(int i=0; i<Nh; ++i) {
        int pid = pids[i];
        ShZY.row(pid + 1) += Ys.row(i); // Patient sum
      }
      
      // 6. ShYY = Y'Y
      mat ShYY = Ys.t() * Ys;
      
      // Return 1-based index for R compatibility
      uvec idx_outcomes_r = idx_outcomes + 1; 
      
      List pat_list = List::create(
        Named("key") = key,
        Named("s") = s,
        Named("idx_outcomes") = idx_outcomes_r,
        Named("ShX") = ShX,
        Named("ShXZ") = ShXZ,
        Named("ShXY") = ShXY,
        Named("ShZ") = ShZ,
        Named("ShZY") = ShZY,
        Named("ShYY") = ShYY,
        Named("Nh") = Nh,
        Named("pzh") = pz
      );
      site_list[key] = pat_list;
    }
    
    Sh[site_name] = site_list;
  }
  
  return Sh;
}


// ----------------------------------------------------------------
// 2. PROFILE LIKELIHOOD CALCULATOR
//    Aggregates summaries to compute log-lik and terms for beta/sigma
// ----------------------------------------------------------------
// [[Rcpp::export]]
List compute_site_stats_cpp(List ShPat, 
                            vec par, 
                            int R, 
                            int px, 
                            std::string corstr,
                            double rho_fixed,
                            mat Corr_full,
                            bool estimate_rho) {
  
  // Global Accumulators
  double lpterm1 = 0;
  double lpterm2 = 0;
  mat bterm1 = zeros(R * px, R * px);
  vec bterm2 = zeros(R * px);
  
  double Nsum_rows = 0;
  double Nobs_total = 0;
  
  CharacterVector site_names = ShPat.names();
  int K = site_names.size();
  
  double theta_u = par[0];
  
  for (int h = 0; h < K; ++h) {
    std::string site = as<std::string>(site_names[h]);
    List S_h = ShPat[site];
    if (S_h.size() == 0) continue;
    
    double theta_vh = par[1 + h];
    
    // --- Priors ---
    List first_pat = S_h[0];
    int pzh = first_pat["pzh"];
    
    vec v0_diag = zeros(pzh);
    v0_diag[0] = theta_u;
    for(int i=1; i<pzh; ++i) v0_diag[i] = theta_vh;
    
    double logdet_Theta_h = R * sum(log(v0_diag));
    mat V0_inv = diagmat(1.0 / v0_diag);
    mat Theta_h_inv = kron(eye(R, R), V0_inv);
    
    // --- Pattern Accumulation ---
    mat SxxR_sum = zeros(R * px, R * px);
    mat Sxz_sum  = zeros(R * px, R * pzh);
    mat Szz_sum  = zeros(R * pzh, R * pzh);
    vec sxy_sum  = zeros(R * px);
    vec szy_sum  = zeros(R * pzh);
    double ytildeY_sum = 0;
    double logdet_Rcorr_sum = 0;
    double Nh_site = 0;
    
    for (int k = 0; k < S_h.size(); ++k) {
      List S = S_h[k];
      
      uvec o = as<uvec>(S["idx_outcomes"]) - 1; 
      int s = S["s"];
      int Nh = S["Nh"];
      Nh_site += Nh;
      Nobs_total += Nh * s;
      
      // Get Summary Matrices
      mat ShX = as<mat>(S["ShX"]);
      mat ShXZ = as<mat>(S["ShXZ"]);
      mat ShXY = as<mat>(S["ShXY"]);
      mat ShZ = as<mat>(S["ShZ"]);
      mat ShZY = as<mat>(S["ShZY"]);
      mat ShYY = as<mat>(S["ShYY"]);
      
      // Handle Correlation (Subset Inversion)
      mat Rinve_s;
      double ld_s = 0;
      
      if (corstr == "independence") {
        Rinve_s = eye(s, s);
        ld_s = 0;
      } else if (corstr == "exchangeable") {
        double rho_curr = (estimate_rho) ? par[par.size() - 1] : rho_fixed;
        double a = 1.0 / (1.0 - rho_curr);
        double b = rho_curr / (1.0 - rho_curr + rho_curr * s);
        Rinve_s = a * (eye(s, s) - b * ones(s, s));
        ld_s = (s - 1) * log(1.0 - rho_curr) + log(1.0 - rho_curr + rho_curr * s);
      } else {
        mat Coo = Corr_full.submat(o, o);
        // Robust Cholesky
        mat Uoo = chol(Coo + eye(s,s)*1e-9);
        Rinve_s = inv(trimatu(Uoo)); 
        Rinve_s = Rinve_s * Rinve_s.t();
        ld_s = 2.0 * sum(log(Uoo.diag()));
      }
      
      logdet_Rcorr_sum += Nh * ld_s;
      
      // Embed
      mat Rinve_embed = zeros(R, R);
      Rinve_embed.submat(o, o) = Rinve_s;
      
      SxxR_sum += kron(Rinve_embed, ShX);
      Sxz_sum  += kron(Rinve_embed, ShXZ);
      Szz_sum  += kron(Rinve_embed, ShZ);
      
      mat ShXY_full = zeros(ShXY.n_rows, R);
      mat ShZY_full = zeros(ShZY.n_rows, R);
      ShXY_full.cols(o) = ShXY * Rinve_s;
      ShZY_full.cols(o) = ShZY * Rinve_s;
      
      sxy_sum += vectorise(ShXY_full);
      szy_sum += vectorise(ShZY_full);
      
      ytildeY_sum += trace(Rinve_s * ShYY);
    }
    
    // --- Woodbury Step ---
    mat A_h = Theta_h_inv + Szz_sum;
    mat A_h_inv;
    double logdet_Ah;
    double sign;
    
    // Fast Symmetric Positive Definite Inverse
    bool success = log_det(logdet_Ah, sign, A_h);
    if(success) {
      A_h_inv = inv_sympd(A_h);
    } else {
      A_h_inv = inv(A_h + eye(size(A_h))*1e-6);
      log_det(logdet_Ah, sign, A_h);
    }
    
    lpterm1 += logdet_Rcorr_sum + logdet_Theta_h + logdet_Ah;
    
    // Projection
    mat term_xz = Sxz_sum * A_h_inv;
    bterm1 += SxxR_sum - term_xz * Sxz_sum.t();
    bterm2 += sxy_sum - term_xz * szy_sum;
    lpterm2 += ytildeY_sum - as_scalar(szy_sum.t() * A_h_inv * szy_sum);
    
    Nsum_rows += Nh_site;
  }
  
  return List::create(
    Named("lpterm1") = lpterm1,
    Named("lpterm2") = lpterm2,
    Named("bterm1") = bterm1,
    Named("bterm2") = bterm2,
    Named("Nobs_total") = Nobs_total,
    Named("Nsum_rows") = Nsum_rows
  );
}


// ----------------------------------------------------------------
// 3. ANALYTIC GRADIENT (Corrected Woodbury Formulation)
// ----------------------------------------------------------------
// [[Rcpp::export]]
vec compute_gradient_theta_cpp(List ShPat, 
                               vec par, 
                               int R, 
                               int px, 
                               std::string corstr,
                               double rho_fixed,
                               mat Corr_full,
                               bool estimate_rho,
                               vec beta, 
                               double sigma2) { 
  
  CharacterVector site_names = ShPat.names();
  int K = site_names.size();
  
  // Gradient vector for Thetas [theta_u, theta_v1 ... theta_vK]
  // We compute d(LogLik)/d(Theta). 
  // R will negate this to minimize Negative LogLik.
  vec grad = zeros(1 + K);
  double theta_u = par[0];
  
  for (int h = 0; h < K; ++h) {
    std::string site = as<std::string>(site_names[h]);
    List S_h = ShPat[site];
    if (S_h.size() == 0) continue;
    
    double theta_vh = par[1 + h];
    
    // --- Reconstruct Matrices (Same logic as profile) ---
    List first_pat = S_h[0];
    int pzh = first_pat["pzh"];
    
    vec v0_diag = zeros(pzh);
    v0_diag[0] = theta_u;
    for(int i=1; i<pzh; ++i) v0_diag[i] = theta_vh;
    
    mat V0_inv = diagmat(1.0 / v0_diag);
    mat Theta_h_inv = kron(eye(R, R), V0_inv);
    
    mat Szz_sum  = zeros(R * pzh, R * pzh);
    vec szy_sum  = zeros(R * pzh);
    mat Sxz_sum  = zeros(R * px, R * pzh); 
    
    for (int k = 0; k < S_h.size(); ++k) {
      List S = S_h[k];
      uvec o = as<uvec>(S["idx_outcomes"]) - 1;
      int s = S["s"];
      
      mat Rinve_s;
      if (corstr == "exchangeable") {
        double rho_curr = (estimate_rho) ? par[par.size()-1] : rho_fixed;
        double a = 1.0/(1.0-rho_curr);
        double b = rho_curr/(1.0-rho_curr+rho_curr*s);
        Rinve_s = a * (eye(s,s) - b*ones(s,s));
      } else if (corstr == "independence") {
        Rinve_s = eye(s, s);
      } else {
        mat Coo = Corr_full.submat(o, o);
        Rinve_s = inv(symmatu(Coo + eye(s,s)*1e-9));
      }
      
      mat Rinve_embed = zeros(R, R);
      Rinve_embed.submat(o, o) = Rinve_s;
      
      mat ShZ = as<mat>(S["ShZ"]);
      mat ShZY = as<mat>(S["ShZY"]);
      mat ShXZ = as<mat>(S["ShXZ"]);
      
      Szz_sum += kron(Rinve_embed, ShZ);
      Sxz_sum += kron(Rinve_embed, ShXZ);
      
      mat ShZY_full = zeros(ShZY.n_rows, R);
      ShZY_full.cols(o) = ShZY * Rinve_s;
      szy_sum += vectorise(ShZY_full);
    }
    
    mat A_h = Theta_h_inv + Szz_sum;
    mat A_h_inv;
    // Robust Inverse
    try {
      A_h_inv = inv_sympd(A_h);
    } catch(...) {
      A_h_inv = inv(A_h + eye(size(A_h))*1e-8);
    }
    
    // --- Gradient Calculation (Woodbury Explicit) ---
    
    // 1. Calculate intermediate vector w = A^-1 * Z' R^-1 (y - Xb)
    // v_ZR = szy - Sxz' beta
    vec v_ZR = szy_sum - Sxz_sum.t() * beta;
    vec w = A_h_inv * v_ZR;
    
    // 2. Precompute Diagonal of A_inv for the Trace term
    vec diag_A_inv = A_h_inv.diag();
    
    // --- Accumulate Site Intercept (theta_u) ---
    // Indices: 0, pzh, 2*pzh, ... (R times)
    double tr_term_u = 0;
    double data_term_u = 0;
    int count_u = 0;
    
    for(int r=0; r<R; ++r) {
      int idx = r * pzh; // Index for site intercept
      
      // Trace part: (1/theta - 1/theta^2 * A_ii)
      tr_term_u += (1.0/theta_u) - (1.0/(theta_u*theta_u)) * diag_A_inv[idx];
      
      // Data part: (w_i^2)
      data_term_u += w[idx] * w[idx];
      count_u++;
    }
    // Apply sigma2 scaling ONLY to data term
    data_term_u = (1.0 / (sigma2 * theta_u * theta_u)) * data_term_u;
    
    // dLL/dTheta = 0.5 * (Data - Trace)
    grad[0] += 0.5 * (data_term_u - tr_term_u);
    
    
    // --- Accumulate Patient Intercepts (theta_vh) ---
    // Indices: All except the site ones
    double tr_term_v = 0;
    double data_term_v = 0;
    
    for(int r=0; r<R; ++r) {
      for(int p=1; p<pzh; ++p) { 
        int idx = r * pzh + p;
        
        tr_term_v += (1.0/theta_vh) - (1.0/(theta_vh*theta_vh)) * diag_A_inv[idx];
        data_term_v += w[idx] * w[idx];
      }
    }
    
    data_term_v = (1.0 / (sigma2 * theta_vh * theta_vh)) * data_term_v;
    grad[1 + h] = 0.5 * (data_term_v - tr_term_v);
  }
  
  return grad;
}