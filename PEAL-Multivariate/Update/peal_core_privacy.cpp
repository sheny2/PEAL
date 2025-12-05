// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <map>
#include <vector>

using namespace Rcpp;
using namespace arma;

std::string make_key_v3(const Row<double>& y_row) {
  std::string key = "";
  for(unsigned int i=0; i<y_row.n_elem; ++i) key += (std::isfinite(y_row[i])) ? "1" : "0";
  return key;
}

// ----------------------------------------------------------------
// PRIVACY-PRESERVING SITE UPDATE
// Inputs:
//   - Old_Sh_List: The anonymous summary list for this site (from Server)
//   - Old_Map_List: The PRIVATE mapping (PatientID -> Index) (from Local Disk)
//   - X_new, Y_new: New raw data
// ----------------------------------------------------------------
// [[Rcpp::export]]
List update_site_privacy_cpp(List Old_Sh_List, List Old_Map_List, 
                             mat X_new, mat Y_new, 
                             CharacterVector patients_new) {
  
  // 1. Load the Private Local Map
  // Map: Key = PatientID, Value = Matrix Index (0-based relative to Z columns)
  std::map<std::string, int> local_map;
  CharacterVector old_pats = Old_Map_List["ids"];
  IntegerVector old_idxs = Old_Map_List["indices"];
  
  int max_idx = 0;
  if(old_pats.size() > 0) {
    for(int i=0; i<old_pats.size(); ++i) {
      local_map[as<std::string>(old_pats[i])] = old_idxs[i];
      if(old_idxs[i] > max_idx) max_idx = old_idxs[i];
    }
  }
  
  // 2. Scan New Data & Update Map (Expand dimensions if needed)
  std::vector<int> row_indices_for_calc;
  std::vector<int> pid_indices_for_calc; // The target column for each row
  std::vector<std::string> row_keys;
  
  int p_counter = (local_map.empty()) ? 1 : max_idx + 1; // Start after last known
  
  for(int i=0; i<X_new.n_rows; ++i) {
    std::string key = make_key_v3(Y_new.row(i));
    bool is_empty = true;
    for(char c : key) if(c=='1') is_empty = false;
    if(is_empty) continue;
    
    std::string pid = as<std::string>(patients_new[i]);
    
    // If Patient is NEW, assign new index
    if(local_map.find(pid) == local_map.end()) {
      local_map[pid] = p_counter++;
    }
    
    row_indices_for_calc.push_back(i);
    pid_indices_for_calc.push_back(local_map[pid]);
    row_keys.push_back(key);
  }
  
  int pz_new = p_counter; // New dimension size (1 site + patients)
  // If no previous data, pz starts at 1 (site intercept)
  if(local_map.empty() && pz_new == 1) pz_new = 1; 
  
  // 3. Update Matrices
  List New_Sh_List = clone(Old_Sh_List);
  
  // Group new rows by pattern
  std::map<std::string, std::vector<int>> pattern_map;
  for(size_t i=0; i<row_keys.size(); ++i) pattern_map[row_keys[i]].push_back(i);
  
  // Iterate over patterns found in NEW data
  for(auto const& [key, rel_indices] : pattern_map) {
    
    // Check if this pattern existed in Old List
    List Pat_Obj;
    bool exists = false;
    
    // We have to search the list by name (or key property)
    // R lists with names are slower to search in C++, assuming straightforward lookup
    if(New_Sh_List.containsElementNamed(key.c_str())) {
      Pat_Obj = New_Sh_List[key];
      exists = true;
    }
    
    // Prepare New Data block
    int Nh_new = rel_indices.size();
    int s = 0; 
    
    // Determine s and indices
    std::vector<uword> temp_idx;
    for(size_t c=0; c<key.length(); ++c) if(key[c] == '1') temp_idx.push_back(c);
    uvec idx_outcomes = conv_to<uvec>::from(temp_idx);
    s = idx_outcomes.n_elem;
    
    mat Xs(Nh_new, X_new.n_cols);
    mat Ys_col(Nh_new, s);
    
    // Aggregators for this batch
    mat Batch_ShX = zeros(X_new.n_cols, X_new.n_cols);
    mat Batch_ShXZ = zeros(X_new.n_cols, pz_new);
    mat Batch_ShZ = zeros(pz_new, pz_new);
    mat Batch_ShXY = zeros(X_new.n_cols, s);
    mat Batch_ShZY = zeros(pz_new, s);
    mat Batch_ShYY = zeros(s, s);
    
    for(int i=0; i<Nh_new; ++i) {
      int r_idx = row_indices_for_calc[rel_indices[i]];
      int p_idx = pid_indices_for_calc[rel_indices[i]]; // The Matrix Column
      
      rowvec x_r = X_new.row(r_idx);
      rowvec y_r_full = Y_new.row(r_idx);
      rowvec y_r = y_r_full.elem(idx_outcomes).t();
      
      Xs.row(i) = x_r;
      Ys_col.row(i) = y_r;
      
      // Incremental Accumulation
      Batch_ShX += x_r.t() * x_r;
      Batch_ShXY += x_r.t() * y_r;
      Batch_ShYY += y_r.t() * y_r;
      
      // Sparse Z Logic
      // Z row has 1 at col 0 (Site) and 1 at col p_idx (Patient)
      Batch_ShXZ.col(0) += x_r.t();
      Batch_ShXZ.col(p_idx) += x_r.t();
      
      Batch_ShZ(0,0) += 1;
      Batch_ShZ(0, p_idx) += 1;
      Batch_ShZ(p_idx, 0) += 1;
      Batch_ShZ(p_idx, p_idx) += 1;
      
      Batch_ShZY.row(0) += y_r;
      Batch_ShZY.row(p_idx) += y_r;
    }
    
    // Combine with Old Data (if exists)
    if(exists) {
      // Need to resize Old Matrices to pz_new before adding
      int old_pz = Pat_Obj["pzh"];
      
      mat Old_ShX = as<mat>(Pat_Obj["ShX"]);
      mat Old_ShXY = as<mat>(Pat_Obj["ShXY"]);
      mat Old_ShYY = as<mat>(Pat_Obj["ShYY"]);
      
      // Resizing Z-related matrices
      mat Old_ShZ = as<mat>(Pat_Obj["ShZ"]);
      mat Old_ShZ_big = zeros(pz_new, pz_new);
      Old_ShZ_big.submat(0,0, old_pz-1, old_pz-1) = Old_ShZ;
      
      mat Old_ShXZ = as<mat>(Pat_Obj["ShXZ"]);
      mat Old_ShXZ_big = zeros(X_new.n_cols, pz_new);
      Old_ShXZ_big.cols(0, old_pz-1) = Old_ShXZ;
      
      mat Old_ShZY = as<mat>(Pat_Obj["ShZY"]);
      mat Old_ShZY_big = zeros(pz_new, s);
      Old_ShZY_big.rows(0, old_pz-1) = Old_ShZY;
      
      // ADD
      Pat_Obj["ShX"] = Old_ShX + Batch_ShX;
      Pat_Obj["ShXY"] = Old_ShXY + Batch_ShXY;
      Pat_Obj["ShYY"] = Old_ShYY + Batch_ShYY;
      Pat_Obj["ShZ"] = Old_ShZ_big + Batch_ShZ;
      Pat_Obj["ShXZ"] = Old_ShXZ_big + Batch_ShXZ;
      Pat_Obj["ShZY"] = Old_ShZY_big + Batch_ShZY;
      Pat_Obj["Nh"] = as<int>(Pat_Obj["Nh"]) + Nh_new;
      Pat_Obj["pzh"] = pz_new;
      
      New_Sh_List[key] = Pat_Obj;
      
    } else {
      // Create New Pattern Entry
      uvec idx_outcomes_r = idx_outcomes + 1;
      List New_Pat = List::create(
        Named("key") = key, Named("s") = s, Named("idx_outcomes") = idx_outcomes_r,
        Named("ShX") = Batch_ShX, Named("ShXZ") = Batch_ShXZ, Named("ShXY") = Batch_ShXY,
              Named("ShZ") = Batch_ShZ, Named("ShZY") = Batch_ShZY, Named("ShYY") = Batch_ShYY,
                    Named("Nh") = Nh_new, Named("pzh") = pz_new
      );
      New_Sh_List[key] = New_Pat;
    }
  }
  
  // 4. Update Remaining Patterns (Padding only)
  // We updated patterns present in New Data. 
  // Patterns NOT in New Data but in Old Data must still be resized to pz_new.
  
  CharacterVector all_keys = New_Sh_List.names();
  for(int k=0; k<all_keys.size(); ++k) {
    std::string key = as<std::string>(all_keys[k]);
    // If we just updated it, pzh is already pz_new. Check.
    List S = New_Sh_List[key];
    int current_pz = S["pzh"];
    
    if(current_pz < pz_new) {
      // This pattern wasn't observed in the new batch, but needs padding
      mat Old_ShZ = as<mat>(S["ShZ"]);
      mat Old_ShZ_big = zeros(pz_new, pz_new);
      Old_ShZ_big.submat(0,0, current_pz-1, current_pz-1) = Old_ShZ;
      S["ShZ"] = Old_ShZ_big;
      
      mat Old_ShXZ = as<mat>(S["ShXZ"]);
      mat Old_ShXZ_big = zeros(X_new.n_cols, pz_new);
      Old_ShXZ_big.cols(0, current_pz-1) = Old_ShXZ;
      S["ShXZ"] = Old_ShXZ_big;
      
      mat Old_ShZY = as<mat>(S["ShZY"]);
      mat Old_ShZY_big = zeros(pz_new, as<mat>(S["ShZY"]).n_cols);
      Old_ShZY_big.rows(0, current_pz-1) = Old_ShZY;
      S["ShZY"] = Old_ShZY_big;
      
      S["pzh"] = pz_new;
      New_Sh_List[key] = S;
    }
  }
  
  // 5. Serialize Updated Map for Return
  std::vector<std::string> new_ids_vec;
  std::vector<int> new_idxs_vec;
  for(auto const& [pid, idx] : local_map) {
    new_ids_vec.push_back(pid);
    new_idxs_vec.push_back(idx);
  }
  
  List Updated_Map = List::create(Named("ids") = new_ids_vec, Named("indices") = new_idxs_vec);
  
  return List::create(Named("summary") = New_Sh_List, Named("map") = Updated_Map);
}

// NOTE: You still need the 'compute_site_stats_cpp' and 'compute_gradient_theta_cpp' 
// from the previous iteration included in this file for the Optimization step.