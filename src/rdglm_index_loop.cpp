#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void rdglm_index_loop(const IntegerVector & group, List & subset) {
  int n = group.size();
  int ngroups = subset.size();
  int g, j, k;
  IntegerVector subset_g;
  IntegerVector group_pos(ngroups); // must default to 0

  for(k = 0; k < n; k++) {
    g = group[k] - 1; //convert from R to C++ indexing
    j = group_pos[g];
    subset_g = subset(g); //extract the appropriate vector from the list
    subset_g[j] = k + 1; // output R indexing
    group_pos[g] = j+1; //where to put the next obs from this group
  }
}
