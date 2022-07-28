#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;



NumericVector loglike_exact_oneterm_C(NumericVector value) {
  int n = value.size();
  
  NumericVector term(n);
  for (int i = 0; i < n; i++){
    NumericVector value_j = clone(value);
    value_j[i] = 1 - value[i];
    
    
    double sum_ctl = 0;
    for (int j = 0; j < n; j++){
      if (j == i) {
        continue;
      }else {
        sum_ctl += value[j];
      }
      
    }
    term[i] = value_j[i] * sum_ctl + .00001;
  }
  
  return term;
}


// [[Rcpp::export]]
double loglike_exact_C(CharacterVector Groupnames, List idIndex, NumericVector score, double c){
  int n = Groupnames.size();
  
  
  double loglike = 0.0;
  for (int i=0; i < n; i++){
    //print(Groupnames[i]);
    String Groupname = Groupnames[i];
    NumericVector idx = idIndex[Groupname];
    
    //print(idx);

    int n_idx = idx.size();
    NumericVector value(n_idx);
    for (int j=0; j < n_idx; j++){
      value[j] = 1 - 1*(score[idx[j] - 1] > c);
    }
    //print(value);
    
    NumericVector term = loglike_exact_oneterm_C(value);
    // print(term);
            
    double term_sum = 0.0;
    for (int m = 0; m < term.size(); m++){
        term_sum += term[m];
    }
    
    double loglike_t = log(term[0]) - log(term_sum);
    
    loglike += loglike_t;
  } 
  
  return(-loglike);
}


// [[Rcpp::export]]
double loglike_kernel_C(CharacterVector Groupnames, List idIndex, NumericVector score, double c,
                        double h_n, double constant){
  int n = Groupnames.size();


  double loglike = 0.0;
  for (int i=0; i < n; i++){
    // print(Groupnames[i]);
    String Groupname = Groupnames[i];
    NumericVector idx = idIndex[Groupname];

    // print(idx);

    int n_idx = idx.size();
    NumericVector value(n_idx);
    for (int j=0; j < n_idx; j++){
      value[j] = 1 - (erf((score[idx[j]-1]-c)/(h_n*constant)/sqrt(2))*0.5 + 0.5);
    }
    // print(value);

    NumericVector term = loglike_exact_oneterm_C(value);
    // print(term);

    double term_sum = 0.0;
    for (int m = 0; m < term.size(); m++){
      term_sum += term[m];
    }


    double loglike_t = log(term[0]) - log(term_sum);

    loglike += loglike_t;
  }

  return(-loglike);
}
