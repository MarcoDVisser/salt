#include <Rcpp.h>

using namespace Rcpp;
//' tav function
//' 
//' Stern F. (1964), Transmission of isotropic radiation across an
//' interface between two dielectrics, Appl. Opt., 3(1):111-113.
//' Allen W.A. (1973), Transmission of isotropic light across a
//' dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
//' 63(6):664-666.
//' @references Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical
//' Properties Model Separating Photosynthetic Pigments, Remote Sensing of
//' Environment
//' 
//' @param(theta) Angle (in degrees!)
//' @param(ref) refractive index
//' @export
// [[Rcpp::export]]
NumericVector resHyp(NumericVector wl, NumericVector R, NumericVector mu, NumericVector sigma){

  NumericVector final(mu.size());
  NumericVector normdens(wl.size());
  NumericVector weight(wl.size());
  double normsum;
    
  for(int i =0; i < mu.size(); i++){
    normsum = 0;
    
    normdens = Rcpp::dnorm(wl,mu[i],sigma[i]);
    weight = R*normdens;
    normsum = sum(normdens);
      
      
      final[i] = sum(weight/normsum);
    
  }


return final;
  
}
