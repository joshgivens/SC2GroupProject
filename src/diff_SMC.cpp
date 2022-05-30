#include <cmath>
# include <Rcpp.h>
using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
NumericVector g(NumericVector x){
  double zeta=100.;
  double sqrt_zeta=sqrt(zeta);
  double omega2=1e-15;
  LogicalVector logicv=(x>(omega2-sqrt(zeta))) & (x<(sqrt(zeta)-omega2));
  NumericVector out_vector=ifelse(logicv,exp(1/(zeta-pow(x,2))),0.);
  return out_vector;
}

double sin_sampler(double x_0, double t){
  double pi = 3.141592653589793238463 ;
  double omega=pi;
  double mean=sin(x_0-omega)*t+x_0;
  double sample=R::rnorm(mean,t);
  return sample;
}

double sin_lik(double x, double x_0, double t){
  double pi = 3.141592653589793238463 ;
  double omega=pi;
  double mean=sin(x_0-omega)*t+x_0;
  double likelihood=R::dnorm4(x,mean,t,FALSE);
  return likelihood;
}



//' SMC Sampler 2
//' 
//' Samples using SMC as in the framework laid out in ...
//'
//' @param g Function which we are taking expectation over
//' @param p_delta true transition kernel between discrete time steps
//' @param M_T_samplers list of samplers for T approximating transition kernels
//' @param M_T_lik list of likelihood for T approximating kernels
//' @param n Number of particles to simulate
//' @param timesteps Number of time-steps to simulate
//' @param varphi varphi function used in approximation
//' @param x_0 starting value
//' @param resample Logical indicating whether to resample or not
//'
//' @return A list containing the simulated particles and their associated weights
//' @export
// [[Rcpp::export]]
List diff_SMC_sin(double x_0,double t,int n, int timesteps, double varphi, bool resample) {
  // initalise Matrix
  NumericMatrix Xmat(n,timesteps);
  NumericMatrix Wmat(n, timesteps);
  NumericMatrix gmat(n, timesteps);
  NumericMatrix gmat_transf(n, timesteps);
  NumericMatrix pmat(n, timesteps);
  NumericVector pvec0(n);
  NumericVector w_new(n);
  
  for (int i=0; i<n; i++){
    Xmat(i,0)=sin_sampler(x_0,t);
    pmat(i,0)=sin_lik(Xmat(i,0),x_0,t);
  }
  gmat(_,0)=g(Xmat(_,0));
  gmat_transf(_,0)=pow(abs(gmat(_,0)),varphi);
  Wmat(_,0)=gmat_transf(_,0);


  for (int time=1; time<timesteps; time++){
    NumericVector x_sampled(n);
    if (resample){
      NumericVector probs=Wmat(_,time-1)/sum(Wmat(_,time-1));
      NumericVector temp_samp=Xmat(_,time-1);
      x_sampled=Rcpp::sample(temp_samp,n,true,probs);
      NumericVector ones(n,1.);
      Wmat(_,time)=ones;
    }
    else{
      x_sampled=Xmat(_,time-1);
      Wmat(_,time)=Wmat(_,time-1);
    }
    for (int i=0;i<n;i++){
      Xmat(i,time)=sin_sampler(Xmat(i,time-1),t);
      pmat(i,time)=sin_lik(Xmat(i,time),Xmat(i,time-1),t); 
    }
    gmat(_,time)=g(Xmat(_,time));
    gmat_transf(_,time)=pow(abs(gmat(_,time)),varphi);
    Wmat(_,time)=gmat_transf(_,time)/gmat_transf(_,time-1);
  }
  
  List outlist = List::create(Named("Xmat")=Xmat, Named("Wmat")=Wmat,
                              Named("gmat")=gmat);
  return outlist;
}