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
//' Samples using SMC as in the framework laid out in \insertCite{doucet}{SC2GroupProject}.
//' This function uses the SMC framework to implement Example 1 of Section 5.3 in 
//' \insertCite{doucet}{SC2GroupProject} given in section 5.2.
//'
//' @param x_0 starting value
//' @param t The length of each time-step. 
//' @param n Number of particles to simulate
//' @param timesteps Number of time-steps to simulate (`t*timesteps` should give the end time.) 
//' @param varphi varphi function used in approximation
//' @param resample Logical indicating whether to resample or not
//'
//' @examples 
//' # Generate 10,000 samples and at 20 time-steps up to time 1
//' out <- diff_SMC_sin(x_0=5, t=1/20, n=10000, timesteps=20, varphi=0.9,TRUE)
//' @return A list with matrices containing various information from each observation, each entry gives an individual observation with 
//' each column being a separate time-point and reach row being a separate  particle
//' \item{Xmat - }{A matrix containing the simulated samples.}
//' \item{Wmat - }{A matrix containing the weights associated with each sample}
//' \item{gmat - }{A matrix containing the value of the function g evaulated at each sample.}
//'  
//' @references
//' \insertRef{doucet}{SC2GroupProject}
//'  
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