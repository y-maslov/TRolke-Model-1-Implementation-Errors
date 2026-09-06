// @(#)root/physics:$Id$
// Author: Jan Conrad
 
/** \class TRolke128_num
    \legacy{TRolke128_num, Consider switching to RooStats.}
    \ingroup Physics
 This class computes confidence intervals for the rate of a Poisson
 process in the presence of uncertain background and/or efficiency.
 
 The treatment and the resulting limits are fully frequentist. The
 limit calculations make use of the profile likelihood method.
 
\author Jan Conrad (CERN) 2004, Updated: Johan Lundberg (CERN) 2009
 
 For a full list of methods and their syntax, and build instructions,
 consult the header file TRolke128_num.h.
 
 Examples/tutorials are found in the separate file Rolke.C
 
### TRolke128_num implements the following Models
 
 The signal is always assumed to be Poisson, with the following
 combinations of models of background and detection efficiency:
 
 If unsure, first consider model 3, 4 or 5.
 
1: SetPoissonBkgBinomEff(x,y,z,tau,m)
~~~
   Background: Poisson
   Efficiency: Binomial
~~~
   when the background is simultaneously measured
   from sidebands (or MC), and
   the signal efficiency was determined from Monte Carlo
 
2: SetPoissonBkgGaussEff(x,y,em,tau,sde)
~~~
   Background: Poisson
   Efficiency: Gaussian
~~~
   when the background is simultaneously measured
   from sidebands (or MC), and
   the efficiency is modeled as Gaussian
 
3: SetGaussBkgGaussEff(x,bm,em,sde,sdb)
~~~
   Background: Gaussian
   Efficiency: Gaussian
~~~
   when background and efficiency can both be
   modeled as Gaussian.
 
4: SetPoissonBkgKnownEff(x,y,tau,e)
~~~
   Background: Poisson
   Efficiency: Known
~~~
   when the background is simultaneously measured
   from sidebands (or MC).
 
5: SetGaussBkgKnownEff(x,bm,sdb,e)
~~~
   Background: Gaussian
   Efficiency: Known
~~~
   when background is Gaussian
 
6: SetKnownBkgBinomEff(x,z,b,m)
~~~
   Background: Known
   Efficiency: Binomial
~~~
   when signal efficiency was determined from Monte Carlo
 
7: SetKnownBkgGaussEff(x,em,sde,b)
~~~
   Background: Known
   Efficiency: Gaussian
~~~
   when background is known and efficiency Gaussian
 
### Parameters and further explanation
 
#### For all models:
~~~
   x = number of observed events in the experiment
~~~
   Efficiency (e or em) is the detection probability for signal.
   A low efficiency hence generally means weaker limits.
   If the efficiency of an experiment (with analysis cuts) is
   dealt with elsewhere, em or e can be set to one.
 
#### For Poisson background measurements (sideband or MC):
~~~
   y = number of observed events in background region
   tau =
       Either: the ratio between signal and background region
       in case background is observed.
       Or: the ratio between observed and simulated live-time
       in case background is determined from MC.
~~~
#### For Gaussian efficiency or background:
~~~
   bm  = estimate of the background
   sdb = corresponding standard deviation
 
   em  = estimate of the efficiency
   sde = corresponding standard deviation
~~~
   If the efficiency scale of dealt with elsewhere,
   set em to 1 and sde to the relative uncertainty.
 
#### For Binomial signal efficiency:
~~~
    m = number of MC events generated
    z = number of MC events observed
~~~
#### For the case of known background expectation or known efficiency:
~~~
    e = true efficiency (considered known)
    b = background expectation value (considered known)
~~~
 
 The confidence level (CL) is set either at construction
 time or with either of SetCL or SetCLSigmas
 
 The TRolke128_num method is very similar to the one used in MINUIT (MINOS).
 
 Two options are offered to deal with cases where the maximum likelihood
 estimate (MLE) is not in the physical region. Version "bounded likelihood"
 is the one used by MINOS if bounds for the physical region are chosen.
 Unbounded likelihood (the default) allows the MLE to be in the
 unphysical region. It has however better coverage.
 For more details consult the reference (see below).
 
 For a description of the method and its properties:
 
 W.Rolke, A. Lopez, J. Conrad and Fred James
 "Limits and Confidence Intervals in presence of nuisance parameters"
  http://lanl.arxiv.org/abs/physics/0403059
  Nucl.Instrum.Meth.A551:493-503,2005
 
#### Should I use TRolke128_num, TFeldmanCousins, TLimit?
 
   1. Does TRolke128_num make TFeldmanCousins obsolete?
      Certainly not. TFeldmanCousins is the fully frequentist construction and
      should be used in case of no (or negligible) uncertainties. It is however
      not capable of treating uncertainties in nuisance parameters. In other
      words, it does not handle background expectations or signal efficiencies
      which are known only with some limited accuracy.
      TRolke128_num is designed for this case and it is shown in the reference above
      that it has good coverage properties for most cases, and can be used
      where FeldmannCousins can't.
 
   2. What are the advantages of TRolke128_num over TLimit?
      TRolke128_num is fully frequentist. TLimit treats nuisance parameters Bayesian.
      For a coverage study of a Bayesian method refer to
      physics/0408039 (Tegenfeldt & J.C). However, this note studies
      the coverage of Feldman&Cousins with Bayesian treatment of nuisance
      parameters. To make a long story short: using the Bayesian method you
      might introduce a small amount of over-coverage (though I haven't shown it
      for TLimit). On the other hand, coverage of course is a not so interesting
      when you consider yourself a Bayesian.
*/
 
#include "TRolke128_num.h"
#include "TMath.h"
#include <quadmath.h>
#include <iostream>
 
// ClassImp(TRolke128_num);
 
////////////////////////////////////////////////////////////////////////////////
/// Constructor with optional Confidence Level argument.
/// 'option' is not used.
 
TRolke128_num::TRolke128_num(__float128 CL, Option_t * /*option*/)
:  fCL(CL),
   fUpperLimit(0.0),
   fLowerLimit(0.0),
   fBounding(false),  // true gives bounded likelihood
   fNumWarningsDeprecated1(0),
   fNumWarningsDeprecated2(0)
{
   SetModelParameters();
}
 
////////////////////////////////////////////////////////////////////////////////
/// Destructor.
 
TRolke128_num::~TRolke128_num()
{
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 1: Background - Poisson, Efficiency - Binomial
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - z   : number of MC events observed
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - m   : number of MC events generated
 
void TRolke128_num::SetPoissonBkgBinomEff(Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         y  ,       //   Int_t y,
         z  ,       //   Int_t z,
         0  ,       //   __float128 bm,
         0  ,       //   __float128 em,
         0  ,       //   __float128 e,
         1  ,       //   Int_t mid,
         0  ,       //   __float128 sde,
         0  ,       //   __float128 sdb,
         tau,       //   __float128 tau,
         0  ,       //   __float128 b,
         m);        //   Int_t m
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 2: Background - Poisson, Efficiency - Gaussian
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - em  : estimate of the efficiency
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - sde : efficiency estimate's standard deviation
 
void TRolke128_num::SetPoissonBkgGaussEff(Int_t x, Int_t y, __float128 em, __float128 tau, __float128 sde)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         y  ,       //   Int_t y,
         0  ,       //   Int_t z,
         0  ,       //   __float128 bm,
         em ,       //   __float128 em,
         0  ,       //   __float128 e,
         2  ,       //   Int_t mid,
         sde,       //   __float128 sde,
         0  ,       //   __float128 sdb,
         tau,       //   __float128 tau,
         0  ,       //   __float128 b,
         0);        //   Int_t m
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 3: Background - Gaussian, Efficiency - Gaussian (x,bm,em,sde,sdb)
///   - x   : number of observed events in the experiment
///   - bm  : estimate of the background
///   - em  : estimate of the efficiency
///   - sde : efficiency estimate's standard deviation
///   - sdb : background estimate's standard deviation
 
void TRolke128_num::SetGaussBkgGaussEff(Int_t x, __float128 bm, __float128 em, __float128 sde, __float128 sdb)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         0  ,       //   Int_t y,
         0  ,       //   Int_t z,
         bm ,       //   __float128 bm,
         em ,       //   __float128 em,
         0  ,       //   __float128 e,
         3  ,       //   Int_t mid,
         sde,       //   __float128 sde,
         sdb,       //   __float128 sdb,
         0  ,       //   __float128 tau,
         0  ,       //   __float128 b,
         0);        //   Int_t m
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 4: Background - Poisson, Efficiency - known     (x,y,tau,e)
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - e   : true efficiency (considered known)
 
void TRolke128_num::SetPoissonBkgKnownEff(Int_t x, Int_t y, __float128 tau, __float128 e)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         y  ,       //   Int_t y,
         0  ,       //   Int_t z,
         0  ,       //   __float128 bm,
         0  ,       //   __float128 em,
         e  ,       //   __float128 e,
         4  ,       //   Int_t mid,
         0  ,       //   __float128 sde,
         0  ,       //   __float128 sdb,
         tau,       //   __float128 tau,
         0  ,       //   __float128 b,
         0);        //   Int_t m
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 5: Background - Gaussian, Efficiency - known    (x,bm,sdb,e
///   - x   : number of observed events in the experiment
///   - bm  : estimate of the background
///   - sdb : background estimate's standard deviation
///   - e   : true efficiency (considered known)
 
void TRolke128_num::SetGaussBkgKnownEff(Int_t x, __float128 bm, __float128 sdb, __float128 e)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         0  ,       //   Int_t y,
         0  ,       //   Int_t z,
         bm ,       //   __float128 bm,
         0  ,       //   __float128 em,
         e  ,       //   __float128 e,
         5  ,       //   Int_t mid,
         0  ,       //   __float128 sde,
         sdb,       //   __float128 sdb,
         0  ,       //   __float128 tau,
         0  ,       //   __float128 b,
         0);        //   Int_t m
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 6: Background - known, Efficiency - Binomial    (x,z,m,b)
///   - x   : number of observed events in the experiment
///   - z   : number of MC events observed
///   - m   : number of MC events generated
///   - b   : background expectation value (considered known)
 
void TRolke128_num::SetKnownBkgBinomEff(Int_t x, Int_t z, Int_t m, __float128 b)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         0  ,       //   Int_t y
         z  ,       //   Int_t z,
         0  ,       //   __float128 bm,
         0  ,       //   __float128 em,
         0  ,       //   __float128 e,
         6  ,       //   Int_t mid,
         0  ,       //   __float128 sde,
         0  ,       //   __float128 sdb,
         0  ,       //   __float128 tau,
         b  ,       //   __float128 b,
         m);        //   Int_t m
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Model 7: Background - known, Efficiency - Gaussian    (x,em,sde,b)
///   - x   : number of observed events in the experiment
///   - em  : estimate of the efficiency
///   - sde : efficiency estimate's standard deviation
///   - b   : background expectation value (considered known)
 
void TRolke128_num::SetKnownBkgGaussEff(Int_t x, __float128 em, __float128 sde, __float128 b)
{
   SetModelParameters(
         x  ,       //   Int_t x,
         0  ,       //   Int_t y
         0  ,       //   Int_t z,
         0  ,       //   __float128 bm,
         em ,       //   __float128 em,
         0  ,       //   __float128 e,
         7  ,       //   Int_t mid,
         sde,       //   __float128 sde,
         0  ,       //   __float128 sdb,
         0  ,       //   __float128 tau,
         b  ,       //   __float128 b,
         0);        //   Int_t m
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculate and get the upper and lower limits for the pre-specified model.
 
bool TRolke128_num::GetLimits(__float128& low, __float128& high)
{
   if ((f_mid<1)||(f_mid>7)) {
      std::cerr << "TRolke128_num - Error: Model id "<< f_mid<<std::endl;
      if (f_mid<1) {
         std::cerr << "TRolke128_num - Please specify a model with e.g. 'SetGaussBkgGaussEff' (read the docs in Rolke.cxx )"<<std::endl;
      }
      return false;
   }
 
   ComputeInterval(f_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
   low = fLowerLimit;
   high = fUpperLimit;
   if (low < high) {
      return true;
   }else{
      std::cerr << "TRolke128_num - Warning: no limits found" <<std::endl;
      return false;
   }
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculate and get upper limit for the pre-specified model.
 
__float128 TRolke128_num::GetUpperLimit()
{
   __float128 low(0), high(0);
   GetLimits(low,high);
   return fUpperLimit;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculate and get lower limit for the pre-specified model.
 
__float128 TRolke128_num::GetLowerLimit()
{
   __float128 low(0), high(0);
   GetLimits(low,high);
   return fLowerLimit;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Return a simple background value (estimate/truth) given the pre-specified model.
 
__float128 TRolke128_num::GetBackground()
{
   __float128 background = 0;
   switch (f_mid) {
      case 1:
      case 2:
      case 4:
         background = f_y / f_tau;
         break;
      case 3:
      case 5:
         background = f_bm;
         break;
      case 6:
      case 7:
         background = f_b;
         break;
      default:
         std::cerr << "TRolke128_num::GetBackground(): Model NR: " <<
         f_mid << " unknown"<<std::endl;
         return 0;
   }
   return background;
}
 
////////////////////////////////////////////////////////////////////////////////
/// get the upper and lower average limits based on the specified model.
/// No uncertainties are considered for the Poisson weights in the averaging sum
 
bool TRolke128_num::GetSensitivity(__float128& low, __float128& high, __float128 pPrecision)
{
   __float128 background = GetBackground();
 
   __float128 weight = 0;
   __float128 weightSum = 0;
 
   int loop_x = 0;
 
   while (true) {
      ComputeInterval(loop_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
      weight = PoissonI128(loop_x, background);
      low += fLowerLimit * weight;
      high += fUpperLimit * weight;
      weightSum += weight;
      if (loop_x > (background + 1)) { // don't stop too early
         if ((weightSum > (1 - pPrecision)) || (weight < 1e-12)) break;
      }
      loop_x++;
   }
   low /= weightSum;
   high /= weightSum;
 
   return (low < high); // could also add more detailed test
}
 
////////////////////////////////////////////////////////////////////////////////
/// get the upper and lower limits for the outcome corresponding to
/// a given quantile.
/// For integral=0.5 this gives the median limits
/// in repeated experiments. The returned out_x is the corresponding
/// (e.g. median) value of x.
/// No uncertainties are considered for the Poisson weights when calculating
/// the Poisson integral.
 
bool TRolke128_num::GetLimitsQuantile(__float128& low, __float128& high, Int_t& out_x, __float128 integral)
{
   __float128 background = GetBackground();
   __float128 weight = 0;
   __float128 weightSum = 0;
   Int_t loop_x = 0;
 
   while (true) {
      weight = PoissonI128(loop_x, background);
      weightSum += weight;
      if (weightSum >= integral) {
         break;
      }
      loop_x++;
   }
 
   out_x = loop_x;
 
   ComputeInterval(loop_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
   low = fLowerLimit;
   high = fUpperLimit;
   return (low < high); // could also add more detailed test
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// get the upper and lower limits for the most likely outcome.
/// The returned out_x is the corresponding value of x
/// No uncertainties are considered for the Poisson weights when finding ML.
 
bool TRolke128_num::GetLimitsML(__float128& low, __float128& high, Int_t& out_x)
{
   __float128 background = GetBackground();
 
   Int_t loop_x = 0; // this can be optimized if needed.
   Int_t loop_max = 1000 + (Int_t)background; //     --||--
 
   __float128 max = PoissonI128(loop_x, background);
   while (loop_x <= loop_max) {
      if (PoissonI128(loop_x + 1, background) < max) {
         break;
      }
      loop_x++;
      max = PoissonI128(loop_x, background);
   }
   if (loop_x >= loop_max) {
      std::cout << "internal error finding maximum of distribution" << std::endl;
      return false;
   }
 
   out_x = loop_x;
 
   ComputeInterval(loop_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
   low = fLowerLimit;
   high = fUpperLimit;
   return (low < high); // could also add more detailed test
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// get the value of x corresponding to rejection of the null hypothesis.
/// This means a lower limit >0 with the pre-specified Confidence Level.
/// Optionally give maxtry; the maximum value of x to try. Of not, or if
/// maxtry<0 an automatic mode is used.
 
bool TRolke128_num::GetCriticalNumber(Int_t& ncrit, Int_t maxtry)
{
   __float128 background = GetBackground();
 
   int j = 0;
   int rolke_ncrit = -1;
   int maxj =maxtry ;
   if(maxtry<1){
     maxj = 1000 + (Int_t)background; // max value to try
   }
   for (j = 0;j < maxj;j++) {
      Int_t rolke_x = j;
      ComputeInterval(rolke_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
      __float128 rolke_ll = fLowerLimit;
      if (rolke_ll > 0) {
         rolke_ncrit = j;
         break;
      }
   }
 
   if (rolke_ncrit == -1) {
     std::cerr << "TRolke128_num GetCriticalNumber : Error: problem finding rolke inverse. Specify a larger maxtry value. maxtry was: " << maxj << ". highest x considered was j "<< j<< std::endl;
      ncrit = -1;
      return false;
   } else {
      ncrit = rolke_ncrit;
      return true;
   }
}
 
////////////////////////////////////////////////////////////////////////////////
/// Deprecated name for SetBounding.
 
void TRolke128_num::SetSwitch(bool bnd) {
   if(fNumWarningsDeprecated1<2){
      std::cerr << "*******************************************" <<std::endl;
      std::cerr << "TRolke128_num - Warning: 'SetSwitch' is deprecated and may be removed from future releases:" <<std::endl;
      std::cerr << " - Use 'SetBounding' instead "<<std::endl;
      std::cerr << "*******************************************" <<std::endl;
      fNumWarningsDeprecated1++;
   }
   SetBounding(bnd);
}
 
////////////////////////////////////////////////////////////////////////////////
/// Dump internals. Print members.
 
void TRolke128_num::Print(Option_t*) const {
   std::cout << "*******************************************" <<std::endl;
   std::cout << "* TRolke128_num::Print() - dump of internals:                " <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "* model id, mid = "<<f_mid <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "*             x = "<<f_x   <<std::endl;
   std::cout << "*            bm = "<< (double) f_bm  <<std::endl;
   std::cout << "*            em = "<<(double) f_em  <<std::endl;
   std::cout << "*           sde = "<< (double) f_sde <<std::endl;
   std::cout << "*           sdb = "<< (double) f_sdb <<std::endl;
   std::cout << "*             y = "<<f_y   <<std::endl;
   std::cout << "*           tau = "<< (double) f_tau <<std::endl;
   std::cout << "*             e = "<< (double) f_e   <<std::endl;
   std::cout << "*             b = "<< (double) f_b   <<std::endl;
   std::cout << "*             m = "<<f_m   <<std::endl;
   std::cout << "*             z = "<<f_z   <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "*            CL = "<< (double) fCL <<std::endl;
   std::cout << "*      Bounding = "<<fBounding <<std::endl;
   std::cout << "*"<<std::endl;
   std::cout << "* calculated on demand only:"<<std::endl;
   std::cout << "*   fUpperLimit = "<< (double) fUpperLimit<<std::endl;
   std::cout << "*   fLowerLimit = "<<(double) fLowerLimit<<std::endl;
   std::cout << "*******************************************" <<std::endl;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Deprecated and error prone model selection interface.
/// It's use is trongly discouraged. 'mid' is the model ID (1 to 7).
/// This method is provided for backwards compatibility/developer use only. */
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - z   : number of MC events observed
///   - bm  : estimate of the background
///   - em  : estimate of the efficiency
///   - e   : true efficiency (considered known)
///   - mid : internal model id (really, you should not use this method at all)
///   - sde : efficiency estimate's standard deviation
///   - sdb : background estimate's standard deviation
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - b   : background expectation value (considered known)
///   - m   : number of MC events generated
 
__float128 TRolke128_num::CalculateInterval(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m){
   if (fNumWarningsDeprecated2<2 ) {
      std::cerr << "*******************************************" <<std::endl;
      std::cerr << "TRolke128_num - Warning: 'CalculateInterval' is deprecated and may be removed from future releases:" <<std::endl;
      std::cerr << " - Use e.g. 'SetGaussBkgGaussEff' and 'GetLimits' instead (read the docs in Rolke.cxx )"<<std::endl;
      std::cerr << "*******************************************" <<std::endl;
      fNumWarningsDeprecated2++;
   }
   SetModelParameters(
         x,
         y,
         z,
         bm,
         em,
         e,
         mid,
         sde,
         sdb,
         tau,
         b,
         m);
   return ComputeInterval(f_x, f_y, f_z, f_bm, f_em, f_e, f_mid, f_sde, f_sdb, f_tau, f_b, f_m);
}
 
////////////////////////////////////////////////////////////////////////////////
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - z   : number of MC events observed
///   - bm  : estimate of the background
///   - em  : estimate of the efficiency
///   - e   : true efficiency (considered known)
///   - mid : internal model id
///   - sde : efficiency estimate's standard deviation
///   - sdb : background estimate's standard deviation
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - b   : background expectation value (considered known)
///   - m   : number of MC events generated
 
void TRolke128_num::SetModelParameters(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m)
{
   f_x   = x   ;
   f_y   = y   ;
   f_z   = z   ;
   f_bm  = bm  ;
   f_em  = em  ;
   f_e   = e   ;
   f_mid = mid ;
   f_sde = sde ;
   f_sdb = sdb ;
   f_tau = tau ;
   f_b   = b   ;
   f_m   = m   ;
}
 
void TRolke128_num::SetModelParameters()
{
/* Clear internal model */
   SetModelParameters(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
   f_mid=0;
}
 
////////////////////////////////////////////////////////////////////////////////
/// ComputeInterval, the internals.
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - z   : number of MC events observed
///   - bm  : estimate of the background
///   - em  : estimate of the efficiency
///   - e   : true efficiency (considered known)
///   - mid : internal model id (really, you should not use this method at all)
///   - sde : efficiency estimate's standard deviation
///   - sdb : background estimate's standard deviation
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - b   : background expectation value (considered known)
///   - m   : number of MC events generated
 
__float128 TRolke128_num::ComputeInterval(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m)
{
   //calculate interval
   Int_t done = 0;
   __float128 limit[2];
 
   limit[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
 
   if (limit[1] > 0) {
      done = 1;
   }
 
   if (! fBounding) {
 
      Int_t trial_x = x;
 
      while (done == 0) {
         trial_x++;
         limit[1] = Interval(trial_x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
         if (limit[1] > 0) done = 1;
      }
   }
 
   return limit[1];
}
 
////////////////////////////////////////////////////////////////////////////////
/// Internal helper function 'Interval'
///   - x   : number of observed events in the experiment
///   - y   : number of observed events in background region
///   - z   : number of MC events observed
///   - bm  : estimate of the background
///   - em  : estimate of the efficiency
///   - e   : true efficiency (considered known)
///   - mid : internal model id (really, you should not use this method at all)
///   - sde : efficiency estimate's standard deviation
///   - sdb : background estimate's standard deviation
///   - tau : ratio parameter (read TRolke128_num.cxx for details)
///   - b   : background expectation value (considered known)
///   - m   : number of MC events generated
 
__float128 TRolke128_num::Interval(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m)
{
   __float128 dchi2 = TMath::ChisquareQuantile(fCL, 1);
   __float128 tempxy[2], limits[2] = {0, 0};
   __float128 slope, fmid, low, flow, high, fhigh, test, ftest, mu0, maximum, target, l, f0;
   __float128 med = 0;
   __float128 maxiter = 1000, acc = 0.00001;
   // __float128 maxiter = 1e6, acc = 1e-7;
   Int_t i;
   Int_t bp = 0;
 
   if ((mid != 3) && (mid != 5)) bm = y;
   if ((mid == 3) || (mid == 5)) {
      if (bm == 0) bm = 0.00001;
   }
 
   if ((mid == 6) || (mid == 7)) {
      if (bm == 0) bm = 0.00001;
   }
 
   if ((mid <= 2) || (mid == 4)) bp = 1;
 
 
   if (bp == 1 && x == 0 && bm > 0) {
      for (i = 0; i < 2; i++) {
         x++;
         tempxy[i] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      }
 
      slope = tempxy[1] - tempxy[0];
      limits[1] = tempxy[0] - slope;
      limits[0] = 0.0;
      if (limits[1] < 0) limits[1] = 0.0;
      goto done;
   }
 
   if (bp != 1 && x == 0) {
 
      for (i = 0; i < 2; i++) {
         x++;
         tempxy[i] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      }
      slope = tempxy[1] - tempxy[0];
      limits[1] = tempxy[0] - slope;
      limits[0] = 0.0;
      if (limits[1] < 0) limits[1] = 0.0;
      goto done;
   }
 
   if (bp != 1  && bm == 0) {
      for (i = 0; i < 2; i++) {
         bm++;
         limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
         tempxy[i] = limits[1];
      }
      slope = tempxy[1] - tempxy[0];
      limits[1] = tempxy[0] - slope;
      if (limits[1] < 0) limits[1] = 0;
      goto done;
   }
 
   if (x == 0 && bm == 0) {
      x++;
      bm++;
      limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      tempxy[0] = limits[1];
      x  = 1;
      bm = 2;
      limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      tempxy[1] = limits[1];
      x  = 2;
      bm = 1;
      limits[1] = Interval(x, y, z, bm, em, e, mid, sde, sdb, tau, b, m);
      limits[1] = 3 * tempxy[0] - tempxy[1] - limits[1];
      if (limits[1] < 0) limits[1] = 0;
      goto done;
   }
 
   mu0 = Likelihood(0, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 1);
   maximum = Likelihood(0, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 2);
   test = 0;
   f0 = Likelihood(test, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
   if (fBounding) {
      if (mu0 < 0) maximum = f0;
   }
 
   target = maximum - dchi2;
   if (f0 > target) {
      limits[0] = 0;
   } else {
      if (mu0 < 0) {
         limits[0] = 0;
         limits[1] = 0;
      }
 
      low   = 0;
      flow  = f0;
      high  = mu0;
      fhigh = maximum;
      for (i = 0; i < maxiter; i++) {
         l = (target - fhigh) / (flow - fhigh);
         if (l < 0.2) l = 0.2;
         if (l > 0.8) l = 0.8;
         med = l * low + (1 - l) * high;
         if (med < 0.01) {
            limits[1] = 0.0;
            goto done;
         }
         fmid = Likelihood(med, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
         if (fmid > target) {
            high  = med;
            fhigh = fmid;
         } else {
            low  = med;
            flow = fmid;
         }
         if ((high - low) < acc*high) break;
      }
      limits[0] = med;
   }
 
   if (mu0 > 0) {
      low  = mu0;
      flow = maximum;
   } else {
      low  = 0;
      flow = f0;
   }
 
   test = low + 1 ;
   ftest = Likelihood(test, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
   if (ftest < target) {
      high  = test;
      fhigh = ftest;
   } else {
      slope = (ftest - flow) / (test - low);
      high  = test + (target - ftest) / slope;
      fhigh = Likelihood(high, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
   }
 
   for (i = 0; i < maxiter; i++) {
      l = (target - fhigh) / (flow - fhigh);
      if (l < 0.2) l = 0.2;
      if (l > 0.8) l = 0.8;
      med  = l * low + (1. - l) * high;
      fmid = Likelihood(med, x, y, z, bm, em, mid, sde, sdb, tau, b, m, 3);
 
      if (fmid < target) {
         high  = med;
         fhigh = fmid;
      } else {
         low  = med;
         flow = fmid;
      }
 
      if (high - low < acc*high) break;
   }
 
   limits[1] = med;
 
done:
 
   // apply known efficiency
   if ((mid == 4) || (mid == 5)) {
      limits[0] /= e;
      limits[1] /= e;
   }
 
   fUpperLimit = limits[1];
   fLowerLimit = Max128(limits[0], 0.0);
 
   return limits[1];
}
 
////////////////////////////////////////////////////////////////////////////////
/// Internal helper function
/// Chooses between the different profile likelihood functions to use for the
/// different models.
/// Returns evaluation of the profile likelihood functions.
 
__float128 TRolke128_num::Likelihood(__float128 mu, Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m, Int_t what)
{
   switch (mid) {
      case 1:
         return EvalLikeMod1(mu, x, y, z, tau, m, what);
      case 2:
         return EvalLikeMod2(mu, x, y, em, sde, tau, what);
      case 3:
         return EvalLikeMod3(mu, x, bm, em, sde, sdb, what);
      case 4:
         return EvalLikeMod4(mu, x, y, tau, what);
      case 5:
         return EvalLikeMod5(mu, x, bm, sdb, what);
      case 6:
         return EvalLikeMod6(mu, x, z, b, m, what);
      case 7:
         return EvalLikeMod7(mu, x, em, sde, b, what);
      default:
         std::cerr << "TRolke128_num::Likelihood(...): Model NR: " <<
         f_mid << " unknown"<<std::endl;
         return 0;
   }
 
   return 0;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 1:
/// Poisson background/ Binomial Efficiency
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod1(__float128 mu, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m, Int_t what)
{
   __float128 f  = 0.0Q;
   __float128 zm = __float128(z) / m;


   // std::cout << x << " " << y << " " << (double) tau << " " << (double) zm <<std::endl;

   // std::cout << (double)(x * (1 - (y / tau) / (__float128) x) / zm) <<std::endl; 
   // std::cout << (double)(x * (1 - (y / (__float128) x) /tau) / zm) <<std::endl;
   // std::cout << (double)((x - y /tau) / zm) <<std::endl;
 
   if (what == 1) {
      f = x * (1.0Q - (y / tau) / (__float128) x) / zm;
   }
 
   if (what == 2) {
      mu = x * (1.0Q - (y / tau) / (__float128) x) / zm;
      __float128 b  = y / tau;
      __float128 e = zm;
      f = LikeMod1(mu, b, e, x, y, z, tau, m);
   }
 
   if (what == 3) {
      if (mu == 0) {
         __float128 b = (x + y) / (1.0Q + tau);
         __float128 e = zm;
         f = LikeMod1(mu, b, e, x, y, z, tau, m);
      } else {
         __float128 e = 0.0Q;
         __float128 b = 0.0Q;
         

         ProfLikeMod1(mu, b, e, x, y, z, tau, m);
         f = LikeMod1(mu, b, e, x, y, z, tau, m);

         // std::cout << std::endl;
         // std::cout << "LikeMod1 old: " << (double) f << std::endl;
         // std::cout << " bOld: " << (double) b  << " eOld: " << (double) e << std::endl;

         ProfLikeMod1Fix(mu, b, e, x, y, z, tau, m);
         f = LikeMod1(mu, b, e, x, y, z, tau, m);
      }
   }
 
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 1:
/// Poisson background/ Binomial Efficiency
 
__float128 TRolke128_num::LikeMod1(__float128 mu, __float128 b, __float128 e, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m)
{
   __float128 s = e*mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
   __float128 bg = tau*b;
   __float128 llb =  -bg;
   if ( y > 0) llb =  y*Log128( bg) - bg - LogFactorial(y);
 
   __float128 lle = 0.0Q;  // binomial log-like
   if (z == 0)         lle = m * Log128(1-e);
   else if ( z == m)   lle = m * Log128(e);
   else                lle =   z * Log128(e) + (m - z)*Log128(1 - e) + LogFactorial(m) - LogFactorial(m-z) - LogFactorial(z);

   __float128 f = 2.0Q*( lls + llb + lle);

   // std::cout << "s: " << (double) s << " lls: " << (double) lls << " bg: " << (double) bg << " llb: " << (double) llb << " lle: " << (double) lle << " f: " << (double) f << std::endl;

   return f;
}
 
 
// this code is non-sense - // need to solve using Minuit
struct LikeFunction1 {
};
 
////////////////////////////////////////////////////////////////////////////////
/// Helper for calculation of estimates of efficiency and background for model 1 (search for extreme b and e)
 
void TRolke128_num::ProfLikeMod1(__float128 mu, __float128 &b, __float128 &e, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m)
{
   __float128 med = 0.0Q, fmid;
   Int_t maxiter = 1000;
   __float128 acc = 0.00001Q;
   __float128 emin = ((m + mu * tau) - Sqrt128((m + mu * tau) * (m + mu * tau) - 4.0Q * mu * tau * z)) / 2.0Q / mu / tau;
 
   __float128 low  = Max128(1e-10Q, emin + 1e-10Q);
   __float128 high = 1.0Q - 1e-10Q;
 
   for (Int_t i = 0; i < maxiter; i++) {
      med = (low + high) / 2.0Q;
 
      fmid = LikeGradMod1(med, mu, x, y, z, tau, m);
 
      if (high < 0.5Q) acc = 0.00001Q * high;
      else           acc = 0.00001Q * (1.0Q - high);
 
      if ((high - low) < acc*high) break;
 
      if (fmid > 0.0Q) low  = med;
      else         high = med;
   }
 
   e = med;
   __float128 eta = __float128(z) / e - __float128(m - z) / (1.0Q - e);
 
   b = __float128(y) / (tau - eta / mu); // YM: why is minus here?
}

void TRolke128_num::ProfLikeMod1Fix(__float128 mu, __float128 &b, __float128 &e, 
      Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m)
{
   __float128 A = - mu * tau * (1.0Q + tau)*(1.0Q + tau);
   __float128 B = mu * (1.0Q + tau) * (- mu * tau * (1.0Q + tau) + 2.0Q * tau * x + y + 
   3.0Q * tau * y + m * (1.0Q + tau));
   __float128 C = mu * (- (x + y) * (tau * x + 2.0Q * y + 3.0Q * tau * y) - 
   m * (1.0Q + tau) * (x + 2.0Q * y) + mu * (1.0Q + tau) * (y + z + tau * (x + 3.0Q * y + z)));
   __float128 D = mu * y * (((__float128) (x + y) * (x + y))  -
    mu * (x + 2.0Q * tau * x + 2.0Q * y + 3.0Q * tau * y + 2.0Q * (1.0Q + tau) * z) + m * (x + y));
   __float128 E = mu * mu * y * y * (x + y + z);

   auto poly4comp = [A, B, C, D, E] (__complex128 b) -> __complex128 {
      __complex128 b_2 = b * b;
      return A * b_2 * b_2 + B * b_2 * b + C * b_2 + D * b + E;
   };

   auto poly4rat = [A, B, C, D, E] (__float128 b) -> __float128 {
      __float128 b_2 = b * b;
      return A * b_2 * b_2 + B * b_2 * b + C * b_2 + D * b + E;
   };

   __float128 A_2 = A * A;
   __float128 B_2 = B * B;

   __float128 alpha = -3.0Q * B_2 / (8.0Q * A_2) + C / A;
   __float128 alpha_2 = alpha * alpha;
   __float128 beta = B_2 * B / (8.0Q * A_2 * A) - B * C / (2.0Q * A_2) + D / A;

   __float128 gamma = - 3.0Q * B_2 * B_2 / (256.0Q * A_2 * A_2) +
    B_2 * C / (16.0Q * A_2 * A) - B * D / (4.0Q * A_2) + E / A; 

   __float128 P = - alpha_2 / 12.0Q - gamma;
   __float128 Q = - alpha_2 * alpha / 108.0Q + alpha * gamma / 3.0Q - beta * beta / 8.0Q;

   __complex128 R = - Q / 2.0Q + csqrtq(Q * Q / 4.0Q + P * P * P / 27.0Q); // taken + sqrt

   __complex128 U = cpowq(R, 1.0Q / 3.0Q);

   __complex128 y_eq =- 5.0Q / 6.0Q * alpha + U - P / (3.0Q * U); // assuming that U !=0;
   __complex128 W = csqrtq(alpha + 2.0Q * y_eq);

   __complex128 back[4] =     {-1.0Q,  -1.0Q, -1.0Q, -1.0Q};
   __complex128 eff[4]  =     {10.0Q,  10.0Q, 10.0Q, 10.0Q};
   __float128 minusPlusS[4] = { 1.0Q,   1.0Q, -1.0Q, -1.0Q};
   __float128 minusPlusT[4] = { 1.0Q,  -1.0Q,  1.0Q, -1.0Q};

   __float128 maxLL = (-1.0Q) * FLT128_MAX; 
   __float128 backResult = 0;
   __float128 effResult = 0;

   // loop through signs
   for (int j = 0; j < 4; j ++) {


      back[j] = - B / (4.0Q * A) + (minusPlusS[j] * W + minusPlusT[j] *
       csqrtq(-(3.0Q * alpha + 2.0Q * y_eq + minusPlusS[j] * 2.0Q * beta / W))) / 2.0Q;
      eff[j] = CalcEff(mu, back[j], x, y, z, tau, m);

      __complex128 poly4compVal = poly4comp(back[j]);

      if (crealq(eff[j]) <= 1.0Q && crealq(eff[j]) >= 0.0Q && crealq(back[j]) >= 0.0Q) {
         // std::cout << "LikeMod1 " << j << ": " << std::endl;
         __float128 currentLL = LikeMod1(mu, crealq(back[j]), crealq(eff[j]), x, y, z, tau, m);
         if (currentLL >= maxLL) {
            backResult = crealq(back[j]);
            effResult = crealq(eff[j]);
            maxLL = currentLL;
            
         } else {
            std::cout << "CurrentLL < maxLL" << std::endl;
         }
         
      } else {
         std::cout << "e or/and b are out of range" << std::endl;
      }

   }

   std::cout << "backResult: " << (double) backResult << " effResult: " << (double) effResult << std::endl;
   b = backResult;
   e = effResult;
}

__complex128 TRolke128_num::CalcEff(__float128 mu, __complex128 b, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m) 
{
   __complex128 R = x + y - b * (1.0Q + tau);
   __complex128 Q = b * (1.0Q + tau) - y;
   return b * R / (mu * Q);
}
 
////////////////////////////////////////////////////////////////////////////////
/// Gradient model likelihood
 
__float128 TRolke128_num::LikeGradMod1(__float128 e, __float128 mu, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m)
{
   __float128 eta, etaprime, bprime, f;
   eta = static_cast<__float128>(z) / e - static_cast<__float128>(m - z) / (1.0Q - e);
   etaprime = (-1.0Q) * (static_cast<__float128>(m - z) / ((1.0Q - e) * (1.0Q - e)) + static_cast<__float128>(z) / (e * e));
   __float128 b = y / (tau - eta / mu); // YM: why is minus here?
   bprime = (b * b * etaprime) / mu / ((__float128) y);  // YM: why is absence of minus here?
   f = (mu + bprime) * (x / (e * mu + b) - 1) + (y / b - tau) * bprime + eta;
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 2:
/// Poisson background/ Gauss Efficiency
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod2(__float128 mu, Int_t x, Int_t y, __float128 em, __float128 sde, __float128 tau, Int_t what)
{
   __float128 v =  sde * sde;
   __float128 coef[4], roots[3];
   __float128 f = 0.0Q;
 
   if (what == 1) {
      f = x * (1.0Q - (y / tau) / (__float128) x) / em;
   }
 
   if (what == 2) {
      mu = x * (1.0Q - (y / tau) / (__float128) x) / em;
      __float128 b = y / tau;
      __float128 e = em;
      f = LikeMod2(mu, b, e, x, y, em, tau, v);
   }
 
   if (what == 3) {
      if (mu == 0) {
         __float128 b = (x + y) / (1.0Q + tau);
         __float128 e = em ;
         f = LikeMod2(mu, b, e, x, y, em, tau, v);
      } else {
         coef[3] = mu;
         coef[2] = mu * mu * v - 2.0Q * em * mu - mu * mu * v * tau;
         coef[1] = (- x) * mu * v - mu * mu * mu * v * v * tau - mu * mu * v * em + em * mu * mu * v * tau + em * em * mu - y * mu * v;
         coef[0] = x * mu * mu * v * v * tau + x * em * mu * v - y * mu * mu * v * v + y * em * mu * v;
 
         RootsCubic128(coef, roots[0], roots[1], roots[2]);
 
         __float128 e = roots[1];
         __float128 b;
         if ( v > 0) b = y / (tau + (em - e) / mu / v);
         else b = y/tau;
         f = LikeMod2(mu, b, e, x, y, em, tau, v);
      }
   }
 
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 2:
/// Poisson background/Gauss Efficiency
 
__float128 TRolke128_num::LikeMod2(__float128 mu, __float128 b, __float128 e, Int_t x, Int_t y, __float128 em, __float128 tau, __float128 v)
{
   __float128 s = e*mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
   __float128 bg = tau*b;
   __float128 llb =  -bg;
   if ( y > 0) llb =  y*Log128( bg) - bg - LogFactorial(y);
   __float128 lle = 0;
   if ( v > 0) lle = - 0.9189385Q - Log128(v) / 2.0Q - (em - e)*(em - e) / v / 2.0Q;
 
   return 2*( lls + llb + lle);
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 3:
/// Gauss  background/ Gauss Efficiency
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod3(__float128 mu, Int_t x, __float128 bm, __float128 em, __float128 sde, __float128 sdb, Int_t what)
{
   __float128 f = 0.;
   __float128  v = sde * sde;
   __float128  u = sdb * sdb;
 
   if (what == 1) {
      f = (x - bm) / em;
   }
 
 
   if (what == 2) {
      mu = (x - bm) / em;
      __float128 b  = bm;
      __float128 e  = em;
      f  = LikeMod3(mu, b, e, x, bm, em, u, v);
   }
 
 
   if (what == 3) {
      if (mu == 0.0) {
         __float128 b = ((bm - u) + Sqrt128((bm - u) * (bm - u) + 4 * x * u)) / 2.;
         __float128 e = em;
         f = LikeMod3(mu, b, e, x, bm, em, u, v);
      } else {
         __float128 e = em;
         __float128 b = bm;
         if ( v > 0) {
            __float128 temp[3];
            temp[0] = mu * mu * v + u;
            temp[1] = mu * mu * mu * v * v + mu * v * u - mu * mu * v * em + mu * v * bm - 2.0Q * u * em;
            temp[2] = mu * mu * v * v * bm - mu * v * u * em - mu * v * bm * em + u * em * em - mu * mu * v * v * x;
            e = (-temp[1] + Sqrt128(temp[1] * temp[1] - 4 * temp[0] * temp[2])) / 2.0Q / temp[0];
            b = bm - (u * (em - e)) / v / mu;
         }
         f = LikeMod3(mu, b, e, x, bm, em, u, v);
      }
   }
 
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 3:
/// Gauss background/Gauss Efficiency
 
__float128 TRolke128_num::LikeMod3(__float128 mu, __float128 b, __float128 e, Int_t x, __float128 bm, __float128 em, __float128 u, __float128 v)
{
   __float128 s = e*mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
   __float128 llb =  0.0Q;
   if ( u > 0) llb = - 0.9189385Q - Log128(u) / 2.0Q - (bm - b)*(bm - b) / u / 2.0Q;
   __float128 lle = 0.0Q;
   if ( v > 0) lle = - 0.9189385Q - Log128(v) / 2.0Q - (em - e)*(em - e) / v / 2.0Q;
 
   return 2.0Q*( lls + llb + lle);
 
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 4:
/// Poiss  background/Efficiency known
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod4(__float128 mu, Int_t x, Int_t y, __float128 tau, Int_t what)
{
   __float128 f = 0.0;
 
   if (what == 1) f = x - y / tau;
   if (what == 2) {
      mu = x - y / tau;
      __float128 b  = y / tau;
      f  = LikeMod4(mu, b, x, y, tau);
   }
   if (what == 3) {
      if (mu == 0.0) {
         __float128 b = __float128(x + y) / (1.0Q + tau);
         f = LikeMod4(mu, b, x, y, tau);
      } else {
         __float128 b = (x + y - (1.0Q + tau) * mu + sqrtq(( (x + y - (1.0Q + tau) * mu) * (x + y - (1.0Q + tau) * mu) + 4.0Q * (1.0Q + tau) * y * mu) ) ) / 2.0Q / (1.0Q + tau);
         f = LikeMod4(mu, b, x, y, tau);
      }
   }
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 4:
/// Poiss background/Efficiency known
 
__float128 TRolke128_num::LikeMod4(__float128 mu, __float128 b, Int_t x, Int_t y, __float128 tau)
{
   __float128 s = mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
   __float128 bg = tau*b;
   __float128 llb =  -bg;
   if ( y > 0) llb =  y*Log128( bg) - bg - LogFactorial(y);
 
   return 2*( lls + llb);
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 5:
/// Gauss  background/Efficiency known
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod5(__float128 mu, Int_t x, __float128 bm, __float128 sdb, Int_t what)
{
   __float128 u = sdb * sdb;
   __float128 f = 0;
 
   if (what == 1) {
      f = x - bm;
   }
   if (what == 2) {
      mu = x - bm;
      __float128 b  = bm;
      f  = LikeMod5(mu, b, x, bm, u);
   }
 
   if (what == 3) {
      __float128 b = ((bm - u - mu) + Sqrt128((bm - u - mu) * (bm - u - mu) - 4 * (mu * u - mu * bm - u * x))) / 2;
      f = LikeMod5(mu, b, x, bm, u);
   }
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 5:
/// Gauss background/Efficiency known
 
__float128 TRolke128_num::LikeMod5(__float128 mu, __float128 b, Int_t x, __float128 bm, __float128 u)
{
   __float128 s = mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
   __float128 llb =  0;
   if ( u > 0) llb = - 0.9189385 - Log128(u) / 2.0Q - (bm - b)*(bm - b) / u / 2;
 
   return 2*( lls + llb);
}
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 6:
/// Background known/Efficiency binomial
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod6(__float128 mu, Int_t x, Int_t z, __float128 b, Int_t m, Int_t what)
{
   __float128 coef[4], roots[3];
   __float128 f = 0.;
   __float128 zm = __float128(z) / m;
 
   if (what == 1) {
      f = (x - b) / zm;
   }
 
   if (what == 2) {
      mu = (x - b) / zm;
      __float128 e  = zm;
      f  = LikeMod6(mu, b, e, x, z, m);
   }
   if (what == 3) {
      __float128 e;
      if (mu == 0) {
         e = zm;
      } else {
         coef[3] = mu * mu;
         coef[2] = mu * b - mu * x - mu * mu - mu * m;
         coef[1] = mu * x - mu * b + mu * z - m * b;
         coef[0] = b * z;
         RootsCubic128(coef, roots[0], roots[1], roots[2]);
         e = roots[1];
      }
      f = LikeMod6(mu, b, e, x, z, m);
   }
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 6:
/// background known/ Efficiency binomial
 
__float128 TRolke128_num::LikeMod6(__float128 mu, __float128 b, __float128 e, Int_t x, Int_t z, Int_t m)
{
   __float128 s = e*mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
 
   __float128 lle = 0;
   if (z == 0)        lle = m * Log128(1-e);
   else if ( z == m)  lle = m * Log128(e);
   else               lle =   z * Log128(e) + (m - z)*Log128(1 - e) + LogFactorial(m) - LogFactorial(m-z) - LogFactorial(z);
 
   return 2* (lls + lle);
}
 
 
////////////////////////////////////////////////////////////////////////////////
/// Calculates the Profile Likelihood for MODEL 7:
/// background known/Efficiency Gauss
///  - what = 1: Maximum likelihood estimate is returned
///  - what = 2: Profile Likelihood of Maximum Likelihood estimate is returned.
///  - what = 3: Profile Likelihood of Test hypothesis is returned
/// otherwise parameters as described in the beginning of the class)
 
__float128 TRolke128_num::EvalLikeMod7(__float128 mu, Int_t x, __float128 em, __float128 sde, __float128 b, Int_t what)
{
   __float128 v = sde * sde;
   __float128 f = 0.;
 
   if (what ==  1) {
      f = (x - b) / em;
   }
 
   if (what == 2) {
      mu = (x - b) / em;
      __float128 e  = em;
      f  = LikeMod7(mu, b, e, x, em, v);
   }
 
   if (what == 3) {
      __float128 e;
      if (mu == 0) {
         e = em;
      } else {
         e = (-(mu * em - b - mu * mu * v) - Sqrt128((mu * em - b - mu * mu * v) * (mu * em - b - mu * mu * v) + 4 * mu * (x * mu * v - mu * b * v + b * em))) / (- mu) / 2;
      }
      f = LikeMod7(mu, b, e, x, em, v);
   }
 
   return f;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Profile Likelihood function for MODEL 6:
/// background known/ Efficiency gaussian
 
__float128 TRolke128_num::LikeMod7(__float128 mu, __float128 b, __float128 e, Int_t x, __float128 em, __float128 v)
{
   __float128 s = e*mu+b;
   __float128 lls = - s;
   if (x > 0) lls = x*Log128(s) - s - LogFactorial(x);
 
   __float128 lle = 0;
   if ( v > 0) lle = - 0.9189385 - Log128(v) / 2 - (em - e)*(em - e) / v / 2;
 
   return 2*( lls + lle);
}
 
////////////////////////////////////////////////////////////////////////////////
/// Evaluate polynomial
 
__float128 TRolke128_num::EvalPolynomial(__float128 x, const Int_t  coef[], Int_t N)
{
   const Int_t   *p;
   p = coef;
   __float128 ans = *p++;
   Int_t i = N;
 
   do
      ans = ans * x  +  *p++;
   while (--i);
 
   return ans;
}
 
////////////////////////////////////////////////////////////////////////////////
/// Evaluate mononomial
 
__float128 TRolke128_num::EvalMonomial(__float128 x, const Int_t coef[], Int_t N)
{
   __float128 ans;
   const Int_t   *p;
 
   p   = coef;
   ans = x + *p++;
   Int_t i = N - 1;
 
   do
      ans = ans * x  + *p++;
   while (--i);
 
   return ans;
}
 
////////////////////////////////////////////////////////////////////////////////
/// LogFactorial function (use the logGamma function via the relation Gamma(n+1) = n!
 
__float128 TRolke128_num::LogFactorial(Int_t n)
{ // YM: seems working at our range
   if (TMath::LnGamma(n+1) < 0) exit(-1);
   return TMath::LnGamma(n+1);
}


   //TMath replacement with float128

__float128 TRolke128_num::Poisson128(__float128 x, __float128 par) {
   if (x<0)
      return 0;
   else if (x == 0.0)
      return 1./expq(par);
   else {
      __float128 lnpoisson = x*logq(par)-par-TMath::LnGamma(x+1.);
      return expq(lnpoisson);
   }
 }


 __float128 TRolke128_num::PoissonI128(__float128 x, __float128 par) {
      Int_t ix = Int_t(x);
      return Poisson128(ix,par);
 }
   // TMath::PoissonI(loop_x, background)

   // TMath::ChisquareQuantile
   // TMath::LnGamma

 __float128 TRolke128_num::Max128(__float128 a, __float128 b) {
   return (a > b) ? a : b;
 }
   // TMath::Max()

 __float128 TRolke128_num::Log128(__float128 a) {
   return logq(a);
 }
   // TMath::Log()
   
   
   
 __float128 TRolke128_num::Sqrt128(__float128 a) {
   return sqrtq(a);
 }
   // TMath::Sqrt()

 Bool_t TRolke128_num::RootsCubic128(const __float128 coef[4],__float128 &a, __float128 &b, __float128 &c) {
   Bool_t complex = kFALSE;
   __float128 r,s,t,p,q,d,ps3,ps33,qs2,u,v,tmp,lnu,lnv,su,sv,y1,y2,y3;
   a    = 0;
   b    = 0;
   c    = 0;
   if (coef[3] == 0) return complex;
   r    = coef[2]/coef[3];
   s    = coef[1]/coef[3];
   t    = coef[0]/coef[3];
   p    = s - (r*r)/3.0Q;
   ps3  = p/3.0Q;
   q    = (2.0Q*r*r*r)/27.0Q - (r*s)/3.0Q + t;
   qs2  = q/2.0Q;
   ps33 = ps3*ps3*ps3;
   d    = ps33 + qs2*qs2;
   if (d>=0) {
      complex = kTRUE;
      d   = sqrtq(d);
      u   = -qs2 + d;
      v   = -qs2 - d;
      tmp = 1.0Q/3.0Q;
      lnu = logq(fabsq(u));
      lnv = logq(fabsq(v));
      // su  = TMath::Sign(1.,u);
      // sv  = TMath::Sign(1.,v);

      if (u < 0) {
         su = -1.0Q;
      } else {
         su = 1;
      }

      if (v < 0) {
         sv = -1.0Q;
      } else {
         sv = 1.0Q;
      }

      u   = su*expq(tmp*lnu);
      v   = sv*expq(tmp*lnv);
      y1  = u + v;
      y2  = -y1/2.0Q;
      y3  = ((u-v)*sqrtq(3.0Q))/2.0Q;
      tmp = r/3.0Q;
      a   = y1 - tmp;
      b   = y2 - tmp;
      c   = y3;
   } else {
      __float128 phi,cphi,phis3,c1,c2,c3,pis3;
      ps3   = -ps3;
      ps33  = -ps33;
      cphi  = -qs2/sqrtq(ps33);
      phi   = acosq(cphi);
      phis3 = phi/3.0Q;
      pis3  = M_PIq/3.0Q;
      c1    = cosq(phis3);
      c2    = cosq(pis3 + phis3);
      c3    = cosq(pis3 - phis3);
      tmp   = sqrtq(ps3);
      y1    = 2.0Q*tmp*c1;
      y2    = -2.0Q*tmp*c2;
      y3    = -2.0Q*tmp*c3;
      tmp = r/3.0Q;
      a   = y1 - tmp;
      b   = y2 - tmp;
      c   = y3 - tmp;
   }
   return complex;
 }
   // TMath::RootsCubic