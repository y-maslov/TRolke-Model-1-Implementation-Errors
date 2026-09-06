//////////////////////////////////////////////////////////////////////////////
//
//  TRolke128_num
//
//  This class computes confidence intervals for the rate of a Poisson
//  in the presence of background and efficiency with a fully frequentist
//  treatment of the uncertainties in the efficiency and background estimate
//  using the profile likelihood method.
//
//      Author: Jan Conrad (CERN) 2004
//      Updated: Johan Lundberg (CERN) 2009
//
//      Copyright CERN 2004,2009           Jan.Conrad@cern.ch,
//                                     Johan.Lundberg@cern.ch
//
//  For information about the statistical meaning of the parameters
//  and the syntax, consult TRolke128_num.cxx
//                  ------------------
//
//  Examples are found in the file Rolke.C
//  --------------------------------------
//
//////////////////////////////////////////////////////////////////////////////
 
#ifndef ROOT_TRolke128_num
#define ROOT_TRolke128_num
 
#include "TObject.h"
#include "TMath.h"
#include <quadmath.h>
 
// Class definition. This class is not intended to be used as a base class.
class TRolke128_num //: public TObject
{
 
private:
   __float128 fCL;         // confidence level as a fraction [0.9 for 90% ]
   __float128 fUpperLimit; // the calculated upper limit
   __float128 fLowerLimit; // the calculated lower limit
   bool fBounding;       // false for unbounded likelihood
                         // true for bounded likelihood
   Int_t fNumWarningsDeprecated1;
   Int_t fNumWarningsDeprecated2;
 
   /* ----------------------------------------------------------------- */
   /* These variables are set by the Set methods for the various models */
   Int_t f_x;
   Int_t f_y;
   Int_t f_z;
   __float128 f_bm;
   __float128 f_em;
   __float128 f_e;
   Int_t f_mid;
   __float128 f_sde;
   __float128 f_sdb;
   __float128 f_tau;
   __float128 f_b;
   Int_t f_m;
 
   /* ----------------------------------------------------------------- */
   /* Internal helper functions and methods */
   // The Calculator
   __float128 Interval(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m);
 
   // LIKELIHOOD ROUTINE
   __float128 Likelihood(__float128 mu, Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m, Int_t what);
 
   //MODEL 1
   // __float128 EvalLikeMod1(__float128 mu, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m, Int_t what);
   __float128 LikeMod1(__float128 mu, __float128 b, __float128 e, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m);
   void     ProfLikeMod1(__float128 mu, __float128 &b, __float128 &e, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m);
   void     ProfLikeMod1Fix(__float128 mu, __float128 &b, __float128 &e, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m);
   __complex128 CalcEff(__float128 mu, __complex128 b, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m);
   __float128 LikeGradMod1(__float128 e, __float128 mu, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m);
 
   //MODEL 2
   __float128 EvalLikeMod2(__float128 mu, Int_t x, Int_t y, __float128 em, __float128 sde, __float128 tau, Int_t what);
 
   __float128 LikeMod2(__float128 mu, __float128 b, __float128 e, Int_t x, Int_t y, __float128 em, __float128 tau, __float128 v);
 
   //MODEL 3
   __float128 EvalLikeMod3(__float128 mu, Int_t x, __float128 bm, __float128 em, __float128 sde, __float128 sdb, Int_t what);
   __float128 LikeMod3(__float128 mu, __float128 b, __float128 e, Int_t x, __float128 bm, __float128 em, __float128 u, __float128 v);
 
   //MODEL 4
   __float128 EvalLikeMod4(__float128 mu, Int_t x, Int_t y, __float128 tau, Int_t what);
   __float128 LikeMod4(__float128 mu, __float128 b, Int_t x, Int_t y, __float128 tau);
 
   //MODEL 5
   __float128 EvalLikeMod5(__float128 mu, Int_t x, __float128 bm, __float128 sdb, Int_t what);
   __float128 LikeMod5(__float128 mu, __float128 b, Int_t x, __float128 bm, __float128 u);
 
   //MODEL 6
   __float128 EvalLikeMod6(__float128 mu, Int_t x, Int_t z, __float128 b, Int_t m, Int_t what);
   __float128 LikeMod6(__float128 mu, __float128 b, __float128 e, Int_t x, Int_t z, Int_t m);
 
   //MODEL 7
   __float128 EvalLikeMod7(__float128 mu, Int_t x, __float128 em, __float128 sde, __float128 b, Int_t what);
   __float128 LikeMod7(__float128 mu, __float128 b, __float128 e, Int_t x, __float128 em, __float128 v);
 
   //MISC
   static __float128 EvalPolynomial(__float128 x, const Int_t coef[], Int_t N);
   static __float128 EvalMonomial(__float128 x, const Int_t coef[], Int_t N);
   __float128 LogFactorial(Int_t n);
 
   __float128 ComputeInterval(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m);
 
   void SetModelParameters(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m);
 
   void SetModelParameters();
 
   __float128 GetBackground();
 
public:
 
   __float128 EvalLikeMod1(__float128 mu, Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m, Int_t what);

   /* Constructor */
   TRolke128_num(__float128 CL = 0.9, Option_t *option = "");
 
   /* Destructor */
   ~TRolke128_num();
 
   /* Get and set the Confidence Level */
   __float128 GetCL() const         {
      return fCL;
   }
   void     SetCL(__float128 CL)  {
      fCL = CL;
   }
 
   /* Set the Confidence Level in terms of Sigmas. */
   void SetCLSigmas(__float128 CLsigmas) {
      fCL = TMath::Erf(CLsigmas / TMath::Sqrt(2.0)) ;
   }
 
   // The Set methods for the different models are described in Rolke.cxx
   // model 1
   void SetPoissonBkgBinomEff(Int_t x, Int_t y, Int_t z, __float128 tau, Int_t m);
 
   // model 2
   void SetPoissonBkgGaussEff(Int_t x, Int_t y, __float128 em, __float128 tau, __float128 sde);
 
   // model 3
   void SetGaussBkgGaussEff(Int_t x, __float128 bm, __float128 em, __float128 sde, __float128 sdb);
 
   // model 4
   void SetPoissonBkgKnownEff(Int_t x, Int_t y, __float128 tau, __float128 e);
 
   // model 5
   void SetGaussBkgKnownEff(Int_t x, __float128 bm, __float128 sdb, __float128 e);
 
   // model 6
   void SetKnownBkgBinomEff(Int_t x, Int_t z, Int_t m, __float128 b);
 
   // model 7
   void SetKnownBkgGaussEff(Int_t x, __float128 em, __float128 sde, __float128 b);
 
   /* Deprecated interface method (read Rolke.cxx). May be removed from future releases */
   __float128 CalculateInterval(Int_t x, Int_t y, Int_t z, __float128 bm, __float128 em, __float128 e, Int_t mid, __float128 sde, __float128 sdb, __float128 tau, __float128 b, Int_t m);
 
   // get the upper and lower limits based on the specified model
   bool GetLimits(__float128& low, __float128& high);
   __float128 GetUpperLimit();
   __float128 GetLowerLimit();
 
   // get the upper and lower average limits
   bool GetSensitivity(__float128& low, __float128& high, __float128 pPrecision = 0.00001);
 
   // get the upper and lower limits for the outcome corresponding to
   // a given quantile.
   bool GetLimitsQuantile(__float128& low, __float128& high, Int_t& out_x, __float128 integral = 0.5);
 
   // get the upper and lower limits for the most likely outcome.
   bool GetLimitsML(__float128& low, __float128& high, Int_t& out_x);
 
   // get the value of x corresponding to rejection of the null hypothesis.
   bool GetCriticalNumber(Int_t& ncrit,Int_t maxtry=-1);
 
   /* Get the bounding mode flag. True activates bounded mode. Read
      TRolke128_num.cxx and the references therein for details. */
   bool GetBounding() const {
      return fBounding;
   }
 
   /* Get the bounding mode flag. True activates bounded mode. Read
      TRolke128_num.cxx and the references therein for details. */
   void SetBounding(const bool bnd) {
      fBounding = bnd;
   }
 
   /* Deprecated name for SetBounding. */
   void SetSwitch(bool bnd) ;
 
   /* Dump internals. Option is not used */
   void Print(Option_t*) const;
 
   //TMath replacement with float128

   __float128 Poisson128(__float128 x, __float128 par);

   __float128 PoissonI128(__float128 x, __float128 par);
   // TMath::PoissonI(loop_x, background);

   // TMath::ChisquareQuantile
   // TMath::LnGamma

   __float128 Max128(__float128 a, __float128 b);
   // TMath::Max()

   __float128 Log128(__float128 a);
   // TMath::Log()
   
   
   
   __float128 Sqrt128(__float128 a);
   // TMath::Sqrt()

   Bool_t RootsCubic128(const __float128 coef[4],__float128 &a, __float128 &b, __float128 &c);
   // TMath::RootsCubic

};
 

//calculate confidence limits using the Rolke method
#endif
 