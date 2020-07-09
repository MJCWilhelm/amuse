//Implementation of rates


#include "rates.h"
#include "Common.h"

using namespace std;



/************************** Source spectra ***************************************/

// Planck curve divided by nu
double PlanckDivNu( double x, double T ){
  return pow( x, 2.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

//Planck curve
double Planck( double x, double T ){
  return pow( x, 3.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

// Integral of Planck curve
double PlanckInt( double _upperBound, double _T){
  return qromb( Planck, 1.0, _upperBound, _T);
}

/************************** Cross sections ***************************************/

// cross-section for H from Verner et al. 1996
double cross_HI( const double& _nu, const double& nu0HI ){

  return (_nu < nu0HI) ? 0.0 : 1.16032e43 * pow( 9.62231e-15 * _nu - 1.0, 2.) * 
    pow( _nu, -4.0185) * pow( 1.0 + 1.7107e-8 * pow( _nu, 0.5), -2.963 ) ;

}

// hydrogen cross section from Verner et al. 1996 times black body curve
double fHI( double x, double T ){

  //  return ( x < 1.0) ? 0.0 : 6.3 * pow( x, -1.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

  return ( x < 1.0) ? 0.0 : 
    1.16032e61 * 
    pow( 31.6379628 * x - 1.0, 2.) * 
    pow( nu0HI * x, -4.0185) * 
    pow( 1.0 + 1.7107e-8 * pow( nu0HI * x , 0.5), -2.963 ) * 
    pow( x, 2.0 ) / 
    ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

}

// Integral of Planck curve times cross section
double crossIntHI( double _upperBound, double _T ){
  return qromb(fHI, 1.0, _upperBound, _T);
}

/************************** Ionisations ***************************************/

//collisional ionisation coefficient from Theuns et al. 1998
double coll_ion_coeff_HI( const double& _T ){

   return 1.17e-10 * sqrt( _T ) * exp( -1.578091e5 / _T ) * pow(1.0 + sqrt(_T/1e5), -1.0);

}


/************************** Recombinations ************************************/

// Recombination coefficient HII for case B from Hui & Gnedin 1997
double recomb_coeff_HII_caseB( const double& tempGas ) {  

  // If no temperature is given, 3K is assumed
  return (tempGas>0.0) ? 2.753e-14 * pow( 2.0 * 1.57807e5/tempGas, 1.5 ) / pow( ( 1.0 + pow( 0.729927007 * 1.57807e5/tempGas, 0.407 )), 2.242 ) : 5.97797e-11;

}

// Recombination coefficient HII for case A from Hui & Gnedin 1997
double recomb_coeff_HII_caseA( const double& tempGas ) {  

  return 1.269e-13 * pow( 2.0 * 1.57807e5/tempGas, 1.503 ) / pow( ( 1.0 + pow( 3.831417625 * 1.57807e5/tempGas, 0.47 )), 1.923 );

}

// Recombination coefficient HII for grain-assisted recombination from Weingartner & Draine 2001
double recomb_coeff_HII_grain( const double& tempGas, const double& n_e, const double& G, const double& phi_pah ) {

  double psi = G/G0 * sqrt(tempGas)/((n_e > 0. ? n_e : 1. )*phi_pah);

  double temp = G > 0. ? 8.074e-6*pow(psi, 1.378)*(1.0 + 5.087e2*pow(tempGas, 1.586e-2)*pow(psi, -0.4723 - 1.102e-5*log(tempGas))) : 0.;

  return 1.e-14 * 12.25 * pow(1. + temp, -1.0);
}

/************************** Molecule Dissociation *****************************/

// Photodissociation of H2 by FUV radiation from Draine & Bertholdi 1996
// G is the local FUV field in Draine units; note that self-shielding and absorption effects are handled by simplex
double photodissoc_coeff_H2( const double& G ) {

  return 3.3e-11 * G/G0;
}

// Collisional dissociation of H2 by collision with H atoms from Glover & Mac Low 2007
double dissoc_coeff_H_atom( const double& densGas, const double& tempGas, const double& dens_cr ) {

  double x = densGas/dens_cr;
  double kHh = 3.52e-9 * exp(-(4.39e4)/tempGas);
  double kHl = 1.2e-16 * sqrt(tempGas/4500.);

  double logkH = (x/(1.0 + x))*log10(kHh) + (1.0/(1.0 + x))*log10(kHl);

  return pow(10., logkH);
}

// Collisional dissociation of H2 by collision with H2 molecules from Glover & Mac Low 2007
double dissoc_coeff_H_molecule( const double& densGas, const double& tempGas, const double& dens_cr ) {

  double x = densGas/dens_cr;
  double kH2h = 1.30e-9 * exp(-(5.33e4)/tempGas);
  double kH2l = interpolate_kH2l(tempGas);

  double logkH2 = (x/(1.0 + x))*log10(kH2h) + (1.0/(1.0 + x))*log10(kH2l);

  return pow(10., logkH2);
}

// Helper function for the above functions, computes the collision-induced dissociation of H2 
// with other H2 by interpolating table 2 from Martin, Keogh & Mandy 1998
double interpolate_kH2l( const double& tempGas ) {

  int N = 11;
  double logT [N] = {log10(1000.), log10(1400.), log10(2000.), log10(3000.), log10(4500.), log10(6000.), log10(7800.), log10(10000.), log10(14000.), log10(20000.), log10(30000.)};
  double logA [N] = {log10(3.8e-41), log10(9.4e-34), log10(5.0e-28), log10(2.4e-23), log10(5.4e-20), log10(3.5e-18), log10(8.1e-17), log10(9.9e-16), log10(1.7e-14), log10(2.0e-13), log10(1.9e-12)};

  double logTgas = log10(tempGas);
  double a, b;

  if (logTgas < logT[0]) { return pow(10., logA[0]); }

  for (int i = 0; i < N-1; i++) {
    if (logTgas > logT[i] && logTgas < logT[i+1]) {
      b = (logA[i+1] - logA[i])/(logT[i+1] - logT[i]);
      a = logA[i+1] - b*logT[i+1];
      return pow(10., a + b*logTgas);
    }
  }

  return pow(10., logA[N-1]);
}

// Helper function for the above functions, to compute dens_cr
double critical_density( const double& xH, const double& xH2, const double& tempGas ) {

  double x = log10(tempGas/1e4);
  double dens_cr_H  = pow(10.0, 4.00  - 0.416*x - 0.327*x*x)/10.;
  double dens_cr_H2 = pow(10.0, 4.845 - 1.3*x + 1.62*x*x);

  return 1.0/(xH/dens_cr_H + xH2/dens_cr_H2);
}

// Dissociation coefficient of CO from Nelson & Langer 1997
double dissoc_coeff_CO ( const double& G ) {

  return G/G0 * 1e-10;
}

/************************** Molecule Formation ********************************/

// Formation coefficient of molecular hydrogen on dust grains from Hollenbach & McKee 1979
double formation_coeff_H2_grain( const double& tempGas, const double& tempDust ) {

  // Values in Hollenbach & McKee 1979
  double NtNi = 1e4; // Sect III-b
  double tempCrit = 65.; // Sect III-b
  //double tempDust = 100.; // Sect III-c

  double fa = tempDust > 0.0 ? 1.0/(1.0 + exp(log(NtNi)*(1.0 - tempCrit/tempDust))) : 1.0;
  return 3e-17 * pow(tempGas/1e2, 0.5) * fa / (1.0 + 0.4*pow(tempGas/1e2 + tempDust/1e2, 0.5) + 0.2*tempGas/1e2 + 0.08*pow(tempGas/1e2, 2));
}

// Formation coefficient of CO from Nelson & Langer 1997
double formation_coeff_CO( const double& densGas, const double& densH2, const double& densO, const double& G ) {

    double beta = 0.;
    if (densH2 > 0. && densO > 0.)
        beta = 1./( 1. + ( G/G0 / densH2 )/( densO/densGas ) );

    return 5e-16 * beta;
}

/************************** Heating *******************************************/

//Planck curve times H cross section times h*(nu-nu_0) and divided by h*nu
//Verner et al. 1996 cross section is used
double enerfHI( double x, double T ){

  return ( x < 1.0) ? 0.0 :
    ( x - 1.0) *
    1.16032e61 * 
    pow( 31.6379628 * x - 1.0, 2.0) * 
    pow( nu0HI * x, -4.0185) * 
    pow( 1.0 + 1.7107e-8 * pow( nu0HI * x , 0.5), -2.963 ) * 
    pow( x, 2.0 ) / 
    ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 ); 

}


// Photoelectric heating from UV irradiated grains and PAHs from Bakes & Tielens 1994 and Wolfire et al. 2003
double photoelectric_heating_coeff( const double& densElectron, const double& tempGas, const double& phi_pah, const double& G ) {

  double psi = G/G0 * sqrt(tempGas)/((densElectron > 0. ? densElectron : 1. )*phi_pah);
  double epsilon = 4.9e-2/(1.0 + 4e-3*pow(psi, 0.73)) + 3.7e-2*pow(tempGas/1e4, 0.7)/(1.0 + 2e-4*psi);

  return 1.3e-24 * epsilon * G/G0;
}

// Heating by photodissociation of H2 from Black & Dalgarno 1977 (see also note of photodissoc_coeff_H2)
double H2_photodissocation_heating_coeff( const double& G ) {

  return 6.4e-13 * photodissoc_coeff_H2( G );
}

// Heating by radiative pumping of H2 from Glover and Mac Low 2007
// We use the scaling faction 8.5, with an additional scaling factor 
// to go from 0.4 eV per deposition for dissociation to whatever the scaling deposition is
double H2_pumping_heating_coeff( const double& densGas, const double& dens_cr, const double& G ) {

  double excitation_energy = 2./(1. + dens_cr/densGas);
  return 8.5 * H2_photodissocation_heating_coeff( G ) * (excitation_energy / 0.4);
}

// Molecular hydrogen formation heating from Glover & Mac Low 2007
double H2_formation_heating_coeff( const double& densGas, const double& tempGas, const double& dens_cr, const double& tempDust ) {

  double RH2 = formation_coeff_H2_grain( tempGas, tempDust );
  return 7.2e-12 * RH2 / (1.0 + dens_cr/densGas);
}

// Cosmic ray heating from Goldsmith & Langer 1978
double cosmic_ray_heating_coeff( const double& cr_ir ) {

  return 3.2e-28 * (cr_ir/1e-17);
}


/************************** Cooling ******************************************/

// Recombination cooling coefficient HII for case B from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseB( const double& tempGas ) {  

  return 3.435e-30 * tempGas * pow( 2.0 * 1.57807e5/tempGas, 1.970 ) / 
    pow( ( 1.0 + pow( 0.88888888889 * 1.57807e5/tempGas, 0.376 )), 3.720 );

}

// Recombination cooling coefficient HII for case A from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseA( const double& tempGas ) {  

  return 1.778e-29 * tempGas * pow( 2.0 * 1.57807e5/tempGas, 1.965 ) / 
    pow( ( 1.0 + pow( 3.696857671 * 1.57807e5/tempGas, 0.502 )), 2.697 );

}

// Collisional ionisation cooling coefficient from Theuns et al. 1998
double coll_ion_cooling_coeff_HI( const double& _T ){
  
  return 2.54e-21 * sqrt( _T ) * exp( -1.578091e5 / _T ) * pow(1.0 + sqrt(_T/1e5), -1.0);

}

// Collisional excitation cooling from Cen 1992
double coll_excit_cooling_coeff_HI( const double& _T ){
  
  return 7.5e-19 * pow(1.0 + sqrt(_T/1.e5), -1.0) * exp( -1.18348e5 / _T);

}

//free-free cooling coefficient with gaunt factor from Theuns et al. 1998
double ff_cooling_coeff( const double& tempGas ){  

  //gaunt factor  
  double g_ff = 1.1 + 0.34 * exp( -0.33333333333 * pow(5.5-log(tempGas)/log(10.), 2.0) );

  return 1.42e-27 * g_ff * sqrt( tempGas );// This must be multiplied with (n_HII + n_HeII + 4 n_HeIII) * n_e

}

// Grain recombination cooling from Wolfire et al. 2003
double recomb_cooling_coeff_HII_grain( const double& densElectron, const double& tempGas, const double& phi_pah, const double& G ) {

  double beta = 0.74/pow(tempGas, 0.068);
  double psi = G/G0*sqrt(tempGas)/((densElectron > 0. ? densElectron : 1. )*phi_pah);

  return 4.65e-30 * phi_pah * pow(tempGas, 0.94) * pow(psi, beta);

}

// Dust-gas energy transfer
double grain_transfer_cooling_coeff( const double& densGas, const double& tempGas, const double& tempDust, const double& amin ) {

  return 1.2e-31 * densGas*densGas * sqrt(tempGas/1000.) * sqrt(1e-6/amin) * (1.0 - 0.8*exp(-75./tempGas)) * (tempGas - tempDust);
}

//H2 rovibrational cooling from Coppola et al. 2019
double H2_cooling_coeff( const double& tempGas, const double& densGas){

  double temp = 0.0;
  double logT = log10(tempGas);
  double logn = log10(densGas);

  //if (tempGas < 100. || tempGas > 4000.) {
  //  cout << "Temperature out of bounds (" << tempGas << " K), extrapolating fit\nExtrapolation likely safe" << endl;
  //}
  if (densGas < 1e-4 || densGas > 1e8) {
    cout << "Density out of bounds (" << densGas << " cm^-3), extrapolating fit\nExtrapolation potentially unsafe" << endl;
  }

  temp += -1.07761178e2;
  temp +=  1.17901741e2*logT;
  temp += -6.61721218e1*pow(logT, 2);
  temp +=  1.67152116e1*pow(logT, 3);
  temp += -1.55475014e0*pow(logT, 4);

  temp +=  8.50704285e0*logn;
  temp += -1.09489742e1*logT*logn;
  temp +=  5.95232129e0*pow(logT, 2)*logn;
  temp += -1.43513388e0*pow(logT, 3)*logn;
  temp +=  1.29639563e-1*pow(logT, 4)*logn;

  temp += -3.08850636e-1*pow(logn, 2);
  temp +=  2.51708582e-1*logT*pow(logn, 2);
  temp += -1.10639748e-1*pow(logT, 2)*pow(logn, 2);
  temp +=  2.84953800e-2*pow(logT, 3)*pow(logn, 2);
  temp += -3.13943372e-3*pow(logT, 4)*pow(logn, 2);

  temp += -3.97357274e-1*pow(logn, 3);
  temp +=  5.47452438e-1*logT*pow(logn, 3);
  temp += -2.92423859e-1*pow(logT, 2)*pow(logn, 3);
  temp +=  7.02040444e-2*pow(logT, 3)*pow(logn, 3);
  temp += -6.36261310e-3*pow(logT, 4)*pow(logn, 3);

  temp +=  3.65589231e-2*pow(logn, 4);
  temp += -4.89055533e-2*logT*pow(logn, 4);
  temp +=  2.57375696e-2*pow(logT, 2)*pow(logn, 4);
  temp += -6.18478233e-3*pow(logT, 3)*pow(logn, 4);
  temp +=  5.65797161e-4*pow(logT, 4)*pow(logn, 4);

  return pow(10., temp);
}

// CO cooling rate from Whitworth & Jaffa 2018
// replaced 'velocity divergence' with 'inverse free-fall time' times a turbulence factor
// (basically Mach number) as code is hydrodynamics-agnostic
double CO_cooling_coeff( const double& densGas, const double& densH2, const double& densCO, const double& tempGas, const double& t_ff, const double& turbulence) {


  double LabdaLO = 2.16e-27 * densH2 * pow(tempGas, 3./2.);
  double LabdaHI = 0.;
  if (densH2 > 0. && densCO > 0.)
    LabdaHI = 2.21e-28 / densH2 * pow(tempGas, 4.0) / ( densCO/densGas * t_ff / (megaParsecToCm * 1e-11 ) ) * turbulence;
  double beta = 1.;
  if (densH2 > 0.)
    beta = 1.23 * pow(densH2, 0.0533) * pow(tempGas, 0.164);

  return pow(pow(LabdaLO, -1.0/beta) + pow(LabdaHI, -1.0/beta), -beta);
}

// Collisional dissociation cooling rate of H2, with both H and H2 from Lepp & Shull 1983
double H2_coll_dissoc_cooling_coeff( const double& densGas, const double& densH2, const double& densH, const double& tempGas ) {

  double dens_cr = critical_density( densH/densGas, 2.0*densH2/densGas, tempGas );

  double LvH = 1.1e-13 * exp(-6744./tempGas);
  double LvL = (densH*dissoc_coeff_H_atom( densGas, tempGas, dens_cr ) + densH2*dissoc_coeff_H_molecule( densGas, tempGas, dens_cr ))*(8.18e-13);

  double x = log10(tempGas/1e4);
  double LrH =   tempGas < 1087. ? pow(10., -19.24 + 0.474*x - 1.247*x*x) : 3.90e-19 * exp(-6118./tempGas);
  double LrL = ( tempGas < 4031. ? pow(10., -22.9  - 0.553*x - 1.148*x*x) : 1.38e-22 * exp(-9243./tempGas) )*(pow(densH2, 0.77) + 1.2*pow(densH, 0.77));

  return LvH/(1. + LvH/LvL) + LrH/(1. + LrH/LrL);
}
