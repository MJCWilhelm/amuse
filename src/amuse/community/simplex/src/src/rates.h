
//this header contains the rates and cross sections needed by SimpleX.cpp

#ifndef RATES_H
#define RATES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <float.h>


using namespace std;

/************************** Source spectra ***************************************/

//! Planck curve divided by frequency
double PlanckDivNu( double x, double T );

//! Planck curve
double Planck( double x, double T );

//! Integral of Planck curve
double PlanckInt( double _upperBound, double _T);


/************************** Cross sections ***************************************/

// The Planck curve and cross section have been scaled for stable numerics
// constants in Planck curve do not matter as they are divided out by normalisation
// Multiply any cross section with 1e-18 to obtain physical values


//! cross-section for H from Verner et al. 1996
double cross_HI( const double& _nu, const double& nu0HI );

//! Planck curve times cross section and divided by nu
double fHI( double x, double T );

//! Integral of Planck curve times cross section
double crossIntHI( double _upperBound, double _T );


/************************** Ionisations ***************************************/


//! collisional ionisation coefficient of hydrogen from Theuns et al. 1998
double coll_ion_coeff_HI( const double& _T );

double Xray_ionisations (const double& Nw, const double& n_e, const double& n_H);


/************************** Recombinations ************************************/

//! Recombination coefficient HII for case B from Hui & Gnedin 1997
double recomb_coeff_HII_caseB( const double& tempGas );

//! Recombination coefficient HII for case A from Hui & Gnedin 1997
double recomb_coeff_HII_caseA( const double& tempGas );

// Recombination coefficient HII for grain-assisted recombination from Weingartner & Draine 2001
double recomb_coeff_HII_grain( const double& tempGas, const double& n_e, const double& G, const double& phi_pah );



/************************** Molecule Dissociation *****************************/

// Photodissociation of H2 by FUV radiation from Draine & Bertholdi 1996
// G is the local FUV field in Draine units; note that self-shielding and absorption effects are handled by simplex
double photodissoc_coeff_H2( const double& G );

// Collisional dissociation of H2 by collision with H atoms from Glover & Mac Low 2007
double dissoc_coeff_H_atom( const double& densGas, const double& tempGas, const double& dens_cr );

// Collisional dissociation of H2 by collision with H2 molecules from Glover & Mac Low 2007
double dissoc_coeff_H_molecule( const double& densGas, const double& tempGas, const double& dens_cr );

// Helper function for the above functions, computes the collision-induced dissociation of H2 
// with other H2 by interpolating table 2 from Martin, Keogh & Mandy 1998
double interpolate_kH2l( const double& tempGas );

// Helper function for the above functions, to compute dens_cr
double critical_density( const double& xH, const double& xH2, const double& tempGas );

// Photodissociation of CO by FUV radiation from Nelson & Langer 1997
double dissoc_coeff_CO ( const double& G );



/************************** Molecule Formation ********************************/

// Formation coefficient of molecular hydrogen on dust grains from Hollenbach & McKee 1979
double formation_coeff_H2_grain( const double& tempGas, const double& tempDust );

// Formation coefficient of CO from Nelson & Langer 1997
double formation_coeff_CO( const double& densGas, const double& densH2, const double& densO, const double& G );



/************************** Heating *******************************************/


// The Planck curve and cross section have been scaled for stable numerics
// ants in Planck curve do not matter as they are divided out by normalisation 
// except for one factor of h
// Multiply any cross section with 1e-18 to obtain physical values

//! Planck curve times H cross section times h*(nu-nu_i) and divided by h*nu
//! Verner et al. 1996 cross section is used
double enerfHI( double x, double T );

// Photoelectric heating from UV irradiated grains and PAHs from Bakes & Tielens 1994 and Wolfire et al. 2003
double photoelectric_heating_coeff( const double& densElectron, const double& tempGas, const double& phi_pah, const double& G );

// Heating by photodissociation of H2 from Black & Dalgarno 1977 (see also note of photodissoc_coeff_H2)
double H2_photodissocation_heating_coeff( const double& G );

// Heating by radiative pumping of H2 from Glover and Mac Low 2007
// We use the scaling faction 8.5, with an additional scaling factor 
// to go from 0.4 eV per deposition for dissociation to whatever the scaling deposition is
double H2_pumping_heating_coeff( const double& densGas, const double& dens_cr, const double& G );

// Molecular hydrogen formation heating from Glover & Mac Low 2007
double H2_formation_heating_coeff( const double& densGas, const double& tempGas, const double& dens_cr, const double& tempDust );

// Cosmic ray heating from Goldsmith & Langer 1978
double cosmic_ray_heating_coeff( const double& cr_ir );

double Xray_heating_coeff(const double& Nw, const double& n_e, const double& n_H);



/************************** Cooling ******************************************/


//! Recombination cooling coefficient HII for case B from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseB( const double& tempGas );

//! Recombination cooling coefficient HII for case A from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseA( const double& tempGas );

//! Collisional ionisation cooling coefficient of hydrogen from Theuns et al. 1998
double coll_ion_cooling_coeff_HI( const double& _T );

//! Collisional excitation cooling coefficient of hydrogen from Cen 1992
double coll_excit_cooling_coeff_HI( const double& _T );

//! Free-free cooling coefficient with gaunt factor from Theuns et al. 1998
double ff_cooling_coeff( const double& tempGas );

// Grain recombination cooling from Wolfire et al. 2003
double recomb_cooling_coeff_HII_grain( const double& densGas, const double& densElectron, const double& tempGas, const double& phi_pah );

// Dust-gas energy transfer
double grain_transfer_cooling_coeff( const double& densGas, const double& tempGas, const double& tempDust, const double& amin );

//H2 rovibrational cooling from Coppola et al. 2019
double H2_cooling_coeff( const double& tempGas, const double& densGas);

// CO cooling rate from Whitworth & Jaffa 2018
double CO_cooling_coeff( const double& densGas, const double& densH2, const double& densCO, const double& tempGas, const double& t_ff, const double& turbulence);

// Collisional dissociation cooling rate of H2, with both H and H2 from Lepp & Shull 1983
double H2_coll_dissoc_cooling_coeff( const double& densGas, const double& densH2, const double& densH, const double& tempGas );


#endif
