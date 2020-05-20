#include "Main.h"
#include "configfile.h"

class AMUSE_SimpleX : public SimpleX {

  friend class SimpleX;
  
 public:
  
 AMUSE_SimpleX(string output_path, string data_path): SimpleX(output_path, data_path) {
    total_time=0.;
    syncflag=-1;
  };
  
  ~AMUSE_SimpleX(){};
  
  int add_vertex(long *id, double x,double y,double z,double rho,
		 double flux,double xion,double uInt,double xmol,double metallicity,double fuv_flux,double xCO,double turbulence);
  int add_site(long *id, double x,double y,double z,double rho,
	       double flux,double xion,double uInt,double xmol,double metallicity,double fuv_flux,double xCO,double turbulence);
  int remove_site(long id);
  int setup_simplex();

  void convert_units();
   
  int initialize(double current_time);
  int reinitialize();
  int evolve(double dt, int sync);
  int get_time(double *t){*t=total_time;return 0;}
  int set_time(double t){total_time=t;return 0;}
  int get_site(int id, double *x,double *y,double *z,double *rho,
	       double *flux,double *xion, double *uInt, double *xmol, double *metallicity, double *fuv_flux, double *xCO, double *turbulence);
  int get_position(int id, double *x, double *y, double *z);
  int get_density(int id, double *rho);
  int get_flux(int id, double *flux);
  int get_fuv_flux(int id, double *fuv_flux);
  int get_mean_intensity(int id, double *mean_intensity);
  int get_fuv_intensity(int id, double *fuv_intensity);
  int get_diffuse_intensity(int id, double *diffuse_intensity);
  int get_effective_surface(int id, double *surface);
  int get_ionisation(int id, double *xion);
  int get_molecular_fraction(int id, double *xmol);
  int get_CO_fraction(int id, double *xCO);
  int get_metallicity(int id, double *metallicity);
  int get_turbulence(int id, double *turbulence);
  int get_internalEnergy(int id, double *uInt);
  int get_dinternalEnergydt(int id, double *dudt);
  int set_site(int id, double x, double y, double z, double rho,
	       double flux, double xion, double uInt, double xmol, double metallicity, double fuv_flux, double xCO, double turbulence);
  int set_position(int id, double x, double y, double z);
  int set_density(int id, double rho);
  int set_flux(int id, double flux);
  int set_fuv_flux(int id, double fuv_flux);
  int set_ionisation(int id, double xion);
  int set_molecular_fraction(int id, double xmol);
  int set_CO_fraction(int id, double xCO);
  int set_metallicity(int id, double metallicity);
  int set_turbulence(int id, double turbulence);
  int set_internalEnergy(int id, double uInt);
  int set_dinternalEnergydt(int id, double dudt);
  int cleanup();
  int get_syncflag(){return syncflag;}
  int set_UNIT_T(double ts){ UNIT_T_MYR=ts;return 0;}
  int get_UNIT_T(double *ts){ *ts=UNIT_T_MYR;return 0;}
  int set_sizeBox(double bs){ sizeBox=bs;return 0;}
  int get_sizeBox(double *bs){ *bs=sizeBox;return 0;}
  int set_cosmic_ray_ionization_rate(double cr_ir){ cosmic_ray_ionization_rate=cr_ir;return 0;}
  int get_cosmic_ray_ionization_rate(double *cr_ir){ *cr_ir=cosmic_ray_ionization_rate;return 0;}
  int set_interstellar_fuv_field(double isrf){ interstellar_fuv_field=isrf;return 0;}
  int get_interstellar_fuv_field(double *isrf){ *isrf=interstellar_fuv_field;return 0;}
  int set_hilbert_order(int ho){ hilbert_order=ho;return 0;}
  int get_hilbert_order(int *ho){ *ho=hilbert_order;return 0;}

  int set_sourceTeff(double ts){ sourceTeff=ts;return 0;}
  int get_sourceTeff(double *ts){ *ts=sourceTeff;return 0;}
  int set_numFreq(int ts){ numFreq=ts;return 0;}
  int get_numFreq(int *ts){ *ts=numFreq;return 0;}
  int set_heat_cool(int ts){ heat_cool=ts;return 0;}
  int get_heat_cool(int *ts){ *ts=heat_cool;return 0;}
  int set_metal_cooling(int ts){ metal_cooling=ts;return 0;}
  int get_metal_cooling(int *ts){ *ts=metal_cooling;return 0;}
  int set_recombination_radiation(int ts){ rec_rad=ts;return 0;}
  int get_recombination_radiation(int *ts){ *ts=rec_rad;return 0;}
  int set_coll_ion(int ts){ coll_ion=ts;return 0;}
  int get_coll_ion(int *ts){ *ts=heat_cool;return 0;}
  int set_blackBody(int ts){ blackBody=ts;return 0;}
  int get_blackBody(int *ts){ *ts=blackBody;return 0;}
  int set_dust(int ts){ dust = ts;return 0;}
  int get_dust(int *ts){ *ts = dust;return 0;}
  int set_molecular(int ts){ molecular = ts;return 0;}
  int get_molecular(int *ts){ *ts = molecular;return 0;}
  int set_fuvprop(int ts){ fuvprop = ts;return 0;}
  int get_fuvprop(int *ts){ *ts = fuvprop;return 0;}
  int set_carb_monox(int ts){ 
    carbmonox = ts; 

    if (carbmonox) {
      xCO_max = 2.*( 3.31e-4 > 6.76e-4 ? 6.76e-4 : 3.31e-4 )/(3.31e-4 + 6.76e-4);

      rho_frac_H = (1. - 3.31e-4 - 6.76e-4)/(1. + 3.31e-4*11. + 6.76e-4*15.);
      rho_frac_C = 12.*3.31e-4/(1. + 3.31e-4*11. + 6.76e-4*15.);
      rho_frac_O = 16.*6.76e-4/(1. + 3.31e-4*11. + 6.76e-4*15.);
    }
    else {
      xCO_max = 0.;

      rho_frac_H = 1.;
      rho_frac_C = 0.;
      rho_frac_O = 0.;
    }

    return 0;}
  int get_carb_monox(int *ts){ *ts = carbmonox;return 0;}
 
  
  int set_numBorderSites(int ts){ borderSites=ts;return 0;}
  int get_numBorderSites(int *ts){ *ts=borderSites;return 0;}


 private:

  void read_parameters();
  double total_time;     
  int syncflag; 
};
