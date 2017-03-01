#include "thermostat.h"
#include <gsl/gsl_randist.h>
#include <math.h>

double thermostat(double temperature, double time_step, void* thermostat_parameters, double* positions, double* velocities, double* forces, double* masses, unsigned int n_dims, unsigned int n_particles) {

  White_Params* params = (White_Params*) thermostat_parameters;

  unsigned int i, j, ndeg;
  double etemp, kenergy, dgdt;
  kenergy = 0;


  //calculate kinetic energy and get 
#pragma omp parallel for default(shared)	\
  private(etemp) reduction(+:kenergy)
  for(i = 0; i < n_particles; i++) {
    etemp = 0;
    for(j = 0; j < n_dims; j++) {
      etemp += velocities[i * n_dims + j] * velocities[i * n_dims + j];      
    }
    kenergy += 0.5 * etemp * masses[i];
  }

  //calculate  dgdt
#pragma omp parallel for default(shared)	\
  private(etemp) reduction(+:dgdt)
  for(i = 0; i < n_particles; i++) {
    etemp = 0;
    for(j = 0; j < n_dims; j++) {
      etemp += velocities[i * n_dims + j] / masses[i] * forces[i * n_dims + j] / masses[i];      
    }
    dgdt += etemp;
  }  
  
  //Removed COM is reason for -n_dims
  ndeg = (n_particles * n_dims - n_dims);

  params->sdot = -(kenergy - 0.5 * temperature * ndeg - params->s) * time_step;
  //params->sdot = (kenergy - 0.5 * temperature * ndeg) * time_step;
  //params->sdot = params->sdot * 9. / 10 + (kenergy - 0.5 * temperature * ndeg) / 10;
  //params->s += (params->sdot + dgdt - params->s) * time_step / params->mass;
  params->s += (params->sdot) * time_step / params->mass;
  //params->s += (params->sdot - params->s) * time_step / params->mass;
  //params->s += (params->sdot - params->s) * time_step / params->mass;


  //update forces
#pragma omp parallel for
  for(i = 0; i < n_particles; i++) {
    for(j = 0; j < n_dims; j++) {
      velocities[i * n_dims + j] *= (1 + params->s / ndeg);
    }
  }


  return (1 + params->s / ndeg);
}

void* build_white(unsigned int seed, double mass) {
  gsl_rng * rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(rng, seed);
  
  White_Params* params = (White_Params*) malloc(sizeof(White_Params));
  params->mass = mass;
  params->rng = rng;
  params->s = 0;
  params->sdot = 0;
  return((void*) params);

}

void free_thermostat(void* thermostat_parameters) {
 
  White_Params* params = (White_Params*) thermostat_parameters;
  free(params);
 
}
