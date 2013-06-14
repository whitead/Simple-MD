#include "force.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

inline double lj(double r,  double epsilon, double sigma) {
  return 4 * epsilon * (6 * pow(sigma / r, 7) - 12 * pow(sigma / r, 13));  
}

//force derivative wrt sigma
inline double dslj(double r,  double epsilon, double sigma) {
  return 4 * epsilon * (42 * pow(sigma / r, 6) - 156 * pow(sigma / r, 12));  
}


inline double lj_trunc_shift(double r, double epsilon, double sigma, double rc, double shift) {
  if(r >= rc) {
    return 0;
  }
  return lj(r, epsilon, sigma) - shift;
}

//force derivative wrt sigma
inline double dslj_trunc_shift(double r, double epsilon, double sigma, double rc, double shift) {
  if(r >= rc) {
    return 0;
  }
  return dslj(r, epsilon, sigma) - shift;
}


double gather_forces(void* parameters, double* positions, double* forces, double* masses, double* box_size, unsigned int n_dims, unsigned int n_particles) {
  return gather_forces_and_deriv(parameters, positions, forces, NULL, masses, box_size, n_dims, n_particles);
}


double gather_forces_and_deriv(void* parameters, 
			       double* positions, 
			       double* forces, 
			       double* force_derivatives, 
			       double* masses,   
			       double* box_size, 
			       unsigned int n_dims, 
			       unsigned int n_particles) {

  const double epsilon = ((Lj_Parameters*) parameters)->epsilon;
  const double sigma = ((Lj_Parameters*) parameters)->sigma;
  Nlist_Parameters* nlist = ((Lj_Parameters*) parameters)->nlist;

  //update neighbor list
  update_nlist(positions, box_size, n_dims, n_particles, nlist);

  unsigned int i, j, k, n;
  int offset;
  double penergy = 0;
  double r, force, diff, dsforce;
  double force_vector[n_dims];
  double rcut = sqrt(nlist->rcut);
  double lj_shift = lj(rcut, epsilon, sigma);  
  double dslj_shift = dslj(rcut, epsilon, sigma);
  

  //zero forces
  #pragma omp parallel for
  for(i = 0; i < n_particles; i++)
    for(k = 0; k < n_dims; k++)
      forces[i * n_dims + k] = 0;

  //zero derivatives forces
  if(force_derivatives) {
#pragma omp parallel for
    for(i = 0; i < n_particles; i++)
      for(k = 0; k < n_dims; k++)
	force_derivatives[i * n_dims + k] = 0;
  }

  //iterate through all particles

  //This seems strange at first,
  //but it really just distributes
  //the work across the threads
  //The only trick is that the neighborlists
  //are not easy to spread out across the workers,
  //hence the conditionals.

#pragma omp parallel default(shared) \
  private(offset, n, i, j, k, r, force, dsforce, force_vector, diff)	\
  reduction(+:penergy)
  {

#ifdef _OPENMP
    offset = -1; //indicator
#else
    offset = 0; //zeroed
#endif   

#pragma omp for
    for(i = 0; i < n_particles; i++) {

#ifdef _OPENMP
      //accumulate the offset now that we know i
      if(offset == -1) {
	offset = 0;
	for(j = 0; j < i; j++) {
	  offset += nlist->nlist_count[j];
	}
      }
#endif

      //iterate through neighbor list
      for(n = offset; n - offset < nlist->nlist_count[i]; n++) {
	j = nlist->nlist[n];
	r = 0;
	
	//distance between particles
	for(k = 0; k < n_dims; k++) {
	  diff = min_image_dist(positions[j * n_dims + k] - positions[i * n_dims + k], box_size[k]);
	  r += diff * diff;
	  force_vector[k] = diff;
	}
	
	r = sqrt(r);
	//LJ force and potential
	force = lj_trunc_shift(r, epsilon, sigma, rcut, lj_shift);

	//calculate sigma derivative
	dsforce = lj_trunc_shift(r, epsilon, sigma, rcut, dslj_shift);
	
	
#ifdef DEBUG
	printf("F(%d - %d, %g) = %g\n", i, j, r, force);
	printf("dF/dSigma(%d - %d, %g) = %g\n", i, j, r, dsforces);
#endif //DEBUG
	
#pragma omp critical (update_forces)
	for(k = 0; k < n_dims; k++)  {
	  forces[i * n_dims + k] += force / r * force_vector[k];
	  forces[j * n_dims + k] -= force / r * force_vector[k];
	  if(force_derivatives) {
	    force_derivatives[i * n_dims + k] += dsforce / r * force_vector[k];
	    force_derivatives[j * n_dims + k] -= dsforce / r * force_vector[k];
	  }
	}
	
	penergy += 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6)) - 4 * epsilon * (pow(sigma / rcut, 12) - pow(sigma / rcut, 6));
      }

      offset += nlist->nlist_count[i];
    }  

  }

  return(penergy);  
}

void* build_lj(double epsilon, double sigma, Nlist_Parameters* nlist) {
  Lj_Parameters init = {.epsilon = epsilon, .sigma = sigma};
  Lj_Parameters* parameters = (Lj_Parameters*) malloc(sizeof(Lj_Parameters));
  memcpy(parameters, &init, sizeof(init));
  parameters->nlist = nlist;
  return (parameters);
}

void free_forces(void* parameters) {

  Lj_Parameters* lj_parameters = (Lj_Parameters*) parameters;
  Nlist_Parameters* nlist = ((Lj_Parameters*) lj_parameters)->nlist;
  free_nlist(nlist);
  free(lj_parameters);

}
