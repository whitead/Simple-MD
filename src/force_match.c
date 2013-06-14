#include "force_match.h"

int main(int argc, char* argv[]) {

  char* pfile = NULL;

  if(argc == 2) {
    pfile = argv[1];
  }

  Run_Params* p = read_parameters(pfile);

  fm_loop(p);

  return 0;
}


void fm_loop(Run_Params* params) {

  int i;

  FILE* traj_file = params->positions_file;
  if(!traj_file) {
    perror("Could not open trajectory file");
  }

  //get trajectory frame file offsets
  double*  positions = params->initial_positions;
  unsigned int frame_number = (params->steps) / params->position_log_period;
  long* traj_file_mapping = load_traj_file_mapping(traj_file, &frame_number, params->n_particles, params->n_dims);
  
  //get nlist file offsets and parameters
  Nlist_Parameters* np =( (Nlist_Parameters*) ( (Lj_Parameters*) params->force_parameters)->nlist);
  FILE* nlist_file = np->output_file;
  if(!nlist_file) {
    perror("Could not open nlist file");
  }

  
  //load nlist from file and set-up variables
  unsigned int* nlist_lengths;
  long* nlist_file_mapping = load_nlist_file_mapping(nlist_file, frame_number, params->n_particles, &nlist_lengths);
  unsigned int max_nlist_size = 0;
  for(i = 0; i < frame_number; i++)
    if(nlist_lengths[i] > max_nlist_size)
      max_nlist_size = nlist_lengths[i];

  np->nlist_count = (unsigned int*) malloc(sizeof(unsigned int) * params->n_particles);
  np->nlist = (unsigned int*) malloc(sizeof(unsigned int) * max_nlist_size);

  //turn off neighbor list rebuilding
  np->do_not_rebuild = 1;

  double* ref_forces = (double*) malloc(sizeof(double) * params->n_particles * params->n_dims);
  double* target_forces = (double*) malloc(sizeof(double) * params->n_particles * params->n_dims);
  double* force_deriv = (double*) malloc(sizeof(double) * params->n_particles * params->n_dims);
  double* grad = (double*) malloc(sizeof(double) * 2); 

  //optimization parameters
  double ball_radius = 2;
  double lipschitz[2] = {0, 0};
  Lj_Parameters* lj_search = params->search_parameters;
  unsigned int* frames = (unsigned int*) malloc(sizeof(unsigned int) * frame_number);
  unsigned int current_frame, sampled_index, grad_index;
  for(i = 0; i < frame_number; i++)
    frames[i] = i;

  gsl_rng * rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (gsl_rng_default);
  long seed = time(NULL); //pseudo random seed
  gsl_rng_set(rng, seed); 

  //conservative play parameters
  double loss;

  for(i = frame_number; i >= 0; i--) {

    //get random index
    if(i > 0) {
      sampled_index = gsl_rng_uniform_int(rng, i);    
      current_frame = frames[sampled_index];
      frames[sampled_index] = frames[i - 1];
    } else {
      current_frame = frames[i];
    }

#ifdef DEBUG
    printf("Sampled frame %d for force_matching\n", current_frame);
#endif

    //load trjacetory frame
    fseek(traj_file, traj_file_mapping[current_frame], SEEK_SET); //go to the frame in the file
    load_matrix_part(traj_file, positions, params->n_particles, params->n_dims); //load the frame

    //load the neighbor list
    fseek(nlist_file, nlist_file_mapping[current_frame], SEEK_SET);
    load_int_matrix_part(nlist_file, np->nlist_count, params->n_particles, 1); //load neighbor list counts
    load_int_matrix_part(nlist_file, np->nlist, nlist_lengths[current_frame], 1); //load neighbor list
    
    //gather forces according to parameters
    gather_forces(params->force_parameters, positions, ref_forces, params->masses, params->box_size, params->n_dims, params->n_particles);
    
    //gather forces according to parameters guess
    gather_forces_and_deriv(lj_search, positions, target_forces, force_deriv, params->masses, params->box_size, params->n_dims, params->n_particles);

    //only update when necessary
    loss = loss_function(ref_forces, target_forces, params->n_particles * params->n_dims);
    if(loss == 0) {
      continue;
    }
    
    printf("Loss = %g\n", loss);
    
    //calculate gradient 
    gradient(ref_forces, target_forces, force_deriv, grad,
	     lj_search->epsilon, 
	     params->n_particles * params->n_dims);

    printf("dF / dS = %g, dF / dE = %g\n", grad[0], grad[1]);
        
    //choose which component to update
    grad_index = gsl_rng_uniform_int(rng, 2);
    //update gradient accumulator
    if(i % 25 == 0) {
      lipschitz[0] = 0;
      lipschitz[1] = 0;     
    }
    lipschitz[grad_index] += grad[grad_index] * grad[grad_index];

    printf("sigma = %g, epsilon = %g update %s\n", lj_search->sigma, lj_search->epsilon, grad_index ? "epsilon" : "sigma");
    if(grad_index == 0)
      lj_search->sigma = lj_search->sigma - (ball_radius * sqrt(2))  / sqrt(lipschitz[0]) * grad[0];
    else
      lj_search->epsilon = lj_search->epsilon - (ball_radius * sqrt(2))  / sqrt(lipschitz[1]) * grad[1];

    //project back onto search space
    if(lj_search->sigma < 0)
      lj_search->sigma = 0.01;
    if(lj_search->epsilon < 0)
      lj_search->epsilon = 0.01;

    printf("sigma' = %g, epsilon' = %g, eta = %g\n", lj_search->sigma, lj_search->epsilon, (ball_radius * sqrt(2))  / sqrt(lipschitz[grad_index]));
  }

  free(grad);
  free(ref_forces);
  free(target_forces);
  free(force_deriv);
  free(frames);
}


double loss_function(double* forces_1, double* forces_2, unsigned int n) {
  unsigned int i;
  double result = 0;
  for(i = 0; i < n; i++)
    result += pow(forces_1[i] - forces_2[i], 2);
  
  return result;
  

}

void gradient(double* ref_forces, double* target_forces, double* delta_sigma, double* grad, double epsilon, unsigned int n) {
  unsigned int i;
  grad[0] = grad[1] = 0;
  for(i = 0; i < n; i++) {
    grad[0] -= 2 * (ref_forces[i] - target_forces[i]) * delta_sigma[i];
    grad[1] -= 2 * (ref_forces[i] - target_forces[i]) * target_forces[i] / epsilon;
  }
}

