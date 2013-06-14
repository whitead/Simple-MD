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

  unsigned int i;

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
  

  for(i = 0; i < frame_number; i++) {

    //load trjacetory frame
    fseek(traj_file, traj_file_mapping[i], SEEK_SET); //go to the frame in the file
    load_matrix_part(traj_file, positions, params->n_particles, params->n_dims); //load the frame

    //load the neighbor list
    fseek(nlist_file, nlist_file_mapping[i], SEEK_SET);
    load_int_matrix_part(nlist_file, np->nlist_count, params->n_particles, 1); //load neighbor list counts
    load_int_matrix_part(nlist_file, np->nlist, nlist_lengths[i], 1); //load neighbor list
    
    //gather forces according to parameters
    gather_forces(params->force_parameters, positions, ref_forces, params->masses, params->box_size, params->n_dims, params->n_particles);
    
    //gather forces according to parameters guess
    gather_forces_and_deriv(params->search_parameters, positions, target_forces, force_deriv, params->masses, params->box_size, params->n_dims, params->n_particles);

    //calculate delta F^2
    printf("loss = %g\n", loss_function(ref_forces, target_forces, params->n_particles * params->n_dims));
    
    //calculate gradient 
    gradient(ref_forces, target_forces, force_deriv, grad,
	     ((Lj_Parameters*)params->search_parameters)->epsilon, 
	     params->n_particles * params->n_dims);

    printf("dF / dS = %g, dF / dE = %g\n", grad[0], grad[1]);
    
    
    //update 
    
  }

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

