#include "force_match.h"

int main(int argc, char* argv[]) {

  char* pfile = NULL;

  if(argc == 1) {
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

  double*  positions = (double*) malloc(sizeof(double) * params->n_particles * params->n_dims); 
  unsigned int frame_number = (params->steps) / params->position_log_period;
  long* traj_file_mapping = load_traj_file_mapping(traj_file, &frame_number, params->n_particles, params->n_dims);
  
  //get nlist file offsets
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

  double* ref_forces = (double*) malloc(sizeof(double) * params->n_particles);
  double* target_forces = (double*) malloc(sizeof(double) * params->n_particles);

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
    gather_forces(params->search_parameters, positions, target_forces, params->masses, params->box_size, params->n_dims, params->n_particles);
    
    //calculate gradient    
    
    //update 
    
  }

}
