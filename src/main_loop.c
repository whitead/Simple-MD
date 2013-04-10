#include "main_loop.h"

int main(int argc, char* argv[]) {

  //process arguemnts to get parameters
  if(argc != 2) {
    printf("Usage: [run file]\n");
    exit(1);
  }

  FILE* arguments = fopen(argv[1], "r");
  if(arguments == NULL) {
    perror("Could not open file\n");
    exit(1);
  }

  Run_Params* p = read_parameters(arguments, NULL);

  //start main loop
  main_loop(p);

  free_run_params(p);

  return(0);

}

int main_loop(Run_Params* params){

  unsigned int i;
  double* positions = params->initial_positions;
  double* velocities = params->initial_velocities;
  double* forces = malloc(sizeof(double) * params->n_dims * params->n_particles);
  double penergy = 0;
  double kenergy = 0;
  double insta_temperature;

  gather_forces(params->force_parameters, positions, forces, params->masses, params->n_dims, params->n_particles);


  printf("%12s %12s %12s %12s %12s %12s\n", "Step", "Time", "T", "PE", "KE", "E");
  
  for(i = 0; i < params->steps; i++) {

    //integrate 1

    integrate_1(params->time_step, positions, velocities, forces, params->masses, params->n_dims, params->n_particles);


    //gather forces
    penergy =  gather_forces(params->force_parameters, positions, forces, params->masses, params->n_dims, params->n_particles);

    //integrate 2
    integrate_2(params->time_step, positions, velocities, forces, params->masses, params->n_dims, params->n_particles);

    //thermostat
    thermostat(params->temperature, params->time_step, params->thermostat_parameters, positions, velocities, params->masses, params->n_dims, params->n_particles);

    //calculate important quantities
    kenergy = calculate_kenergy(velocities, params->masses, params->n_dims, params->n_particles);
    insta_temperature = kenergy * 2 / (params->n_particles * params->n_dims);
    
    
    //print
    if(i % params->position_log_period == 0)
      log_array(params->positions_file, positions, params->n_dims, params->n_particles, false);

    if(i % params->velocity_log_period == 0)
      log_array(params->velocities_file, velocities, params->n_dims, params->n_particles, true);

    if(i % params->force_log_period == 0)
      log_array(params->forces_file, forces, params->n_dims, params->n_particles, true);

    if(i % params->print_period == 0) {
      printf("%12d %12g %12g %12g %12g %12g\n", i, i * params->time_step, insta_temperature, penergy, kenergy, penergy + kenergy);
    }
  }

  if(params->positions_file != NULL)
    fclose(params->positions_file);
  if(params->velocities_file != NULL)
    fclose(params->velocities_file);
  if(params->forces_file != NULL)
    fclose(params->forces_file);

  free(forces);
  
  return MD_SUCCESS;
}

