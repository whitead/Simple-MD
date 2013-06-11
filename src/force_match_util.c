#include "force_match_util.h"

void load_matrix_part(FILE* mfile, double* matrix, unsigned int nrow, unsigned int ncol) {

  if(mfile != NULL) {

    unsigned int i, j;

    for(i = 0; i < nrow; i++) {
      for(j = 0; j < ncol; j++) {
	if(fscanf(mfile, "%lg", &(matrix[i*ncol + j])) == 0) {
	  fprintf(stderr, "Incorrect number of rows or columns"
		  "at i = %d, and j=%d, nrow=%d, ncol=%d\n", i, j, nrow, ncol);
	  exit(1);
	}
      }  
    }
    
  }else {
    perror("Could not open file\n");
  }

}

void load_int_matrix_part(FILE* mfile, unsigned int* matrix, unsigned int nrow, unsigned int ncol) {

  if(mfile != NULL) {

    unsigned int i, j;

    for(i = 0; i < nrow; i++) {
      for(j = 0; j < ncol; j++) {
	if(fscanf(mfile, "%ud", &(matrix[i*ncol + j])) == 0) {
	  fprintf(stderr, "Incorrect number of rows or columns"
		  "at i = %d, and j=%d, nrow=%d, ncol=%d\n", i, j, nrow, ncol);
	  exit(1);
	}
      }  
    }
    
  }else {
    perror("Could not open file\n");
  }

}

long* load_traj_file_mapping(FILE* trajectory, unsigned int* frame_number, unsigned int n_particles, unsigned int n_dims) {

  long* mapping = (long*) malloc(sizeof(long) * *frame_number);
  unsigned int i, j;
  double temp;

  //read all frames and store binary positions of the beginning of the frames
  for(i = 0; i < max_frames; i++) {
    mapping[i] = ftell(trajectory);
    for(j = 0; j < n_particles * n_dims; j++) {
      if(fscanf(trajectory, "%lg", &temp) == 0) {
	if(feof(trajectory)) {
#ifdef DEBUG
	  for(j = 0; j < max_frames; j++) {
	    printf("%d: %ld\n", j, mapping[j]);
	  }
#endif
	  realloc(mapping, sizeof(long) * i);
	  return mapping;
	}
      }
    }
  }

  return mapping;
}


long* load_nlist_file_mapping(FILE* nlist, unsigned int max_frames, unsigned int n_particles, unsigned int** nlist_length) {

  long* mapping = (long*) malloc(sizeof(long) * max_frames);
  *nlist_length = (unsigned int*) malloc(sizeof(unsigned int) * max_frames);
  unsigned int i, j;
  int temp, time_stamp;;

  //Get the binary position of the neighbor lists and get the neighbor list lenghts
  for(i = 0; i < max_frames; i++) {
    mapping[i] = ftell(nlist);
    //read in the length of the neighbor list and when the next frame update is
    if(fscanf(nlist, "%d ", &(*nlist_length)[i]) == 0 ||
       fscanf(nlist, "%d\n", &time_stamp) == 0) {
      if(feof(nlist)) {
#ifdef DEBUG
	for(j = 0; j < i; j++)
	  printf("%d\n", (*nlist_length)[i]);
#endif
	realloc(mapping, sizeof(long) * i);
	realloc(*nlist_length, sizeof(int) * i);
	return mapping;
      }
    } else {
      perror("Could not read nlist file");
    }

    //skip over the mapping until the next update
    for(j = i + 1; j < next_update; j++) {
      mapping[j] = ftell(nlist);
      *(nlist_length)[j] = *(nlist_lenght)[i];
#ifdef DEBUG
      printf("Filling frame %d nlist with data from frame %d\n", j, i);
#endif
    }

        
    for(j = 0; j < n_particles + (*nlist_length)[i]; j++) {
      if(fscanf(nlist, "%d", &temp) == 0) {
	if(feof(nlist)) {
	  fprintf(stderr, "Error: incomplete neighbor list in file\n");
	  exit(1);
	}
      }
    }
  }

  return mapping;
}
