#include "force_match_util.h"

void load_matrix_part(FILE* mfile, double* matrix, unsigned int nrow, unsigned int ncol) {

  if(mfile != NULL) {

    unsigned int i, j;

    for(i = 0; i < nrow; i++) {
      for(j = 0; j < ncol; j++) {
	if(fscanf(mfile, "%lg", &(matrix[i*ncol + j])) == EOF) {
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
	if(fscanf(mfile, "%d", &(matrix[i*ncol + j])) == EOF) {
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
  for(i = 0; i < *frame_number; i++) {
    mapping[i] = ftell(trajectory);
    for(j = 0; j < n_particles * n_dims; j++) {
      if(fscanf(trajectory, "%lg", &temp) == EOF) {
	if(feof(trajectory)) {
#ifdef DEBUG
	  for(j = 0; j < *frame_number; j++) {
	    printf("%d: %ld\n", j, mapping[j]);
	  }
#endif
	  mapping = (long *) realloc(mapping, sizeof(long) * i);
	  return mapping;
	}
      }
    }
  }

  return mapping;
}


long* load_nlist_file_mapping(FILE* nlist, unsigned int frame_number, unsigned int n_particles, unsigned int** nlist_length) {

  long* mapping = (long*) malloc(sizeof(long) * frame_number);
  *nlist_length = (unsigned int*) malloc(sizeof(unsigned int) * frame_number);
  unsigned int i, j, time_stamp, temp;

  //Get the binary position of the neighbor lists and get the neighbor list lenghts
  for(i = 0; i < frame_number; i++) {

    //read in the length of the neighbor list and when this neighbor list was created
    if(fscanf(nlist, "%u ", &temp) == EOF) {
      if(feof(nlist)) { //all done reading file
	break;
      } else {
	perror("Could not read nlist file");
	exit(1);
      }
    }

    if(fscanf(nlist, "%u\n", &time_stamp) == EOF) {
      perror("Could not read nlist file");
      exit(1);
    }

    //check to see if we're ending early
    if(time_stamp >= frame_number) {
      //reduce index, so that when we backfill after the loop it uses the previous neighbor list
      i--;
      break;
    }

    //store the list, now that we know when it was created
    mapping[time_stamp] = ftell(nlist);
    (*nlist_length)[time_stamp] = temp;

      

    //backfill mapping up until the time stamp
#ifdef DEBUG
    printf("time_stamp = %u, length = %u\n\n", time_stamp, temp);
#endif

    for(j = i; j < time_stamp && j < frame_number; j++) {
      mapping[j] = mapping[i - 1];
      (*nlist_length)[j] = (*nlist_length)[i - 1];
#ifdef DEBUG
      printf("Filling frame %d nlist with data from frame %d\n", j, i - 1);
#endif
    }
    
    i = time_stamp;
        
    for(j = 0; j < n_particles + (*nlist_length)[i]; j++) {
      if(fscanf(nlist, "%d", &temp) == EOF) {
	if(feof(nlist)) {
	  fprintf(stderr, "Error: incomplete neighbor list in file\n");
	  exit(1);
	}
      }
    }
  }

  //if we didn't reach the last frame, fill with last known neighbor list
  for(j = i; j < frame_number; j++) {
    mapping[j] = mapping[i - 1];
    (*nlist_length)[j] = (*nlist_length)[i - 1];
#ifdef DEBUG
    printf("Filling frame %d nlist with data from frame %d\n", j, i - 1);
#endif
  }


#ifdef DEBUG
  for(j = 0; j < frame_number; j++)
    printf("mapping[%u] = %ld\n", j, mapping[j]);
#endif

  return mapping;
}
