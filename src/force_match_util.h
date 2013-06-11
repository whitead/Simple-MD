/*
 * Load and stores the binary position of each trajectory frame. The frame_number
 * should be the maximum frame number and will be set to the actual after the call completes
 */
long* load_traj_file_mapping(FILE* trajectory, unsigned int* frame_number, unsigned int n_particles, unsigned int n_dims);

/*
 * Load the file_mapping for a neighbor list file. It says where to find the neighbor list
 * for each trajectory frame
 *
 */
long* load_nlist_file_mapping(FILE* nlist, unsigned int max_frames, unsigned int** nlist_length);

/*
 * Load a matrix in a file 
 */
void load_matrix_part(FILE* mfile, double* matrix, unsigned int nrow, unsigned int ncol);
