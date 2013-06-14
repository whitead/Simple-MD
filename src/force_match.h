#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "force.h"
#include "force_match_util.h"

void fm_loop(Run_Params* params);

double loss_function(double* forces_1, double* forces_2, unsigned int n);


void gradient(double* forces_ref, double* forces_target, double* delta_sigma, double* grad, double epsilon,unsigned int n);
