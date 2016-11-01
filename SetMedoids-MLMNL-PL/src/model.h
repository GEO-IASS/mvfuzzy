#ifndef _MODEL_H_
#define _MODEL_H_

#pragma once

//#include <stdio.h>
//#include <stdbool.h>

#include "matrix.h"
#include "constraint.h"

#define BUFF_SIZE 1024

extern bool debug;
extern bool verbose;

// Variables for 'save_env'
st_matrix sav_weights;
st_matrix sav_memb;
size_t **sav_medoids;

st_matrix weights;
st_matrix memb;
size_t **medoids;
size_t medoids_ncol;
// constraints for each object
constraint **constr;

// Allocates memory for weights, medoids and memb
void model_init(int objc, int clustc, int dmatrixc, int medoids_card,
        int *labels, int classc, double sample_perc);

// Saves current weights, medoids and memb
void save_env();

// Prints saved environment. At least one call to 'save_env' must
// have been made
void print_env();

// Main loop
double run(st_matrix *dmatrix, int max_iter, double epsilon,
        double theta, double mfuz, double alpha);

// Computes and prints various statistical indexes
void compute_idxs(st_matrix *dmatrix, int *labels, double mfuz);

// Frees weights, medoids and memb
void model_free();

void print_memb(st_matrix *memb);

#endif /* _MODEL_H_ */

