#ifndef _MODEL_H_
#define _MODEL_H_

#pragma once

//#include <stdio.h>
//#include <stdbool.h>

#include "matrix.h"

#define BUFF_SIZE 1024

extern bool debug;
extern bool verbose;

// Variables for 'save_env'
st_matrix sav_weights;
st_matrix sav_memb;
size_t ***sav_medoids;

st_matrix weights;
st_matrix memb;
size_t ***medoids;
size_t medoids_ncol;

// Allocates memory for weights, medoids and memb
void model_init(int objc, int clustc, int dmatrixc, int medoids_card);

// Saves current weights, medoids and memb
void save_env();

// Prints saved environment. At least one call to 'save_env' must
// have been made
void print_env();

// Main loop
double run(st_matrix *dmatrix, int max_iter, double epsilon,
        double theta, double mfuz);

// Frees weights, medoids and memb
void model_free();

void print_memb(st_matrix *memb);

//typedef struct progpar {
//    int objc;
//    int classc;
//    int *labels;
//    int dmatrixc;
//    char **dmtx_file_name;
//    char *out_file_name;
//    int clustc;
//    int insts;
//    int max_iter;
//    double mfuz;
//    double epsilon;
//    int medoids_card;
//    double theta;
//} progpar;

//bool read_par(FILE *cfgfile, int *objc, int *classc, int *labels,
//                int *dmatrixc, char **dmtx_file_name,
//                char *out_file_name, int *clustc, int *insts,
//                int *max_iter, double *mfuz, double *epsilon,
//                int *medoids_card, double *theta);

#endif /* _MODEL_H_ */
