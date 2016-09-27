#ifndef _LWEIGHT_H_
#define _LWEIGHT_H_

#pragma once

#include "matrix.h"

void print_weights(st_matrix *weights);

void init_weights(st_matrix *weights);

void update_weights_md(st_matrix *weights, st_matrix *memb, 
        size_t **medoids, size_t medoids_card, st_matrix *dmatrix,
        double theta, double mfuz);

#endif /* _LWEIGHT_H_ */
