#ifndef _MEDOIDS_H_
#define _MEDOIDS_H_

#pragma once

#include "matrix.h"

void init_medoids(size_t **medoids, int objc, int clustc,
                    int medoids_card);

void print_medoids(size_t **medoids, int clustc, int medoids_card);


void update_medoids_lw(size_t **medoids, int medoids_card,
                    st_matrix *memb, st_matrix *lweights,
                    st_matrix *dmatrix, double mfuz);

#endif /* _MEDOIDS_H_ */
