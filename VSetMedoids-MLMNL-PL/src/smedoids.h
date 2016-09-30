#ifndef _SMEDOIDS_H_
#define _SMEDOIDS_H_

#pragma once

#include "matrix.h"

void init_medoids(size_t ***medoids, int objc, int clustc,
                    int dmatrixc, int medoids_card);

void print_medoids(size_t ***medoids, int clustc, int dmatrixc,
        int medoids_card);


void update_medoids_lw(size_t ***medoids, int medoids_card,
                    st_matrix *memb, int dmatrixc,
                    st_matrix *dmatrix, double mfuz);

#endif /* _SMEDOIDS_H_ */
