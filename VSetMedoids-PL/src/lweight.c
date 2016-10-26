#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "lweight.h"

void print_weights(st_matrix *weights) {
    print_header("Weights", HEADER_SIZE);
    size_t j;
    size_t k;
    double prod;
    double val;
    for(k = 0; k < weights->nrow; ++k) {
        prod = 1.0;
        for(j = 0; j < weights->ncol; ++j) {
            val = get(weights, k, j);
            if(dlt(val, 0.0)) {
                printf("!");
            }
			printf("%lf ", val);
			prod *= val;
        }
        printf("[%lf]", prod);
        if(!deq(prod,1.0)) {
            printf(" =/= 1.0?");
        }
        printf("\n");
    }
}

void init_weights(st_matrix *weights) {
    setall(weights, 1.0);
}

void update_weights_md(st_matrix *weights, st_matrix *memb, 
        size_t **medoids, size_t medoids_card, st_matrix *dmatrix,
        double theta, double mfuz) {
    size_t objc = memb->nrow;
    size_t dmatrixc = weights->ncol;
    size_t clustc = weights->nrow;
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    size_t belowc;
    double v_vals[dmatrixc];
    double chi;
    double above;
    double val;
    double dsum;
    for(k = 0; k < clustc; ++k) {
        belowc = 0;
        for(j = 0; j < dmatrixc; ++j) {
            val = 0.0;
            for(i = 0; i < objc; ++i) {
                dsum = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    dsum += get(&dmatrix[j], i, medoids[k][e]);
                }
                val += pow(get(memb, i, k), mfuz) * dsum;
            }
            v_vals[j] = val;
            if(!dgt(val, theta)) {
                ++belowc;
            }
        }
        chi = 1.0;
        above = 1.0;
        val = 1.0 / ((double) (dmatrixc - belowc));
        for(j = 0; j < dmatrixc; ++j) {
            if(dgt(v_vals[j], theta)) {
                above *= pow(v_vals[j], val);
            } else {
                chi *= pow(get(weights, k, j), val);
            }
        }
        val = (1.0 / chi) * above;
        for(j = 0; j < dmatrixc; ++j) {
            if(v_vals[j] > theta) {
                set(weights, k, j, val / v_vals[j]);
            }
        }
    }
}

void update_weights_smd(st_matrix *weights, st_matrix *memb, 
        size_t ***medoids, size_t medoids_card, st_matrix *dmatrix,
        double theta, double mfuz) {
    size_t objc = memb->nrow;
    size_t dmatrixc = weights->ncol;
    size_t clustc = weights->nrow;
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    size_t belowc;
    double v_vals[dmatrixc];
    double chi;
    double above;
    double val;
    double dsum;
    for(k = 0; k < clustc; ++k) {
        belowc = 0;
        chi = 1.0;
        above = 1.0;
        for(j = 0; j < dmatrixc; ++j) {
            val = 0.0;
            for(i = 0; i < objc; ++i) {
                dsum = 0.0;
                for(e = 0; e < medoids_card; ++e) {
                    dsum += get(&dmatrix[j], i, medoids[k][j][e]);
                }
                val += pow(get(memb, i, k), mfuz) * dsum;
            }
            v_vals[j] = val;
            if(dgt(val, theta)) {
                above *= val;
            } else {
                chi *= get(weights, k, j);
                ++belowc;
            }
        }
        chi = 1.0 / chi;
        val = 1.0 / ((double) (dmatrixc - belowc));
        val = pow(chi * above, val);
        for(j = 0; j < dmatrixc; ++j) {
            if(v_vals[j] > theta) {
                set(weights, k, j, val / v_vals[j]);
            }
        }
    }
}
