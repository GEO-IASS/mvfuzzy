#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "util.h"
#include "smedoids.h"
#include "lweight.h"
#include "model.h"
#include "stex.h"

void model_init(int objc, int clustc, int dmatrixc,
        int medoids_card, int *labels, int classc) {
    init_st_matrix(&weights, clustc, dmatrixc);
    init_st_matrix(&sav_weights, clustc, dmatrixc);
    init_st_matrix(&memb, objc, clustc);
    init_st_matrix(&sav_memb, objc, clustc);
    medoids = malloc(sizeof(size_t **) * clustc);
    sav_medoids = malloc(sizeof(size_t **) * clustc);
    size_t j;
    size_t k;
    for(k = 0; k < clustc; ++k) {
        medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
        sav_medoids[k] = malloc(sizeof(size_t *) * dmatrixc);
        for(j = 0; j < dmatrixc; ++j) {
            medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
            sav_medoids[k][j] = malloc(sizeof(size_t) * medoids_card);
        }
    }
    medoids_ncol = medoids_card;
    groups = asgroups(labels, objc, classc);
    if(debug) {
        print_groups(groups);
    }
}

void update_memb(st_matrix *dmatrix, double mfuzval) {
    size_t e;
    size_t h;
    size_t i;
    size_t j;
    size_t k;
    double vals[memb.ncol];
    double dsum;
    double new_memb;
    int zeroc;
    for(i = 0; i < memb.nrow; ++i) {
        zeroc = 0;
        for(k = 0; k < memb.ncol; ++k) {
            vals[k] = 0.0;
            for(j = 0; j < weights.ncol; ++j) {
                dsum = 0.0;
                for(e = 0; e < medoids_ncol; ++e) {
                    dsum += get(&dmatrix[j], i, medoids[k][j][e]);
                }
                vals[k] += get(&weights, k, j) * dsum;
            }
            if(deq(vals[k], 0.0)) {
                ++zeroc;
            }
        }
        if(zeroc) {
            new_memb = 1.0 / zeroc;
            for(k = 0; k < memb.ncol; ++k) {
                if(deq(vals[k], 0.0)) {
                    set(&memb, i, k, new_memb);
                } else {
                    set(&memb, i, k, 0.0);
                }
            }
        } else {
            for(k = 0; k < memb.ncol; ++k) {
                new_memb = 0.0;
                for(h = 0; h < memb.ncol; ++h) {
                    new_memb += pow(vals[k] / vals[h], mfuzval);
                }
                set(&memb, i, k, 1.0 / new_memb);
            }
        }
    }
}

void print_memb(st_matrix *memb) {
	print_header("Membership", HEADER_SIZE);
	size_t i;
	size_t k;
	double sum;
	for(i = 0; i < memb->nrow; ++i) {
		printf("%u: ", i);
		sum = 0.0;
		for(k = 0; k < memb->ncol; ++k) {
			printf("%lf ", get(memb, i, k));
            sum += get(memb, i, k);
		}
		printf("[%lf]", sum);
		if(!deq(sum, 1.0)) {
			printf("*\n");
		} else {
			printf("\n");
		}
	}
}

double adequacy(st_matrix *dmatrix, double mfuz) {
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double adeq = 0.0;
    double dsum;
    double wsum;
    for(k = 0; k < memb.ncol; ++k) {
        for(i = 0; i < memb.nrow; ++i) {
            wsum = 0.0;
            for(j = 0; j < weights.ncol; ++j) {
                dsum = 0.0;
                for(e = 0; e < medoids_ncol; ++e) {
                    dsum += get(&dmatrix[j], i, medoids[k][j][e]);
                }
                wsum += get(&weights, k, j) * dsum;
            }
            adeq += pow(get(&memb, i, k), mfuz) * wsum;
        }
    }
    return adeq;
}

void save_env() {
    mtxcpy(&sav_weights, &weights);
    mtxcpy(&sav_memb, &memb);
    size_t k;
    for(k = 0; k < memb.ncol; ++k) {
        mtxcpy_size_t(sav_medoids[k], medoids[k], weights.ncol,
                medoids_ncol);
    }
}

void print_env() {
    print_weights(&weights);
    print_medoids(medoids, memb.ncol, weights.ncol, medoids_ncol);
    print_memb(&memb);
}

double run(st_matrix *dmatrix, int max_iter, double epsilon,
        double theta, double mfuz, double sample_perc) {
    double mfuzval = 1.0 / (mfuz - 1.0);
    int objc = memb.nrow;
    int clustc = memb.ncol;
    int medoids_card = medoids_ncol;
    int dmatrixc = weights.ncol;
    size_t classc = groups->nrow;

    // generating constraints
    int **sample = gen_sample(groups, sample_perc);
    print_header(" Sample ", HEADER_SIZE);
    print_groups_(sample, classc);
    constraint **constr = gen_constraints(sample, classc, objc);
    print_constraints(constr, objc);
    size_t k;
    for(k = 0; k < classc; ++k) {
        free(sample[k]);
    }
    free(sample);

    init_weights(&weights);
    print_weights(&weights);
    init_medoids(medoids, objc, clustc, dmatrixc, medoids_card);
    print_medoids(medoids, clustc, dmatrixc, medoids_card);
    update_memb(dmatrix, mfuzval);
    print_memb(&memb);
    double prev_adeq = adequacy(dmatrix, mfuz);
    printf("\nAdequacy: %.15lf\n", prev_adeq);
    double cur_adeq;
    double prev_step_adeq = prev_adeq; // for debug
    double cur_step_adeq; // for debug
    double adeq_diff;
    size_t it;
    for(it = 0; it < max_iter; ++it) {
        printf("\nIteration %d\n", it);
        update_medoids_lw(medoids, medoids_card, &memb, dmatrixc,
                    dmatrix, mfuz);
        if(verbose) {
            print_medoids(medoids, clustc, dmatrixc, medoids_card);
        }
        if(debug) {
            cur_step_adeq = adequacy(dmatrix, mfuz);
            adeq_diff = prev_step_adeq - cur_step_adeq;
            if(adeq_diff < 0.0) {
                printf("[Warn] current step adequacy is greater than "
                        "previous (%.15lf)\n", - adeq_diff);
            }
            prev_step_adeq = cur_step_adeq;
        }
        update_weights_smd(&weights, &memb, medoids, medoids_card,
                dmatrix, theta, mfuz);
        if(verbose) {
            print_weights(&weights);
        }
        if(debug) {
            cur_step_adeq = adequacy(dmatrix, mfuz);
            adeq_diff = prev_step_adeq - cur_step_adeq;
            if(adeq_diff < 0.0) {
                printf("[Warn] current step adequacy is greater than "
                        "previous (%.15lf)\n", - adeq_diff);
            }
            prev_step_adeq = cur_step_adeq;
        }
        update_memb(dmatrix, mfuzval);
        if(verbose) {
            print_memb(&memb);
        }
        cur_adeq = adequacy(dmatrix, mfuz);
        adeq_diff = prev_adeq - cur_adeq;
        printf("\nAdequacy: %.15lf (%.15lf)\n", cur_adeq, adeq_diff);
        if(debug) {
            if(adeq_diff < 0.0) {
                printf("[Warn] current iteration adequacy is greater "
                        "than previous (%.15lf)\n", - adeq_diff);
            }
            adeq_diff = prev_step_adeq - cur_adeq; // cur_step_adeq
            if(adeq_diff < 0.0) {
                printf("[Warn] current step adequacy is greater than "
                        "previous (%.15lf)\n", - adeq_diff);
            }
            prev_step_adeq = cur_adeq;
        }
        if(adeq_diff < epsilon) {
            printf("Adequacy coefficient difference reached in %d "
                    "iterations.\n", it);
            break;
        }
        prev_adeq = cur_adeq;
    }
    size_t i;
    if(constr) {
        for(i = 0; i < objc; ++i) {
            if(constr[i]) {
                constraint_free(constr[i]);
            }
        }
        free(constr);
    }
    printf("\nClustering process finished.\n");
    print_weights(&weights);
    print_medoids(medoids, clustc, dmatrixc, medoids_card);
    print_memb(&memb);
    return cur_adeq;
}

void model_free() {
    if(medoids) {
        size_t j;
        size_t k;
        for(k = 0; k < memb.ncol; ++k) {
            if(medoids[k]) {
                for(j = 0; j < weights.ncol; ++j) {
                    if(medoids[k][j]) {
                        free(medoids[k][j]);
                    }
                }
                free(medoids[k]);
            }
        }
        free(medoids);
    }
    if(sav_medoids) {
        size_t j;
        size_t k;
        for(k = 0; k < memb.ncol; ++k) {
            if(sav_medoids[k]) {
                for(j = 0; j < weights.ncol; ++j) {
                    if(sav_medoids[k][j]) {
                        free(sav_medoids[k][j]);
                    }
                }
                free(sav_medoids[k]);
            }
        }
        free(sav_medoids);
    }
    free_st_matrix(&weights);
    free_st_matrix(&sav_weights);
    free_st_matrix(&memb);
    free_st_matrix(&sav_memb);
    if(groups) {
        free_st_matrix(groups);
        free(groups);
    }
}

