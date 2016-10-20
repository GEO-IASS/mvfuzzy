#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "util.h"
#include "smedoids.h"
#include "lweight.h"
#include "model.h"
#include "stex.h"

void model_init(int objc, int clustc, int dmatrixc,
        int medoids_card, int *labels, int classc,
        double sample_perc) {
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

    // generating constraints
    st_matrix *groups = asgroups(labels, objc, classc);
    if(debug) {
        print_groups(groups);
    }
    int **sample = gen_sample(groups, sample_perc);
    print_header(" Sample ", HEADER_SIZE);
    print_groups_(sample, classc);
    free_st_matrix(groups);
    free(groups);
    constr = gen_constraints(sample, classc, objc);
    print_constraints(constr, objc);
    for(k = 0; k < classc; ++k) {
        free(sample[k]);
    }
    free(sample);
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
            new_memb = 1.0 / (double) zeroc;
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

void update_memb_constr_(st_matrix *dmatrix, double alpha) {
    size_t c;
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    bool set_v[memb.ncol];
    double mtx_a[memb.ncol];
    double mtx_b[memb.ncol];
    double sum_num;
    double sum_den;
    double val;
    int obj;
    double gamma;
    bool changed;
    for(i = 0; i < memb.nrow; ++i) {
        for(k = 0; k < memb.ncol; ++k) {
            mtx_a[k] = 0.0;
            for(j = 0; j < weights.ncol; ++j) {
                val = 0.0;
                for(e = 0; e < medoids_ncol; ++e) {
                    val += get(&dmatrix[j], i, medoids[k][j][e]);
                }
                mtx_a[k] += get(&weights, k, j) * val;
            }
            mtx_a[k] *= 2.0;
        }
        for(k = 0; k < memb.ncol; ++k) {
            set_v[k] = true;
        }
    }
    do {
        changed = false;
        for(i = 0; i < memb.nrow; ++i) {
            sum_num = 1.0;
            sum_den = 0.0;
            for(k = 0; k < memb.ncol; ++k) {
                val = 0.0;
                if(constr[i]) {
                    for(e = 0; e < constr[i]->ml->size; ++e) {
                        obj = constr[i]->ml->get[e];
                        for(c = 0; c < memb.ncol; ++c) {
                            if(c != k) {
                                val += get(&memb, obj, c);
                            }
                        }
                    }
                    for(e = 0; e < constr[i]->mnl->size; ++e) {
                        obj = constr[i]->mnl->get[e];
                        val += get(&memb, obj, k);
                    }
                    val *= alpha;
                }
                mtx_b[k] = val;
                if(set_v[k]) {
                    sum_num += mtx_b[k] / mtx_a[k];
                    sum_den += 1.0 / mtx_a[k];
                }
            }
            gamma = sum_num / sum_den;
            for(k = 0; k < memb.ncol; ++k) {
                val = (gamma - mtx_b[k]) / mtx_a[k];
                if(dgt(val, 0.0)) {
                    changed = !deq(val, get(&memb, i, k));
                    set(&memb, i, k, val);
                } else {
                    set(&memb, i, k, 0.0);
                    set_v[k] = false;
                    changed = true;
                }
            }
        }
    } while(changed);
}

void update_memb_constr(st_matrix *dmatrix, double alpha) {
    size_t c;
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    double mtx_a[memb.nrow][memb.ncol];
    bool set_a[memb.nrow][memb.ncol];
    size_t set_a_c[memb.nrow];
    double val;
    double dsum;
    for(i = 0; i < memb.nrow; ++i) {
        set_a_c[i] = 0;
        for(k = 0; k < memb.ncol; ++k) {
            val = 0.0;
            for(j = 0; j < weights.ncol; ++j) {
                dsum = 0.0;
                for(e = 0; e < medoids_ncol; ++e) {
                    dsum += get(&dmatrix[j], i, medoids[k][j][e]);
                }
                val += get(&weights, k, j) * dsum;
            }
            mtx_a[i][k] = val;
            if(deq(val, 0.0)) {
                set_a[i][k] = true;
                set_a_c[i] += 1;
            } else {
                set_a[i][k] = false;
            }
        }
    }
    bool set_v[memb.nrow][memb.ncol];
    double mtx_b[memb.nrow][memb.ncol];
    int obj;
    for(i = 0; i < memb.nrow; ++i) {
        if(!set_a_c[i]) {
            for(k = 0; k < memb.ncol; ++k) {
                set_v[i][k] = true;
                mtx_a[i][k] *= 2.0;
                val = 0.0;
                if(constr[i]) {
                    for(e = 0; e < constr[i]->ml->size; ++e) {
                        obj = constr[i]->ml->get[e];
                        for(c = 0; c < memb.ncol; ++c) {
                            if(c != k) {
                                val += get(&memb, obj, c);
                            }
                        }
                    }
                    for(e = 0; e < constr[i]->mnl->size; ++e) {
                        obj = constr[i]->mnl->get[e];
                        val += get(&memb, obj, k);
                    }
//                    val *= 2.0;
                }
                mtx_b[i][k] = alpha * val;
            }
        }
    }
    double sum_num;
    double sum_den;
    double gamma;
    bool test;
    do {
        test = false;
        for(i = 0; i < memb.nrow; ++i) {
            if(!set_a_c[i]) {
                sum_num = 0.0;
                sum_den = 0.0;
                for(k = 0; k < memb.ncol; ++k) {
                    if(set_v[i][k]) {
                        sum_num += mtx_b[i][k] / mtx_a[i][k];
                        sum_den += 1.0 / mtx_a[i][k];
                    }
                }
                gamma = (1.0 + sum_num) / sum_den;
                for(k = 0; k < memb.ncol; ++k) {
                    if(set_v[i][k]) {
                        val = (gamma - mtx_b[i][k]) / mtx_a[i][k];
                        if(dgt(val, 0.0)) {
                            set(&memb, i, k, val);
                        } else {
                            set(&memb, i, k, 0.0);
                            set_v[i][k] = false;
                            test = true;
                        }
                    }
                }
            }
        }
    } while(test);
    for(i = 0; i < memb.nrow; ++i) {
        if(set_a_c[i]) {
            val = 1.0 / (double) set_a_c[i];
            for(k = 0; k < memb.ncol; ++k) {
                if(set_a[i][k]) {
                    set(&memb, i, k, val);
                } else {
                    set(&memb, i, k, 0.0);
                }
            }
        }
    }
}

void update_memb_constr_old(st_matrix *dmatrix, double alpha) {
    size_t c;
    size_t e;
    size_t i;
    size_t j;
    size_t k;
    // indicates whether set_a[i][k] is inside A
    bool set_a[memb.nrow][memb.ncol];
    // cardinality of set_a[i]
    size_t set_a_c[memb.nrow];
    double mtx_a[memb.nrow][memb.ncol];
    double val;
    double dsum;
    for(i = 0; i < memb.nrow; ++i) {
        set_a_c[i] = 0;
        for(k = 0; k < memb.ncol; ++k) {
            val = 0.0;
            for(j = 0; j < weights.ncol; ++j) {
                dsum = 0.0;
                for(e = 0; e < medoids_ncol; ++e) {
                    dsum += get(&dmatrix[j], i, medoids[k][j][e]);
                }
                val += get(&weights, k, j) * dsum;
            }
            mtx_a[i][k] = val;
            if(deq(val, 0.0)) {
                set_a[i][k] = true;
                set_a_c[i] += 1;
            } else {
                set_a[i][k] = false;
            }
        }
    }
    // indicates whether set_v[i][k] is inside V
    bool set_v[memb.nrow][memb.ncol];
    double mtx_b[memb.nrow][memb.ncol];
    int obj;
    for(i = 0; i < memb.nrow; ++i) {
        if(!set_a_c[i]) {
            for(k = 0; k < memb.ncol; ++k) {
                set_v[i][k] = true;
                mtx_a[i][k] *= 2.0;
                val = 0.0;
                if(constr[i]) {
                    for(e = 0; e < constr[i]->ml->size; ++e) {
                        obj = constr[i]->ml->get[e];
                        for(c = 0; c < memb.ncol; ++c) {
                            if(c != k) {
                                val += get(&memb, obj, c);
                            }
                        }
                    }
                    for(e = 0; e < constr[i]->mnl->size; ++e) {
                        obj = constr[i]->mnl->get[e];
                        val += get(&memb, obj, k);
                    }
//                    val *= 2.0;
                }
                mtx_b[i][k] = alpha * val;
            }
        }
    }
    double sum_num;
    double sum_den;
    double gamma;
    bool test;
    do {
        test = false;
        for(i = 0; i < memb.nrow; ++i) {
            if(!set_a_c[i]) {
                sum_num = 0.0;
                sum_den = 0.0;
                for(k = 0; k < memb.ncol; ++k) {
                    if(set_v[i][k]) {
                        sum_num += mtx_b[i][k] / mtx_a[i][k];
                        sum_den += 1.0 / mtx_a[i][k];
                    }
                }
                gamma = (1.0 + sum_num) / sum_den;
                for(k = 0; k < memb.ncol; ++k) {
                    if(set_v[i][k]) {
                        val = (gamma - mtx_b[i][k]) / mtx_a[i][k];
                        if(!dgt(val, 0.0)) {
                            set(&memb, i, k, 0.0);
                            set_v[i][k] = false;
                            test = true;
                        } else {
                            set(&memb, i, k, val);
                        }
                    }
                }
            }
        }
    } while(test);
    for(i = 0; i < memb.nrow; ++i) {
        if(set_a_c[i]) {
            val = 1.0 / set_a_c[i];
            for(k = 0; k < memb.ncol; ++k) {
                if(set_a[i][k]) {
                    set(&memb, i, k, val);
                } else {
                    set(&memb, i, k, 0.0);
                }
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

double constr_adequacy(st_matrix *dmatrix, double alpha) {
    size_t c;
    size_t e;
    size_t i;
    size_t k;
    int obj;
    double adeq = 0.0;
    for(i = 0; i < memb.nrow; ++i) {
        if(constr[i]) {
            for(e = 0; e < constr[i]->ml->size; ++e) {
                obj = constr[i]->ml->get[e];
                for(c = 0; c < memb.ncol; ++c) {
                    for(k = 0; k < memb.ncol; ++k) {
                        if(c != k) {
                            adeq += get(&memb, i, c) *
                                        get(&memb, obj, k);
                        }
                    }
                }
            }
            for(e = 0; e < constr[i]->mnl->size; ++e) {
                obj = constr[i]->mnl->get[e];
                for(c = 0; c < memb.ncol; ++c) {
                    adeq += get(&memb, i, c) * get(&memb, obj, c);
                }
            }
        }
    }
    return adeq;
}

double prev_adeq_unconstr;
double prev_adeq_constr;

double adequacy(st_matrix *dmatrix, double mfuz, double alpha) {
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
    double cadeq = constr_adequacy(dmatrix, alpha);
    if(debug) {
        printf("[Debug]Adequacy: %.15lf %.15lf\n", adeq, cadeq);
        if(prev_adeq_unconstr && dgt(adeq, prev_adeq_unconstr)) {
            printf("[Debug] current unconstrained adequacy is greater"
                    " than previous.\n");
        }
        if(prev_adeq_constr && dgt(cadeq, prev_adeq_constr)) {
            printf("[Debug] current constrained adequacy is greater "
                    "than previous.\n");
        }
        prev_adeq_unconstr = adeq;
        prev_adeq_constr = cadeq;
    }
    return adeq + cadeq;
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
        double theta, double mfuz, double alpha) {
    double mfuzval = 1.0 / (mfuz - 1.0);
    int objc = memb.nrow;
    int clustc = memb.ncol;
    int medoids_card = medoids_ncol;
    int dmatrixc = weights.ncol;
    prev_adeq_unconstr = 0.0;
    prev_adeq_constr = 0.0;

    init_weights(&weights);
    print_weights(&weights);
    init_medoids(medoids, objc, clustc, dmatrixc, medoids_card);
    print_medoids(medoids, clustc, dmatrixc, medoids_card);
//    update_memb_constr(dmatrix, mfuzval, alpha);
    update_memb(dmatrix, mfuzval);
    print_memb(&memb);
    double prev_adeq = adequacy(dmatrix, mfuz, alpha);
    printf("\nAdequacy: %.15lf\n", prev_adeq);
    double cur_adeq = 0.0;
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
            cur_step_adeq = adequacy(dmatrix, mfuz, alpha);
            // no need to check for first it in this model
            if(it) {
                adeq_diff = prev_step_adeq - cur_step_adeq;
                if(dlt(adeq_diff, 0.0)) {
                    printf("[Warn] current step adequacy is greater than "
                            "previous (%.15lf)\n", - adeq_diff);
                }
            }
            prev_step_adeq = cur_step_adeq;
        }
        update_weights_smd(&weights, &memb, medoids, medoids_card,
                dmatrix, theta, mfuz);
        if(verbose) {
            print_weights(&weights);
        }
        if(debug) {
            cur_step_adeq = adequacy(dmatrix, mfuz, alpha);
            adeq_diff = prev_step_adeq - cur_step_adeq;
            if(dlt(adeq_diff, 0.0)) {
                printf("[Warn] current step adequacy is greater than "
                        "previous (%.15lf)\n", - adeq_diff);
            }
            prev_step_adeq = cur_step_adeq;
        }
        update_memb_constr(dmatrix, alpha);
//        update_memb(dmatrix, mfuzval);
        if(verbose) {
            print_memb(&memb);
        }
        cur_adeq = adequacy(dmatrix, mfuz, alpha);
        if(debug) {
            cur_step_adeq = cur_adeq;
            adeq_diff = prev_step_adeq - cur_step_adeq;
            if(dlt(adeq_diff, 0.0)) {
                printf("[Warn] current step adequacy is greater than "
                        "previous (%.15lf)\n", - adeq_diff);
            }
            prev_step_adeq = cur_step_adeq;
        }
        adeq_diff = prev_adeq - cur_adeq;
        printf("\nAdequacy: %.15lf (%.15lf)\n", cur_adeq, adeq_diff);
        // no need to check for first it in this model
        if(debug && it) {
            if(dlt(adeq_diff, 0.0)) {
                printf("[Warn] current iteration adequacy is greater "
                        "than previous (%.15lf)\n", - adeq_diff);
            }
        }
        if(fabs(adeq_diff) < epsilon) {
            printf("Adequacy coefficient difference reached in %d "
                    "iterations.\n", it);
            break;
        }
        prev_adeq = cur_adeq;
    }
    printf("\nClustering process finished.\n");
    print_weights(&weights);
    print_medoids(medoids, clustc, dmatrixc, medoids_card);
    print_memb(&memb);
    return cur_adeq;
}

void model_free() {
    size_t i;
    size_t j;
    size_t k;
    if(medoids) {
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
    if(constr) {
        for(i = 0; i < memb.nrow; ++i) {
            if(constr[i]) {
                constraint_free(constr[i]);
            }
        }
        free(constr);
    }
    free_st_matrix(&weights);
    free_st_matrix(&sav_weights);
    free_st_matrix(&memb);
    free_st_matrix(&sav_memb);
}

