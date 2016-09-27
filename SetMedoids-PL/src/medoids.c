#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "util.h"
#include "medoids.h"

void init_medoids(size_t **medoids, int objc, int clustc,
                    int medoids_card) {
	size_t i;
	size_t e;
	size_t k;
	int obj;
	bool chosen[objc];
	for(k = 0; k < clustc; ++k) {
		for(i = 0; i < objc; ++i) {
			chosen[i] = false;
		}
		for(e = 0; e < medoids_card; ++e) {
			do {
				obj = rand() % objc;
			} while(chosen[obj]);
			medoids[k][e] = obj;
			chosen[obj] = true;
		}
	}
}

void print_medoids(size_t **medoids, int clustc, int medoids_card) {
    print_header("Medoids", HEADER_SIZE);
    size_t e;
    size_t k;
    for(k = 0; k < clustc; ++k) {
        for(e = 0; e < medoids_card; ++e) {
            printf("%d ", medoids[k][e]);
        }
        printf("\n");
    }
}

typedef struct objnval {
	size_t obj;
	double val;
} objnval;

static int objnval_cmp(const void *p1, const void *p2) {
	const objnval *a = (const objnval *) p1;
	const objnval *b = (const objnval *) p2;
	return (a->val > b->val) - (a->val < b->val);
}

void update_medoids_lw(size_t **medoids, int medoids_card,
                    st_matrix *memb, st_matrix *lweights,
                    st_matrix *dmatrix, double mfuz) {
	size_t h;
	size_t i;
	size_t j;
    size_t k;
    size_t objc = memb->nrow;
    size_t clustc = memb->ncol;
    objnval candidates[objc];
    double val;
    double dsum;
    for(k = 0; k < clustc; ++k) {
        for(h = 0; h < objc; ++h) {
            candidates[h].obj = h;
            val = 0.0;
            for(i = 0; i < objc; ++i) {
                dsum = 0.0;
                for(j = 0; j < lweights->ncol; ++j) {
                    dsum += get(lweights, k, j) *
                                get(&dmatrix[j], i, h);
                }
                val += pow(get(memb, i, k), mfuz) * dsum;
            }
            candidates[h].val = val;
        }
        qsort(candidates, objc, sizeof(objnval), objnval_cmp);
        for(h = 0; h < medoids_card; ++h) {
            medoids[k][h] = candidates[h].obj;
        }
    }
}

//void update_medoids() {
//	size_t h;
//	size_t i;
//	size_t j;
//	size_t k;
//	objnval candidates[objc];
//	double sumweights;
//	for(k = 0; k < clustc; ++k) {
//		for(h = 0; h < objc; ++h) {
//			candidates[h].obj = h;
//			candidates[h].val = 0.0;
//			for(i = 0; i < objc; ++i) {
//				sumweights = 0.0;
//				for(j = 0; j < dmatrixc; ++j) {
//					sumweights += weights[k][j] * dmatrix[j][i][h];
//				}
//				candidates[h].val += pow(memb[i][k], mfuz) * sumweights;
//			}
//		}
//		qsort(candidates, objc, sizeof(objnval), objnval_cmp);
//        for(h = 0; h < medoids_card; ++h) {
//            medoids[k][h] = candidates[h].obj;
//        }
//	}
//}
