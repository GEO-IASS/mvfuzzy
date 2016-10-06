// TODO:
//  - implement statistic indexes
//  - take constraints into account for all computations
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "util.h"
#include "matrix.h"
#include "model.h"

#define BUFF_SIZE 1024

bool verbose;
bool debug;
//st_matrix* load_dmtx(char **dmtx_file_name, int dmatrixc, int objc) {
//    st_matrix *dmatrix = malloc(sizeof(st_matrix) * dmatrixc);
//    size_t j;
//    for(j = 0; j < dmatrixc; ++j) {
//        init_st_matrix(dmatrix[j], objc, objc);
//        if(!load_data(dmtx_file_name[j], &dmatrix[j])) {
//            printf("[Warn] unexpected EOF on '%s'\n",
//                    dmtx_file_name[j]);
//        }
//    }
//    return dmatrix;
//}

int main(int argc, char **argv) {
    debug = true;
    verbose = true;
    FILE *cfgfile = fopen(argv[1], "r");
    if(!cfgfile) {
        printf("[Error] Could not open config file.\n");
        return 1;
    }
    size_t i;
	size_t j;
    // reading config file start
    int objc;
    fscanf(cfgfile, "%d", &objc);
    if(objc <= 0) {
        printf("[Error] objc <= 0.\n");
        return 2;
    }
    // reading labels
    int classc;
    int labels[objc];
    fscanf(cfgfile, "%d", &classc);
    for(i = 0; i < objc; ++i) {
        fscanf(cfgfile, "%d", &labels[i]);
    }
    // reading labels end
    int dmatrixc;
    fscanf(cfgfile, "%d", &dmatrixc);
    if(dmatrixc <= 0) {
        printf("[Error] dmatrixc <= 0.\n");
        return 2;
    }
    char dmtx_file_name[dmatrixc][BUFF_SIZE];
    for(j = 0; j < dmatrixc; ++j) {
        fscanf(cfgfile, "%s", dmtx_file_name[j]);
    }
    char out_file_name[BUFF_SIZE];
    fscanf(cfgfile, "%s", out_file_name);
    int clustc;
    fscanf(cfgfile, "%d", &clustc);
    if(clustc <= 0) {
        printf("[Error] clustc <= 0.\n");
        return 2;
    }
    int insts;
    fscanf(cfgfile, "%d", &insts);
    if(insts <= 0) {
        printf("[Error] insts <= 0.\n");
        return 2;
    }
    int max_iter;
    fscanf(cfgfile, "%d", &max_iter);
    double mfuz;
    fscanf(cfgfile, "%lf", &mfuz);
    if(!dgt(mfuz, 0.0)) {
        printf("[Error] mfuz <= 0.\n");
        return 2;
    }
    double epsilon;
    fscanf(cfgfile, "%lf", &epsilon);
    if(dlt(epsilon, 0.0)) {
        printf("[Error] epsilon < 0.\n");
        return 2;
    }
    int medoids_card;
    fscanf(cfgfile, "%d", &medoids_card);
    if(medoids_card <= 0) {
        printf("[Error] medoids_card <= 0.\n");
        return 2;
    }
    double theta;
    fscanf(cfgfile, "%lf", &theta);
    if(dlt(theta, 0.0)) {
        printf("[Error] theta < 0.\n");
        return 2;
    }
    double sample_perc;
    fscanf(cfgfile, "%lf", &sample_perc);
    if(dlt(sample_perc, 0.0)) {
        printf("[Error] sample_perc < 0.\n");
        return 2;
    }
    double alpha;
    fscanf(cfgfile, "%lf", &alpha);
    if(dlt(alpha, 0.0)) {
        printf("[Error] alpha < 0.\n");
        return 2;
    }
    fclose(cfgfile);
    // reading config file end 
    freopen(out_file_name, "w", stdout);
	double mfuzval = 1.0 / (mfuz - 1.0);
    printf("######Config summary:######\n");
    printf("Number of objects: %d\n", objc);
    printf("Number of clusters: %d\n", clustc);
    printf("Number of instances: %d\n", insts);
    printf("Number of iterations: %d\n", max_iter);
    printf("Parameter m: %.15lf\n", mfuz);
    printf("Epsilon: %.15lf\n", epsilon);
    printf("Medoids cardinality: %d\n", medoids_card);
    printf("Theta: %.15lf\n", theta);
    printf("Sample perc: %.15f\n", sample_perc);
    printf("Alpha: %.15lf\n", alpha);
    printf("###########################\n");
    // allocating memory start
    st_matrix *dmatrix = malloc(sizeof(st_matrix) * dmatrixc);
    for(j = 0; j < dmatrixc; ++j) {
        init_st_matrix(&dmatrix[j], objc, objc);
        if(!load_data(dmtx_file_name[j], &dmatrix[j])) {
            printf("[Warn] unexpected EOF on '%s'\n",
                    dmtx_file_name[j]);
        }
    }
    // allocating memory end
//    st_matrix *dmatrix = load_dmtx(dmtx_file_name, dmatrixc, objc);
//    char buf[BUFF_SIZE];
//    for(j = 0; j < dmatrixc; ++j) {
//        sprintf(buf, "out-%d.txt", j);
//        freopen(buf, "w", stdout);
//        print_st_matrix(&dmatrix[j], 4, false);
//    }
    srand(time(NULL));
//    srand(1234);
    size_t best_inst;
    double cur_adeq;
    double best_adeq;
    model_init(objc, clustc, dmatrixc, medoids_card, labels, classc,
            sample_perc);
    for(i = 1; i <= insts; ++i) {
        printf("Instance %d:\n", i);
        cur_adeq = run(dmatrix, max_iter, epsilon, theta, mfuz, alpha);
        if(i == 1 || cur_adeq < best_adeq) {
            save_env();
            best_adeq = cur_adeq;
            best_inst = i;
        }
        printf("\n");
    }
    printf("Process finished.\n");
    printf("Best adequacy: %.15lf (on instance %d)\n", best_adeq,
            best_inst);
    print_env();
    model_free();
    if(dmatrix) {
        for(j = 0; j < dmatrixc; ++j) {
            free_st_matrix(&dmatrix[j]);
        }
        free(dmatrix);
    }
    return 0;
}
