# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>

# include "global.h"
# include "rand.h"

int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;

int main(int argc, char **argv) {
    int i;
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *gp;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

    popsize = 4;
    ngen = 1000;
    nobj = 2;
    ncon = 2;
    nreal = 2;
    nbin = 0;
    pcross_real = 0.6;
    pmut_real = 1/nreal;
    eta_c = 5;
    eta_m = 5;

    choice = 1;

    float nrealLimites[2][2] = {{0,1},{0,1}};


    if (argc < 2) {
        printf("\n Usage ./nsga2r random_seed \n");
        exit(1);
    }
    seed = (double) atof(argv[1]);
    if (seed <= 0.0 || seed >= 1.0) {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    fpt1 = fopen("initial_pop.out", "w");
    fpt2 = fopen("final_pop.out", "w");
    fpt3 = fopen("best_pop.out", "w");
    fpt4 = fopen("all_pop.out", "w");
    fpt5 = fopen("params.out", "w");
    fprintf(fpt1, "# This file contains the data of initial population\n");
    fprintf(fpt2, "# This file contains the data of final population\n");
    fprintf(fpt3, "# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4, "# This file contains the data of all generations\n");
    fprintf(fpt5, "# This file contains information about inputs as read by the program\n");

    if (nreal != 0) {
        min_realvar = (double *) malloc(nreal * sizeof(double));
        max_realvar = (double *) malloc(nreal * sizeof(double));
        for (i = 0; i < nreal; i++) {
            min_realvar[i] = nrealLimites[i][0];
            max_realvar[i] = nrealLimites[i][1];
        }
    }

    if (choice > 0) {
        gp = popen(GNUPLOT_COMMAND, "w");
        if (gp == NULL) {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }
        if (nobj == 2) {
            obj1 = 1;
            obj2 = 2;
            obj3 = -1;
        } else {
            if (choice == 2) {
                obj1 = 1;
                obj2 = 2;
                obj3 = -1;
            } else {
                obj1 = 1;
                obj2 = 2;
                obj3 = 3;
                angle1 = 45;
                angle2 = 90;
            }
        }
    }

    printf("\n Input data successfully entered, now performing initialization \n");

    fprintf(fpt5, "\n Population size = %d", popsize);
    fprintf(fpt5, "\n Number of generations = %d", ngen);
    fprintf(fpt5, "\n Number of objective functions = %d", nobj);
    fprintf(fpt5, "\n Number of constraints = %d", ncon);
    fprintf(fpt5, "\n Number of real variables = %d", nreal);
    if (nreal != 0) {
        for (i = 0; i < nreal; i++) {
            fprintf(fpt5, "\n Lower limit of real variable %d = %e", i + 1, min_realvar[i]);
            fprintf(fpt5, "\n Upper limit of real variable %d = %e", i + 1, max_realvar[i]);
        }
        fprintf(fpt5, "\n Probability of crossover of real variable = %e", pcross_real);
        fprintf(fpt5, "\n Probability of mutation of real variable = %e", pmut_real);
        fprintf(fpt5, "\n Distribution index for crossover = %e", eta_c);
        fprintf(fpt5, "\n Distribution index for mutation = %e", eta_m);
    }
    fprintf(fpt5, "\n Number of binary variables = %d", nbin);
    if (nbin != 0) {
        for (i = 0; i < nbin; i++) {
            fprintf(fpt5, "\n Number of bits for binary variable %d = %d", i + 1, nbits[i]);
            fprintf(fpt5, "\n Lower limit of binary variable %d = %e", i + 1, min_binvar[i]);
            fprintf(fpt5, "\n Upper limit of binary variable %d = %e", i + 1, max_binvar[i]);
        }
        fprintf(fpt5, "\n Probability of crossover of binary variable = %e", pcross_bin);
        fprintf(fpt5, "\n Probability of mutation of binary variable = %e", pmut_bin);
    }
    fprintf(fpt5, "\n Seed for random number generator = %e", seed);
    bitlength = 0;
    if (nbin != 0) {
        for (i = 0; i < nbin; i++) {
            bitlength += nbits[i];
        }
    }
    fprintf(fpt1,
            "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",
            nobj, ncon, nreal, bitlength);
    fprintf(fpt2,
            "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",
            nobj, ncon, nreal, bitlength);
    fprintf(fpt3,
            "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",
            nobj, ncon, nreal, bitlength);
    fprintf(fpt4,
            "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",
            nobj, ncon, nreal, bitlength);
    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    parent_pop = (population *) malloc(sizeof(population));
    child_pop = (population *) malloc(sizeof(population));
    mixed_pop = (population *) malloc(sizeof(population));
    allocate_memory_pop(parent_pop, popsize);
    allocate_memory_pop(child_pop, popsize);
    allocate_memory_pop(mixed_pop, 2 * popsize);
    randomize();
    initialize_pop(parent_pop);
    printf("\n Initialization done, now performing first generation");
    decode_pop(parent_pop);
    evaluate_pop(parent_pop);
    assign_rank_and_crowding_distance(parent_pop);
    report_pop(parent_pop, fpt1);
    fprintf(fpt4, "# gen = 1\n");
    report_pop(parent_pop, fpt4);
    printf("\n gen = 1");
    printf("OPA?");
    fflush(stdout);
    if (choice != 0) onthefly_display(parent_pop, gp, 1);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    sleep(1);
    for (i = 2; i <= ngen; i++) {
        selection(parent_pop, child_pop);
        mutation_pop(child_pop);
        decode_pop(child_pop);
        evaluate_pop(child_pop);
        merge(parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort(mixed_pop, parent_pop);
        /* Comment following four lines if information for all
        generations is not desired, it will speed up the execution */
        fprintf(fpt4, "# gen = %d\n", i);
        report_pop(parent_pop, fpt4);
        fflush(fpt4);
        if (choice != 0) onthefly_display(parent_pop, gp, i);
        printf("\n gen = %d", i);
    }
    printf("\n Generations finished, now reporting solutions");
    report_pop(parent_pop, fpt2);
    report_feasible(parent_pop, fpt3);
    if (nreal != 0) {
        fprintf(fpt5, "\n Number of crossover of real variable = %d", nrealcross);
        fprintf(fpt5, "\n Number of mutation of real variable = %d", nrealmut);
    }
    if (nbin != 0) {
        fprintf(fpt5, "\n Number of crossover of binary variable = %d", nbincross);
        fprintf(fpt5, "\n Number of mutation of binary variable = %d", nbinmut);
    }
    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    fclose(fpt5);
    if (choice != 0) {
        pclose(gp);
    }
    if (nreal != 0) {
        free(min_realvar);
        free(max_realvar);
    }
    if (nbin != 0) {
        free(min_binvar);
        free(max_binvar);
        free(nbits);
    }
    deallocate_memory_pop(parent_pop, popsize);
    deallocate_memory_pop(child_pop, popsize);
    deallocate_memory_pop(mixed_pop, 2 * popsize);
    free(parent_pop);
    free(child_pop);
    free(mixed_pop);
    printf("\n Routine successfully exited \n");
    return (0);
}
