#include <iostream>
#include <cmath>
#include <unistd.h>

using namespace std;

# define INF 1.0e14
# define EPS 1.0e-14
# define E  2.71828182845905
# define PI 3.14159265358979
# define GNUPLOT_COMMAND "gnuplot -persist"

typedef struct {
    int rank;
    double constr_violation;
    double *xreal;
    int **gene;
    double *xbin;
    double *obj;
    double *constr;
    double crowd_dist;
} individual;

typedef struct {
    individual *ind;
} population;

typedef struct lists {
    int index;
    struct lists *parent;
    struct lists *child;
} list;

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

double seed;
double oldrand[55];
int jrand;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

const int tamanhoVetores = 25;
const float medias[25] = {22.07f,22.75f,13.02f,43.82f,50.08f,35.27f,20.43f,25.53f,9.57f,15.26f,38.86f,13.63f,22.82f,13.46f,15.53f,7.84f,11.86f,14.15f,29.20f,21.50f,34.80f,17.36f,29.31f,9.55f,7.27f};
const float desvio[25] = {2.11f,0.87f,0.80f,2.41f,2.41f,2.52f,1.78f,1.03f,0.45f,1.55f,3.06f,1.73f,2.66f,0.61f,1.13f,0.63f,0.60f,0.87f,6.06f,1.06f,1.48f,1.20f,1.31f,0.90f,0.57f};
const float correl[25][25] = {
        {1.000f,0.759f,0.819f,0.436f,0.723f,0.056f,0.990f,0.697f,-0.490f,-0.868f,0.623f,-0.916f,-0.674f,0.688f,0.818f,0.887f,-0.697f,0.850f,-0.876f,0.620f,-0.553f,0.776f,-0.449f,-0.579f,0.829f},
        {0.759f,1.000f,0.634f,0.107f,0.432f,-0.038f,0.736f,0.594f,-0.451f,-0.802f,0.509f,-0.826f,-0.692f,0.165f,0.546f,0.803f,-0.700f,0.674f,-0.841f,0.391f,-0.548f,0.721f,-0.566f,-0.774f,0.513f},
        {0.819f,0.634f,1.000f,0.197f,0.936f,0.461f,0.825f,0.884f,-0.121f,-0.670f,0.857f,-0.673f,-0.307f,0.499f,0.678f,0.833f,-0.354f,0.593f,-0.717f,0.485f,-0.264f,0.473f,-0.025f,-0.296f,0.712f},
        {0.436f,0.107f,0.197f,1.000f,0.150f,-0.475f,0.456f,0.046f,-0.511f,-0.316f,-0.120f,-0.408f,-0.569f,0.674f,0.697f,0.317f,-0.545f,0.674f,-0.339f,0.702f,-0.565f,0.486f,-0.458f,-0.078f,0.704f},
        {0.723f,0.432f,0.936f,0.150f,1.000f,0.548f,0.727f,0.832f,0.031f,-0.546f,0.828f,-0.520f,-0.129f,0.548f,0.566f,0.704f,-0.149f,0.489f,-0.555f,0.435f,-0.104f,0.320f,0.184f,-0.076f,0.619f},
        {0.056f,-0.038f,0.461f,-0.475f,0.548f,1.000f,0.104f,0.459f,0.627f,0.126f,0.703f,0.127f,0.606f,-0.074f,-0.088f,0.161f,0.497f,-0.322f,0.119f,-0.297f,0.574f,-0.364f,0.673f,0.341f,-0.036f},
        {0.990f,0.736f,0.825f,0.456f,0.727f,0.104f,1.000f,0.680f,-0.442f,-0.826f,0.645f,-0.899f,-0.628f,0.695f,0.825f,0.880f,-0.683f,0.847f,-0.842f,0.611f,-0.522f,0.739f,-0.424f,-0.546f,0.846f},
        {0.697f,0.594f,0.884f,0.046f,0.832f,0.459f,0.680f,1.000f,-0.024f,-0.545f,0.786f,-0.538f,-0.192f,0.375f,0.584f,0.707f,-0.212f,0.408f,-0.607f,0.321f,-0.107f,0.389f,0.013f,-0.244f,0.601f},
        {-0.490f,-0.451f,-0.121f,-0.511f,0.031f,0.627f,-0.442f,-0.024f,1.000f,0.693f,0.187f,0.619f,0.819f,-0.282f,-0.416f,-0.343f,0.789f,-0.621f,0.636f,-0.411f,0.786f,-0.687f,0.640f,0.582f,-0.353f},
        {-0.868f,-0.802f,-0.670f,-0.316f,-0.546f,0.126f,-0.826f,-0.545f,0.693f,1.000f,-0.447f,0.902f,0.787f,-0.433f,-0.629f,-0.791f,0.773f,-0.780f,0.913f,-0.522f,0.695f,-0.776f,0.494f,0.673f,-0.599f},
        {0.623f,0.509f,0.857f,-0.120f,0.828f,0.703f,0.645f,0.786f,0.187f,-0.447f,1.000f,-0.531f,-0.012f,0.244f,0.485f,0.735f,-0.123f,0.309f,-0.538f,0.233f,-0.009f,0.130f,0.117f,-0.215f,0.520f},
        {-0.916f,-0.826f,-0.673f,-0.408f,-0.520f,0.127f,-0.899f,-0.538f,0.619f,0.902f,-0.531f,1.000f,0.797f,-0.478f,-0.770f,-0.879f,0.852f,-0.854f,0.943f,-0.599f,0.726f,-0.787f,0.663f,0.763f,-0.753f},
        {-0.674f,-0.692f,-0.307f,-0.569f,-0.129f,0.606f,-0.628f,-0.192f,0.819f,0.787f,-0.012f,0.797f,1.000f,-0.410f,-0.613f,-0.621f,0.891f,-0.839f,0.796f,-0.673f,0.856f,-0.851f,0.833f,0.767f,-0.562f},
        {0.688f,0.165f,0.499f,0.674f,0.548f,-0.074f,0.695f,0.375f,-0.282f,-0.433f,0.244f,-0.478f,-0.410f,1.000f,0.671f,0.482f,-0.356f,0.685f,-0.422f,0.609f,-0.297f,0.527f,-0.213f,-0.108f,0.731f},
        {0.818f,0.546f,0.678f,0.697f,0.566f,-0.088f,0.825f,0.584f,-0.416f,-0.629f,0.485f,-0.770f,-0.613f,0.671f,1.000f,0.759f,-0.641f,0.780f,-0.704f,0.705f,-0.539f,0.597f,-0.486f,-0.358f,0.978f},
        {0.887f,0.803f,0.833f,0.317f,0.704f,0.161f,0.880f,0.707f,-0.343f,-0.791f,0.735f,-0.879f,-0.621f,0.482f,0.759f,1.000f,-0.621f,0.774f,-0.872f,0.634f,-0.539f,0.680f,-0.436f,-0.606f,0.758f},
        {-0.697f,-0.700f,-0.354f,-0.545f,-0.149f,0.497f,-0.683f,-0.212f,0.789f,0.773f,-0.123f,0.852f,0.891f,-0.356f,-0.641f,-0.621f,1.000f,-0.829f,0.806f,-0.561f,0.863f,-0.782f,0.834f,0.792f,-0.612f},
        {0.850f,0.674f,0.593f,0.674f,0.489f,-0.322f,0.847f,0.408f,-0.621f,-0.780f,0.309f,-0.854f,-0.839f,0.685f,0.780f,0.774f,-0.829f,1.000f,-0.815f,0.777f,-0.747f,0.823f,-0.657f,-0.633f,0.789f},
        {-0.876f,-0.841f,-0.717f,-0.339f,-0.555f,0.119f,-0.842f,-0.607f,0.636f,0.913f,-0.538f,0.943f,0.796f,-0.422f,-0.704f,-0.872f,0.806f,-0.815f,1.000f,-0.636f,0.742f,-0.751f,0.583f,0.781f,-0.686f},
        {0.620f,0.391f,0.485f,0.702f,0.435f,-0.297f,0.611f,0.321f,-0.411f,-0.522f,0.233f,-0.599f,-0.673f,0.609f,0.705f,0.634f,-0.561f,0.777f,-0.636f,1.000f,-0.673f,0.557f,-0.452f,-0.360f,0.719f},
        {-0.553f,-0.548f,-0.264f,-0.565f,-0.104f,0.574f,-0.522f,-0.107f,0.786f,0.695f,-0.009f,0.726f,0.856f,-0.297f,-0.539f,-0.539f,0.863f,-0.747f,0.742f,-0.673f,1.000f,-0.670f,0.730f,0.669f,-0.497f},
        {0.776f,0.721f,0.473f,0.486f,0.320f,-0.364f,0.739f,0.389f,-0.687f,-0.776f,0.130f,-0.787f,-0.851f,0.527f,0.597f,0.680f,-0.782f,0.823f,-0.751f,0.557f,-0.670f,1.000f,-0.684f,-0.689f,0.578f},
        {-0.449f,-0.566f,-0.025f,-0.458f,0.184f,0.673f,-0.424f,0.013f,0.640f,0.494f,0.117f,0.663f,0.833f,-0.213f,-0.486f,-0.436f,0.834f,-0.657f,0.583f,-0.452f,0.730f,-0.684f,1.000f,0.781f,-0.442f},
        {-0.579f,-0.774f,-0.296f,-0.078f,-0.076f,0.341f,-0.546f,-0.244f,0.582f,0.673f,-0.215f,0.763f,0.767f,-0.108f,-0.358f,-0.606f,0.792f,-0.633f,0.781f,-0.360f,0.669f,-0.689f,0.781f,1.000f,-0.320f},
        {0.829f,0.513f,0.712f,0.704f,0.619f,-0.036f,0.846f,0.601f,-0.353f,-0.599f,0.520f,-0.753f,-0.562f,0.731f,0.978f,0.758f,-0.612f,0.789f,-0.686f,0.719f,-0.497f,0.578f,-0.442f,-0.320f,1.000f}
};

/*const int tamanhoVetores = 3;
const float medias[3] = {11.75, 18.25, 9.75};
const float desvio[3] = {5.402f, 9.443f, 3.491f};
const float correl[3][3] = {{1, -0.748f, 0.367f}, {-0.748f, 1, -0.528f}, {0.367f, -0.528f, 1}};*/

double calcularRetorno(const double *participacoes) {
    float out = 0;
    for (int i = 0; i <= tamanhoVetores; i++) {
        out += participacoes[i] * medias[i];
    }
    return out;
}

double calcularRisco(const double *participacoes) {
    double out = 0;

    auto *WxS = new double[tamanhoVetores]();
    auto *vec = new double[tamanhoVetores]();

    for (int i = 0; i < tamanhoVetores; ++i) {
        WxS[i] = participacoes[i] * desvio[i];
    }

    for (int i = 0; i < tamanhoVetores; ++i) { // linha
        for (int j = 0; j < tamanhoVetores; ++j) { // coluna
            vec[i] += correl[i][j] * WxS[j];
        }
    }

    for (int i = 0; i < tamanhoVetores; ++i) {
        out += WxS[i] * vec[i];
    }

    out = sqrt(out);

    return out;
}

void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr){
    obj[0] = - calcularRetorno(xreal);
    obj[1] = calcularRisco(xreal);

    double sumXreal = 0;
    for (int i = 0; i < tamanhoVetores; ++i) {
        sumXreal+= xreal[i];
    }

    constr[0] = 1 - sumXreal;
    constr[1] = sumXreal - 0.9;
}

void inicializaTeste(){
    int i = 0;

    seed = 0.001f;
    popsize = 100;
    ngen = 1000;
    nobj = 2;
    ncon = 2;
    nreal = tamanhoVetores;
    pcross_real = 0.9; // 0.6 - 1.0
    eta_c = 5; // 5 - 20
    eta_m = 10; // 5 - 50
    nbin = 0;
    choice = 1;
    obj1 = 1;
    obj2 = 2;
    obj3 = -1;

    pmut_real = 1 / nreal;

    min_realvar = (double *) malloc(nreal * sizeof(double));
    max_realvar = (double *) malloc(nreal * sizeof(double));
    for (i = 0; i < nreal; i++) {
        min_realvar[i] = 0;
        max_realvar[i] = 1;
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void advance_random() {
    int j1;
    double new_random;
    for (j1 = 0; j1 < 24; j1++) {
        new_random = oldrand[j1] - oldrand[j1 + 31];
        if (new_random < 0.0) {
            new_random = new_random + 1.0;
        }
        oldrand[j1] = new_random;
    }
    for (j1 = 24; j1 < 55; j1++) {
        new_random = oldrand[j1] - oldrand[j1 - 24];
        if (new_random < 0.0) {
            new_random = new_random + 1.0;
        }
        oldrand[j1] = new_random;
    }
}

void warmup_random(double seed) {
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = seed;
    new_random = 0.000000001;
    prev_random = seed;
    for (j1 = 1; j1 <= 54; j1++) {
        ii = (21 * j1) % 54;
        oldrand[ii] = new_random;
        new_random = prev_random - new_random;
        if (new_random < 0.0) {
            new_random += 1.0;
        }
        prev_random = oldrand[ii];
    }
    advance_random();
    advance_random();
    advance_random();
    jrand = 0;
    return;
}

void randomize() {
    int j1;
    for (j1 = 0; j1 <= 54; j1++) {
        oldrand[j1] = 0.0;
    }
    jrand = 0;
    warmup_random(seed);
    return;
}

double randomperc() {
    jrand++;
    if (jrand >= 55) {
        jrand = 1;
        advance_random();
    }
    return ((double) oldrand[jrand]);
}

int rnd(int low, int high) {
    int res;
    if (low >= high) {
        res = low;
    } else {
        res = low + (randomperc() * (high - low + 1));
        if (res > high) {
            res = high;
        }
    }
    return (res);
}

double rndreal(double low, double high) {
    return (low + (high - low) * randomperc());
}

void allocate_memory_ind(individual *ind) {
    int j;
    if (nreal != 0) {
        ind->xreal = (double *) malloc(nreal * sizeof(double));
    }
    if (nbin != 0) {
        ind->xbin = (double *) malloc(nbin * sizeof(double));
        ind->gene = (int **) malloc(nbin * sizeof(int *));
        for (j = 0; j < nbin; j++) {
            ind->gene[j] = (int *) malloc(nbits[j] * sizeof(int));
        }
    }
    ind->obj = (double *) malloc(nobj * sizeof(double));
    if (ncon != 0) {
        ind->constr = (double *) malloc(ncon * sizeof(double));
    }
}

void allocate_memory_pop(population *pop, int size) {
    int i;
    pop->ind = (individual *) malloc(size * sizeof(individual));
    for (i = 0; i < size; i++) {
        allocate_memory_ind(&(pop->ind[i]));
    }
}

void deallocate_memory_ind(individual *ind) {
    int j;
    if (nreal != 0) {
        free(ind->xreal);
    }
    if (nbin != 0) {
        for (j = 0; j < nbin; j++) {
            free(ind->gene[j]);
        }
        free(ind->xbin);
        free(ind->gene);
    }
    free(ind->obj);
    if (ncon != 0) {
        free(ind->constr);
    }
}

void deallocate_memory_pop(population *pop, int size) {
    int i;
    for (i = 0; i < size; i++) {
        deallocate_memory_ind(&(pop->ind[i]));
    }
    free(pop->ind);
}

double maximum(double a, double b) {
    if (a > b) {
        return (a);
    }
    return (b);
}

double minimum(double a, double b) {
    if (a < b) {
        return (a);
    }
    return (b);
}

void realcross(individual *parent1, individual *parent2, individual *child1, individual *child2) {
    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    if (randomperc() <= pcross_real) {
        nrealcross++;
        for (i = 0; i < nreal; i++) {
            if (randomperc() <= 0.5) {
                if (fabs(parent1->xreal[i] - parent2->xreal[i]) > EPS) {
                    if (parent1->xreal[i] < parent2->xreal[i]) {
                        y1 = parent1->xreal[i];
                        y2 = parent2->xreal[i];
                    } else {
                        y1 = parent2->xreal[i];
                        y2 = parent1->xreal[i];
                    }
                    yl = min_realvar[i];
                    yu = max_realvar[i];
                    rand = randomperc();
                    beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(eta_c + 1.0));
                    if (rand <= (1.0 / alpha)) {
                        betaq = pow((rand * alpha), (1.0 / (eta_c + 1.0)));
                    } else {
                        betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                    }
                    c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                    alpha = 2.0 - pow(beta, -(eta_c + 1.0));
                    if (rand <= (1.0 / alpha)) {
                        betaq = pow((rand * alpha), (1.0 / (eta_c + 1.0)));
                    } else {
                        betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                    }
                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
                    if (c1 < yl)
                        c1 = yl;
                    if (c2 < yl)
                        c2 = yl;
                    if (c1 > yu)
                        c1 = yu;
                    if (c2 > yu)
                        c2 = yu;
                    if (randomperc() <= 0.5) {
                        child1->xreal[i] = c2;
                        child2->xreal[i] = c1;
                    } else {
                        child1->xreal[i] = c1;
                        child2->xreal[i] = c2;
                    }
                } else {
                    child1->xreal[i] = parent1->xreal[i];
                    child2->xreal[i] = parent2->xreal[i];
                }
            } else {
                child1->xreal[i] = parent1->xreal[i];
                child2->xreal[i] = parent2->xreal[i];
            }
        }
    } else {
        for (i = 0; i < nreal; i++) {
            child1->xreal[i] = parent1->xreal[i];
            child2->xreal[i] = parent2->xreal[i];
        }
    }
    return;
}

void bincross(individual *parent1, individual *parent2, individual *child1, individual *child2) {
    int i, j;
    double rand;
    int temp, site1, site2;
    for (i = 0; i < nbin; i++) {
        rand = randomperc();
        if (rand <= pcross_bin) {
            nbincross++;
            site1 = rnd(0, nbits[i] - 1);
            site2 = rnd(0, nbits[i] - 1);
            if (site1 > site2) {
                temp = site1;
                site1 = site2;
                site2 = temp;
            }
            for (j = 0; j < site1; j++) {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
            for (j = site1; j < site2; j++) {
                child1->gene[i][j] = parent2->gene[i][j];
                child2->gene[i][j] = parent1->gene[i][j];
            }
            for (j = site2; j < nbits[i]; j++) {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
        } else {
            for (j = 0; j < nbits[i]; j++) {
                child1->gene[i][j] = parent1->gene[i][j];
                child2->gene[i][j] = parent2->gene[i][j];
            }
        }
    }
    return;
}

void crossover(individual *parent1, individual *parent2, individual *child1, individual *child2) {
    if (nreal != 0) {
        realcross(parent1, parent2, child1, child2);
    }
    if (nbin != 0) {
        bincross(parent1, parent2, child1, child2);
    }
    return;
}

void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right) {
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left < right) {
        index = rnd(left, right);
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        pivot = pop->ind[obj_array[right]].obj[objcount];
        i = left - 1;
        for (j = left; j < right; j++) {
            if (pop->ind[obj_array[j]].obj[objcount] <= pivot) {
                i += 1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }
        index = i + 1;
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        q_sort_front_obj(pop, objcount, obj_array, left, index - 1);
        q_sort_front_obj(pop, objcount, obj_array, index + 1, right);
    }
    return;
}

void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size) {
    q_sort_front_obj(pop, objcount, obj_array, 0, obj_array_size - 1);
    return;
}

void assign_crowding_distance(population *pop, int *dist, int **obj_array, int front_size) {
    int i, j;
    for (i = 0; i < nobj; i++) {
        for (j = 0; j < front_size; j++) {
            obj_array[i][j] = dist[j];
        }
        quicksort_front_obj(pop, i, obj_array[i], front_size);
    }
    for (j = 0; j < front_size; j++) {
        pop->ind[dist[j]].crowd_dist = 0.0;
    }
    for (i = 0; i < nobj; i++) {
        pop->ind[obj_array[i][0]].crowd_dist = INF;
    }
    for (i = 0; i < nobj; i++) {
        for (j = 1; j < front_size - 1; j++) {
            if (pop->ind[obj_array[i][j]].crowd_dist != INF) {
                if (pop->ind[obj_array[i][front_size - 1]].obj[i] == pop->ind[obj_array[i][0]].obj[i]) {
                    pop->ind[obj_array[i][j]].crowd_dist += 0.0;
                } else {
                    pop->ind[obj_array[i][j]].crowd_dist +=
                            (pop->ind[obj_array[i][j + 1]].obj[i] - pop->ind[obj_array[i][j - 1]].obj[i]) /
                            (pop->ind[obj_array[i][front_size - 1]].obj[i] - pop->ind[obj_array[i][0]].obj[i]);
                }
            }
        }
    }
    for (j = 0; j < front_size; j++) {
        if (pop->ind[dist[j]].crowd_dist != INF) {
            pop->ind[dist[j]].crowd_dist = (pop->ind[dist[j]].crowd_dist) / nobj;
        }
    }
    return;
}

void assign_crowding_distance_indices(population *pop, int c1, int c2) {
    int **obj_array;
    int *dist;
    int i, j;
    int front_size;
    front_size = c2 - c1 + 1;
    if (front_size == 1) {
        pop->ind[c1].crowd_dist = INF;
        return;
    }
    if (front_size == 2) {
        pop->ind[c1].crowd_dist = INF;
        pop->ind[c2].crowd_dist = INF;
        return;
    }
    obj_array = (int **) malloc(nobj * sizeof(int *));
    dist = (int *) malloc(front_size * sizeof(int));
    for (i = 0; i < nobj; i++) {
        obj_array[i] = (int *) malloc(front_size * sizeof(int));
    }
    for (j = 0; j < front_size; j++) {
        dist[j] = c1++;
    }
    assign_crowding_distance(pop, dist, obj_array, front_size);
    free(dist);
    for (i = 0; i < nobj; i++) {
        free(obj_array[i]);
    }
    free(obj_array);
    return;
}

void assign_crowding_distance_list(population *pop, list *lst, int front_size) {
    int **obj_array;
    int *dist;
    int i, j;
    list *temp;
    temp = lst;
    if (front_size == 1) {
        pop->ind[lst->index].crowd_dist = INF;
        return;
    }
    if (front_size == 2) {
        pop->ind[lst->index].crowd_dist = INF;
        pop->ind[lst->child->index].crowd_dist = INF;
        return;
    }
    obj_array = (int **) malloc(nobj * sizeof(int *));
    dist = (int *) malloc(front_size * sizeof(int));
    for (i = 0; i < nobj; i++) {
        obj_array[i] = (int *) malloc(front_size * sizeof(int));
    }
    for (j = 0; j < front_size; j++) {
        dist[j] = temp->index;
        temp = temp->child;
    }
    assign_crowding_distance(pop, dist, obj_array, front_size);
    free(dist);
    for (i = 0; i < nobj; i++) {
        free(obj_array[i]);
    }
    free(obj_array);
    return;
}

void decode_ind(individual *ind) {
    int j, k;
    int sum;
    if (nbin != 0) {
        for (j = 0; j < nbin; j++) {
            sum = 0;
            for (k = 0; k < nbits[j]; k++) {
                if (ind->gene[j][k] == 1) {
                    sum += pow(2, nbits[j] - 1 - k);
                }
            }
            ind->xbin[j] =
                    min_binvar[j] + (double) sum * (max_binvar[j] - min_binvar[j]) / (double) (pow(2, nbits[j]) - 1);
        }
    }
    return;
}

void decode_pop(population *pop) {
    int i;
    if (nbin != 0) {
        for (i = 0; i < popsize; i++) {
            decode_ind(&(pop->ind[i]));
        }
    }
    return;
}

void onthefly_display(population *pop, FILE *gp, int ii) {
    int i;
    int flag;
    FILE *fpt;
    fpt = fopen("plot.out", "w");
    flag = 0;
    for (i = 0; i < popsize; i++) {
        if (pop->ind[i].constr_violation == 0) {
            if (choice != 3)
                fprintf(fpt, "%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1]);
            else
                fprintf(fpt, "%e\t%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1],
                        pop->ind[i].obj[obj3 - 1]);
            fflush(fpt);
            flag = 1;
        }
    }
    if (flag == 0) {
        // printf("\n No feasible soln in this pop, hence no display");
    } else {
        if (choice != 3){
            fprintf(gp, "set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n", ii);
        }else{
            fprintf(gp, "set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n", ii, angle1, angle2);
        }
        fflush(gp);
    }
    fclose(fpt);
    return;
}

int check_dominance(individual *a, individual *b) {
    int i;
    int flag1;
    int flag2;
    flag1 = 0;
    flag2 = 0;
    if (a->constr_violation < 0 && b->constr_violation < 0) {
        if (a->constr_violation > b->constr_violation) {
            return (1);
        } else {
            if (a->constr_violation < b->constr_violation) {
                return (-1);
            } else {
                return (0);
            }
        }
    } else {
        if (a->constr_violation < 0 && b->constr_violation == 0) {
            return (-1);
        } else {
            if (a->constr_violation == 0 && b->constr_violation < 0) {
                return (1);
            } else {
                for (i = 0; i < nobj; i++) {
                    if (a->obj[i] < b->obj[i]) {
                        flag1 = 1;

                    } else {
                        if (a->obj[i] > b->obj[i]) {
                            flag2 = 1;
                        }
                    }
                }
                if (flag1 == 1 && flag2 == 0) {
                    return (1);
                } else {
                    if (flag1 == 0 && flag2 == 1) {
                        return (-1);
                    } else {
                        return (0);
                    }
                }
            }
        }
    }
}

void evaluate_ind(individual *ind) {
    int j;
    test_problem(ind->xreal, ind->xbin, ind->gene, ind->obj, ind->constr);
    if (ncon == 0) {
        ind->constr_violation = 0.0;
    } else {
        ind->constr_violation = 0.0;
        for (j = 0; j < ncon; j++) {
            if (ind->constr[j] < 0.0) {
                ind->constr_violation += ind->constr[j];
            }
        }
    }
    return;
}

void evaluate_pop(population *pop) {
    int i;
    for (i = 0; i < popsize; i++) {
        evaluate_ind(&(pop->ind[i]));
    }
    return;
}

void insert(list *node, int x) {
    list *temp;
    if (node == NULL) {
        printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = (list *) malloc(sizeof(list));
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL) {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}

list *del(list *node) {
    list *temp;
    if (node == NULL) {
        printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child != NULL) {
        temp->child->parent = temp;
    }
    free(node);
    return (temp);
}

void q_sort_dist(population *pop, int *dist, int left, int right) {
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left < right) {
        index = rnd(left, right);
        temp = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = pop->ind[dist[right]].crowd_dist;
        i = left - 1;
        for (j = left; j < right; j++) {
            if (pop->ind[dist[j]].crowd_dist <= pivot) {
                i += 1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            }
        }
        index = i + 1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        q_sort_dist(pop, dist, left, index - 1);
        q_sort_dist(pop, dist, index + 1, right);
    }
    return;
}

void quicksort_dist(population *pop, int *dist, int front_size){
    q_sort_dist (pop, dist, 0, front_size-1);
    return;
}

void copy_ind(individual *ind1, individual *ind2) {
    int i, j;
    ind2->rank = ind1->rank;
    ind2->constr_violation = ind1->constr_violation;
    ind2->crowd_dist = ind1->crowd_dist;
    if (nreal != 0) {
        for (i = 0; i < nreal; i++) {
            ind2->xreal[i] = ind1->xreal[i];
        }
    }
    if (nbin != 0) {
        for (i = 0; i < nbin; i++) {
            ind2->xbin[i] = ind1->xbin[i];
            for (j = 0; j < nbits[i]; j++) {
                ind2->gene[i][j] = ind1->gene[i][j];
            }
        }
    }
    for (i = 0; i < nobj; i++) {
        ind2->obj[i] = ind1->obj[i];
    }
    if (ncon != 0) {
        for (i = 0; i < ncon; i++) {
            ind2->constr[i] = ind1->constr[i];
        }
    }
    return;
}

void crowding_fill(population *mixed_pop, population *new_pop, int count, int front_size, list *elite) {
    int *dist;
    list *temp;
    int i, j;
    assign_crowding_distance_list(mixed_pop, elite->child, front_size);
    dist = (int *) malloc(front_size * sizeof(int));
    temp = elite->child;
    for (j = 0; j < front_size; j++) {
        dist[j] = temp->index;
        temp = temp->child;
    }
    quicksort_dist (mixed_pop, dist, front_size);
    for (i=count, j=front_size-1; i<popsize; i++, j--)
    {
        copy_ind(&mixed_pop->ind[dist[j]], &new_pop->ind[i]);
    }
    free (dist);
    return;
}

void fill_nondominated_sort(population *mixed_pop, population *new_pop) {
    int flag;
    int i, j;
    int end;
    int front_size;
    int archieve_size;
    int rank = 1;
    list *pool;
    list *elite;
    list *temp1, *temp2;
    pool = (list *) malloc(sizeof(list));
    elite = (list *) malloc(sizeof(list));
    front_size = 0;
    archieve_size = 0;
    pool->index = -1;
    pool->parent = NULL;
    pool->child = NULL;
    elite->index = -1;
    elite->parent = NULL;
    elite->child = NULL;
    temp1 = pool;
    for (i = 0; i < 2 * popsize; i++) {
        insert(temp1, i);
        temp1 = temp1->child;
    }
    i = 0;
    do {
        temp1 = pool->child;
        insert(elite, temp1->index);
        front_size = 1;
        temp2 = elite->child;
        temp1 = del(temp1);
        temp1 = temp1->child;
        do {
            temp2 = elite->child;
            if (temp1 == NULL) {
                break;
            }
            do {
                end = 0;
                flag = check_dominance(&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
                if (flag == 1) {
                    insert(pool, temp2->index);
                    temp2 = del(temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0) {
                    temp2 = temp2->child;
                }
                if (flag == -1) {
                    end = 1;
                }
            } while (end != 1 && temp2 != NULL);
            if (flag == 0 || flag == 1) {
                insert(elite, temp1->index);
                front_size++;
                temp1 = del(temp1);
            }
            temp1 = temp1->child;
        } while (temp1 != NULL);
        temp2 = elite->child;
        j = i;
        if ((archieve_size + front_size) <= popsize) {
            do {
                copy_ind(&mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size += 1;
                temp2 = temp2->child;
                i += 1;
            } while (temp2 != NULL);
            assign_crowding_distance_indices(new_pop, j, i - 1);
            rank += 1;
        } else {
            crowding_fill(mixed_pop, new_pop, i, front_size, elite);
            archieve_size = popsize;
            for (j = i; j < popsize; j++) {
                new_pop->ind[j].rank = rank;
            }
        }
        temp2 = elite->child;
        do {
            temp2 = del(temp2);
            temp2 = temp2->child;
        } while (elite->child != NULL);
    } while (archieve_size < popsize);
    while (pool != NULL) {
        temp1 = pool;
        pool = pool->child;
        free(temp1);
    }
    while (elite != NULL) {
        temp1 = elite;
        elite = elite->child;
        free(temp1);
    }
    return;
}

void initialize_ind(individual *ind) {
    int j, k;
    if (nreal != 0) {
        for (j = 0; j < nreal; j++) {
            ind->xreal[j] = rndreal(min_realvar[j], max_realvar[j]);
        }
    }
    if (nbin != 0) {
        for (j = 0; j < nbin; j++) {
            for (k = 0; k < nbits[j]; k++) {
                if (randomperc() <= 0.5) {
                    ind->gene[j][k] = 0;
                } else {
                    ind->gene[j][k] = 1;
                }
            }
        }
    }
    return;
}

void initialize_pop(population *pop) {
    int i;
    for (i = 0; i < popsize; i++) {
        initialize_ind(&(pop->ind[i]));
    }
    return;
}

void merge(population *pop1, population *pop2, population *pop3) {
    int i, k;
    for (i = 0; i < popsize; i++) {
        copy_ind(&(pop1->ind[i]), &(pop3->ind[i]));
    }
    for (i = 0, k = popsize; i < popsize; i++, k++) {
        copy_ind(&(pop2->ind[i]), &(pop3->ind[k]));
    }
    return;
}

void real_mutate_ind(individual *ind) {
    int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
    for (j = 0; j < nreal; j++) {
        if (randomperc() <= pmut_real) {
            y = ind->xreal[j];
            yl = min_realvar[j];
            yu = max_realvar[j];
            delta1 = (y - yl) / (yu - yl);
            delta2 = (yu - y) / (yu - yl);
            rnd = randomperc();
            mut_pow = 1.0 / (eta_m + 1.0);
            if (rnd <= 0.5) {
                xy = 1.0 - delta1;
                val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (eta_m + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            } else {
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (eta_m + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }
            y = y + deltaq * (yu - yl);
            if (y < yl)
                y = yl;
            if (y > yu)
                y = yu;
            ind->xreal[j] = y;
            nrealmut += 1;
        }
    }
    return;
}

void bin_mutate_ind(individual *ind) {
    int j, k;
    double prob;
    for (j = 0; j < nbin; j++) {
        for (k = 0; k < nbits[j]; k++) {
            prob = randomperc();
            if (prob <= pmut_bin) {
                if (ind->gene[j][k] == 0) {
                    ind->gene[j][k] = 1;
                } else {
                    ind->gene[j][k] = 0;
                }
                nbinmut += 1;
            }
        }
    }
    return;
}

void mutation_ind(individual *ind) {
    if (nreal != 0) {
        real_mutate_ind(ind);
    }
    if (nbin != 0) {
        bin_mutate_ind(ind);
    }
    return;
}

void mutation_pop(population *pop) {
    int i;
    for (i = 0; i < popsize; i++) {
        mutation_ind(&(pop->ind[i]));
    }
    return;
}

void assign_rank_and_crowding_distance(population *new_pop) {
    int flag;
    int i;
    int end;
    int front_size;
    int rank = 1;
    list *orig;
    list *cur;
    list *temp1, *temp2;
    orig = (list *) malloc(sizeof(list));
    cur = (list *) malloc(sizeof(list));
    front_size = 0;
    orig->index = -1;
    orig->parent = NULL;
    orig->child = NULL;
    cur->index = -1;
    cur->parent = NULL;
    cur->child = NULL;
    temp1 = orig;
    for (i = 0; i < popsize; i++) {
        insert(temp1, i);
        temp1 = temp1->child;
    }
    do {
        if (orig->child->child == NULL) {
            new_pop->ind[orig->child->index].rank = rank;
            new_pop->ind[orig->child->index].crowd_dist = INF;
            break;
        }
        temp1 = orig->child;
        insert(cur, temp1->index);
        front_size = 1;
        temp2 = cur->child;
        temp1 = del(temp1);
        temp1 = temp1->child;
        do {
            temp2 = cur->child;
            do {
                end = 0;
                flag = check_dominance(&(new_pop->ind[temp1->index]), &(new_pop->ind[temp2->index]));
                if (flag == 1) {
                    insert(orig, temp2->index);
                    temp2 = del(temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0) {
                    temp2 = temp2->child;
                }
                if (flag == -1) {
                    end = 1;
                }
            } while (end != 1 && temp2 != NULL);
            if (flag == 0 || flag == 1) {
                insert(cur, temp1->index);
                front_size++;
                temp1 = del(temp1);
            }
            temp1 = temp1->child;
        } while (temp1 != NULL);
        temp2 = cur->child;
        do {
            new_pop->ind[temp2->index].rank = rank;
            temp2 = temp2->child;
        } while (temp2 != NULL);
        assign_crowding_distance_list(new_pop, cur->child, front_size);
        temp2 = cur->child;
        do {
            temp2 = del(temp2);
            temp2 = temp2->child;
        } while (cur->child != NULL);
        rank += 1;
    } while (orig->child != NULL);
    free(orig);
    free(cur);
    return;
}

void report_pop(population *pop, FILE *fpt) {
    int i, j, k;
    for (i = 0; i < popsize; i++) {
        for (j = 0; j < nobj; j++) {
            fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
        }
        if (ncon != 0) {
            for (j = 0; j < ncon; j++) {
                fprintf(fpt, "%e\t", pop->ind[i].constr[j]);
            }
        }
        if (nreal != 0) {
            for (j = 0; j < nreal; j++) {
                fprintf(fpt, "%e\t", pop->ind[i].xreal[j]);
            }
        }
        if (nbin != 0) {
            for (j = 0; j < nbin; j++) {
                for (k = 0; k < nbits[j]; k++) {
                    fprintf(fpt, "%d\t", pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt, "%e\t", pop->ind[i].constr_violation);
        fprintf(fpt, "%d\t", pop->ind[i].rank);
        fprintf(fpt, "%e\n", pop->ind[i].crowd_dist);
    }
    return;
}

void report_feasible(population *pop, FILE *fpt) {
    int i, j, k;
    for (i = 0; i < popsize; i++) {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank == 1) {
            for (j = 0; j < nobj; j++) {
                fprintf(fpt, "%3.5f\t", pop->ind[i].obj[j]);
            }

            fprintf(fpt, "\t\t");

            for (j = 0; j < ncon; j++) {
                fprintf(fpt, "%3.5f\t", pop->ind[i].constr[j]);
            }

            fprintf(fpt, "\t\t");

            for (j = 0; j < nreal; j++) {
                fprintf(fpt, "%3.5f\t", pop->ind[i].xreal[j]);
            }

            //fprintf(fpt, "%e\t", pop->ind[i].constr_violation);
            //fprintf(fpt, "%d\t", pop->ind[i].rank);
            //fprintf(fpt, "%e", pop->ind[i].crowd_dist);
            fprintf(fpt, "\n");
        }
    }
    return;
}

individual *tournament(individual *ind1, individual *ind2) {
    int flag;
    flag = check_dominance(ind1, ind2);
    if (flag == 1) {
        return (ind1);
    }
    if (flag == -1) {
        return (ind2);
    }
    if (ind1->crowd_dist > ind2->crowd_dist) {
        return (ind1);
    }
    if (ind2->crowd_dist > ind1->crowd_dist) {
        return (ind2);
    }
    if ((randomperc()) <= 0.5) {
        return (ind1);
    } else {
        return (ind2);
    }
}

void selection(population *old_pop, population *new_pop) {
    int *a1, *a2;
    int temp;
    int i;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *) malloc(popsize * sizeof(int));
    a2 = (int *) malloc(popsize * sizeof(int));
    for (i = 0; i < popsize; i++) {
        a1[i] = a2[i] = i;
    }
    for (i = 0; i < popsize; i++) {
        rand = rnd(i, popsize - 1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rnd(i, popsize - 1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    for (i = 0; i < popsize; i += 4) {
        parent1 = tournament(&old_pop->ind[a1[i]], &old_pop->ind[a1[i + 1]]);
        parent2 = tournament(&old_pop->ind[a1[i + 2]], &old_pop->ind[a1[i + 3]]);
        crossover(parent1, parent2, &new_pop->ind[i], &new_pop->ind[i + 1]);
        parent1 = tournament(&old_pop->ind[a2[i]], &old_pop->ind[a2[i + 1]]);
        parent2 = tournament(&old_pop->ind[a2[i + 2]], &old_pop->ind[a2[i + 3]]);
        crossover(parent1, parent2, &new_pop->ind[i + 2], &new_pop->ind[i + 3]);
    }
    free(a1);
    free(a2);
    return;
}

int main() {

    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *gp;

    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    inicializaTeste();

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    fpt1 = fopen("initial_pop.out", "w");
    fpt2 = fopen("final_pop.out", "w");
    fpt3 = fopen("best_pop.out", "w");
    fpt4 = fopen("all_pop.out", "w");
    fpt5 = fopen("params.out", "w");

    gp = popen(GNUPLOT_COMMAND, "w");

    printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5, "\n Population size = %d", popsize);
    fprintf(fpt5, "\n Number of generations = %d", ngen);
    fprintf(fpt5, "\n Number of objective functions = %d", nobj);
    fprintf(fpt5, "\n Number of constraints = %d", ncon);
    fprintf(fpt5, "\n Number of real variables = %d", nreal);

    for (int i = 0; i < nreal; i++) {
        fprintf(fpt5, "\n Lower limit of real variable %d = %e", i + 1, min_realvar[i]);
        fprintf(fpt5, "\n Upper limit of real variable %d = %e", i + 1, max_realvar[i]);
    }

    fprintf(fpt5, "\n Probability of crossover of real variable = %e", pcross_real);
    fprintf(fpt5, "\n Probability of mutation of real variable = %e", pmut_real);
    fprintf(fpt5, "\n Distribution index for crossover = %e", eta_c);
    fprintf(fpt5, "\n Distribution index for mutation = %e", eta_m);

    fprintf(fpt5, "\n Number of binary variables = %d", nbin);

    fprintf(fpt5, "\n Seed for random number generator = %e", seed);

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

    fflush(stdout);

    onthefly_display(parent_pop, gp, 1);

    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);

    sleep(1);

    for (int i = 2; i <= ngen; i++) {
        selection(parent_pop, child_pop);

        mutation_pop(child_pop);

        decode_pop(child_pop);

        evaluate_pop(child_pop);

        merge(parent_pop, child_pop, mixed_pop);

        fill_nondominated_sort(mixed_pop, parent_pop);
        // Comment following four lines if information for all generations is not desired, it will speed up the execution
        //fprintf(fpt4,"# gen = %d\n",i);
        //report_pop(parent_pop,fpt4);
        //fflush(fpt4);
        onthefly_display(parent_pop, gp, i);
        // printf("\n gen = %d", i);
    }

    printf("\n Generations finished, now reporting solutions");
    report_pop(parent_pop, fpt2);
    report_feasible(parent_pop, fpt3);

    fprintf(fpt5, "\n Number of crossover of real variable = %d", nrealcross);
    fprintf(fpt5, "\n Number of mutation of real variable = %d", nrealmut);

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

    pclose(gp);

    free(min_realvar);
    free(max_realvar);

    deallocate_memory_pop(parent_pop, popsize);
    deallocate_memory_pop(child_pop, popsize);
    deallocate_memory_pop(mixed_pop, 2 * popsize);

    free(parent_pop);
    free(child_pop);
    free(mixed_pop);

    printf("\n Routine successfully exited \n");

    return 0;
}