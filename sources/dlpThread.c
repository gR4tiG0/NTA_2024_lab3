#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <pthread.h>

const char* TR_PRIMES[] = {"2", "3", "5", "7", "11", "13", "17", "19", "23", "29", "31", "37", "41", "43", "47"};
const int NUM_PRIMES = sizeof(TR_PRIMES) / sizeof(TR_PRIMES[0]);


void mpz_power(mpz_t a, mpz_t b, mpz_t p, mpz_t result) {
    mpz_set_ui(result, 1);
    mpz_t base, temp, dec;
    mpz_init_set(base, a);
    mpz_init(temp);
    mpz_init_set(dec, b);
    mpz_mod(base, base, p); // base = base % p


    while (mpz_cmp_ui(dec, 0) > 0) {
        if (mpz_odd_p(dec)) {
            mpz_mul(temp, result, base); // temp = result * base
            mpz_mod(result, temp, p); // result = temp % p
        }
        mpz_tdiv_q_ui(dec, dec, 2); // b = b / 2
        mpz_mul(temp, base, base); // temp = base * base
        mpz_mod(base, temp, p); // base = temp % p
    }

    mpz_clear(base);
    mpz_clear(temp);
    mpz_clear(dec);
}

uint64_t power(uint64_t a, uint64_t b, uint64_t p) {
    mpz_t mpz_a, mpz_b, mpz_p, result;
    mpz_init_set_ui(mpz_a, a);
    mpz_init_set_ui(mpz_b, b);
    mpz_init_set_ui(mpz_p, p);
    mpz_init_set_ui(result,1);

    mpz_power(mpz_a, mpz_b, mpz_p, result);

    uint64_t res = mpz_get_ui(result);
    mpz_clear(mpz_a);
    mpz_clear(mpz_b);
    mpz_clear(mpz_p);
    mpz_clear(result);
    return res;
}

void xgcd(mpz_t a, mpz_t b, mpz_t g, mpz_t x, mpz_t y) {
    if (mpz_cmp_ui(a, 0) == 0) {
        mpz_set_ui(x, 0);
        mpz_set_ui(y, 1);
        mpz_set(g, b);
    } else {
        mpz_t x0, y0, temp;
        mpz_init(x0);
        mpz_init(y0);
        mpz_init(temp);

        mpz_mod(temp, b, a);
        xgcd(temp, a, g, x0, y0);

        mpz_fdiv_q(temp, b, a);
        mpz_mul(temp, temp, x0);
        mpz_sub(x, y0, temp);
        mpz_set(y, x0);

        mpz_clear(x0);
        mpz_clear(y0);
        mpz_clear(temp);
    }
}

void mpz_inv(mpz_t a, mpz_t p, mpz_t result) {
    mpz_t x, y, g;
    mpz_init(x);
    mpz_init(y);
    mpz_init(g);

    xgcd(a, p, g, x, y);
    if (mpz_cmp_ui(g, 1) != 0) {
        // gmp_printf("a = %Zd, p = %Zd\n", a, p);
        printf("Inverse does not exist\n");
        mpz_set_ui(result, 0);
    } else {
        if (mpz_cmp_si(x, 0) < 0) {
            mpz_add(x, x, p);
        }
        mpz_mod(result, x, p);
    }

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(g);
}

uint64_t inv(uint64_t a, uint64_t p) {
    mpz_t mpz_a, mpz_p, result;
    mpz_init_set_ui(mpz_a, a);
    mpz_init_set_ui(mpz_p, p);
    mpz_init(result);

    mpz_inv(mpz_a, mpz_p, result);

    uint64_t res = mpz_get_ui(result);
    mpz_clear(mpz_a);
    mpz_clear(mpz_p);
    mpz_clear(result);
    return res;
}

void rhoPolard(mpz_t n, mpz_t x0, mpz_t result) {
    mpz_t x, y, d, temp;
    mpz_inits(x, y, d, temp, NULL);
    mpz_set(x, x0);
    mpz_set(y, x0);
    mpz_set_ui(d, 1);

    while (mpz_cmp_ui(d, 1) == 0) {
        mpz_mul(temp, x, x);
        mpz_add_ui(temp, temp, 1);
        mpz_mod(x, temp, n);

        mpz_mul(temp, y, y);
        mpz_add_ui(temp, temp, 1);
        mpz_mod(temp, temp, n);
        mpz_mul(temp, temp, temp);
        mpz_add_ui(temp, temp, 1);
        mpz_mod(y, temp, n);

        mpz_sub(temp, x, y);
        mpz_abs(temp, temp);
        mpz_gcd(d, temp, n);
    }

    mpz_set(result, d);

    mpz_clears(x, y, d, temp, NULL);
}

void rpFactor(mpz_t n, mpz_t result) {
    mpz_t x0;
    mpz_init_set_ui(x0, 2);
    rhoPolard(n, x0, result);
    mpz_clear(x0);
}

uint64_t factor(uint64_t n) {
    mpz_t mpz_n, result;
    mpz_init_set_ui(mpz_n, n);
    mpz_init(result);

    rpFactor(mpz_n, result);

    uint64_t res = mpz_get_ui(result);
    mpz_clear(mpz_n);
    mpz_clear(result);
    return res;
}

int mrIter(mpz_t base, mpz_t num) {
    int power, iterations;
    mpz_t baseRaised, numMinusOne, temp;
    mpz_init(numMinusOne);
    mpz_sub_ui(numMinusOne, num, 1);

    power = 0;
    mpz_init_set(temp, numMinusOne);
    while (mpz_even_p(temp)) {
        mpz_fdiv_q_2exp(temp, temp, 1);
        power++;
    }

    mpz_init(baseRaised);
    mpz_powm(baseRaised, base, temp, num);

    if (mpz_cmp_ui(baseRaised, 1) == 0) {
        mpz_clears(baseRaised, temp, numMinusOne, NULL);
        return 1;
    }

    for(iterations = 0; iterations < power - 1; iterations++) {
        if (mpz_cmp(baseRaised, numMinusOne) == 0) {
            mpz_clears(baseRaised, temp, numMinusOne, NULL);
            return 1;
        }
        mpz_powm_ui(baseRaised, baseRaised, 2, num);
    }

    if (mpz_cmp(baseRaised, numMinusOne) == 0) {
        mpz_clears(baseRaised, temp, numMinusOne, NULL);
        return 1;
    }

    mpz_clears(baseRaised, temp, numMinusOne, NULL);
    return 0;
}
int isMillerRabin(mpz_t num, gmp_randstate_t randState) {
    if (mpz_cmp_ui(num, 1) == 0) {
        return 1;
    }
    mpz_t randomNum;
    int repeat;
    mpz_init(randomNum);
    // gmp_printf("we are here %Zd\n", num);
    for(repeat = 0; repeat < 20; repeat++) {
        
        do {
            mpz_urandomm(randomNum, randState, num);
        } while (mpz_sgn(randomNum) == 0);
        // if (repeat == 0) printf("we are here %i\n", repeat);
        if (mrIter(randomNum, num) == 0) {
            mpz_clear(randomNum);
            return 0;
        }
    }
    mpz_clear(randomNum);
    return 1;
}
void lastResort(mpz_t n, mpz_t result) {
    mpz_t i;
    mpz_init_set_ui(i, 2);
    mpz_t temp;
    mpz_init(temp);
    mpz_set(result, n);

    while(1) {
        mpz_mod(temp, n, i);
        if (mpz_cmp_ui(temp, 0) == 0) {
            mpz_set(result, i);
            mpz_clears(i, temp, NULL);
            return;
        }
        mpz_add_ui(i, i, 1);
    }
    mpz_clears(i, temp, NULL);
}

void trivialFactor(mpz_t result, const mpz_t n) {
    mpz_t temp;
    mpz_init(temp);

    for (int i = 0; i < NUM_PRIMES; i++) {
        mpz_set_str(temp, TR_PRIMES[i], 10);
        if (mpz_divisible_p(n, temp)) {
            mpz_set(result, temp);
            mpz_clear(temp);
            return;
        }
    }

    mpz_set(result, n);
    mpz_clear(temp);
}



uint64_t isPrime(uint64_t n) {
    mpz_t mpz_n;
    mpz_init_set_ui(mpz_n, n);
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    int res = isMillerRabin(mpz_n, rand_state);
    gmp_randclear(rand_state);
    mpz_clear(mpz_n);
    return res;
}

typedef struct {
    mpz_t factor;
    mpz_t power;
} FactorPowerPair;



FactorPowerPair* factorize(mpz_t mpz_n, int* numFactors) {
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    // gmp_printf("Number to factorize: %Zd\n", mpz_n);
    FactorPowerPair* factors = malloc(sizeof(FactorPowerPair) * mpz_sizeinbase(mpz_n, 2));
    *numFactors = 0;

    if (isMillerRabin(mpz_n, rand_state)) {
        mpz_init_set(factors[*numFactors].factor, mpz_n);
        mpz_init_set_ui(factors[*numFactors].power, 1);
        (*numFactors)++;
    } else {
        while (1) {
            mpz_t factor;
            mpz_init(factor);
            rpFactor(mpz_n, factor);
            // trivialFactor(factor, mpz_n);

            mpz_t count;
            mpz_init_set_ui(count, 1);
            // printf("Factor: ");
            // gmp_printf("%Zd\n", factor);
            while (!isMillerRabin(factor, rand_state)){
                // gmp_printf("got not common factor: %Zd\n", factor);
                mpz_t step_res;
                mpz_init(step_res);
                trivialFactor(step_res, factor);
                if (mpz_cmp(step_res, factor) == 0) {
                    // gmp_printf("step_res == factor; %Zd == %Zd\n", step_res, factor);
                    rpFactor(factor, step_res);
                    // gmp_printf("new factor: %Zd\n", step_res);
                }
                if (mpz_cmp(step_res, factor) == 0) {
                    lastResort(factor, step_res);
                }
                mpz_set(factor, step_res);
                mpz_clear(step_res);   
                // gmp_printf("Factor: %Zd\n", factor);
            }
            
            // gmp_printf("factor: %Zd\n", factor);

            // Check if factor is already in the factors array
            int found = 0;
            for (int i = 0; i < *numFactors; i++) {
                if (mpz_cmp(factors[i].factor, factor) == 0) {
                    mpz_add_ui(factors[i].power, factors[i].power, 1);
                    found = 1;
                    break;
                }
            }

            // If factor is not in the factors array, add it
            if (!found) {
                mpz_init_set(factors[*numFactors].factor, factor);
                mpz_init_set(factors[*numFactors].power, count);
                (*numFactors)++;
            }

            mpz_divexact(mpz_n, mpz_n, factor);
            if (isMillerRabin(mpz_n, rand_state)) {
                // printf("n is prime\n");
                // gmp_printf("n: %Zd\n", mpz_n);
                int found = 0;
                for (int i = 0; i < *numFactors; i++) {
                    if (mpz_cmp(factors[i].factor, mpz_n) == 0) {
                        mpz_add_ui(factors[i].power, factors[i].power, 1);
                        found = 1;
                        break;
                    }
                }

                // If factor is not in the factors array, add it
                if (!found) {
                    mpz_init_set(factors[*numFactors].factor, mpz_n);
                    mpz_init_set_ui(factors[*numFactors].power, 1);
                    (*numFactors)++;
                }
                break;
            }

            mpz_clear(count);
        }
    }

    gmp_randclear(rand_state);

    return factors;
}

void factorizeAndPrint(uint64_t n) {
    mpz_t mpz_n;
    mpz_init_set_ui(mpz_n, n);

    int numFactors;
    FactorPowerPair* factors = factorize(mpz_n, &numFactors);

    gmp_printf("%lu = ", n);
    for (int i = 0; i < numFactors; i++) {
        gmp_printf("%Zd^%Zd", factors[i].factor, factors[i].power);
        if (i < numFactors - 1) {
            gmp_printf(" * ");
        }
    }
    gmp_printf("\n");

    for (int i = 0; i < numFactors; i++) {
        mpz_clears(factors[i].factor, factors[i].power, NULL);
    }
    free(factors);
    mpz_clear(mpz_n);
}

void getMax(FactorPowerPair *factors, size_t fLen, mpz_t maxF) {
    for (size_t i = 0; i < fLen; i++) {
        if (mpz_cmp(factors[i].factor, maxF) > 0) {
            // gmp_printf("maxF: %Zd\n", maxF);
            mpz_set(maxF, factors[i].factor);
        }
    }
}

//lab3



void gaussian_elim(mpz_t **matrix, mpz_t *solution_array, mpz_t modulus, int height, int width, mpz_t *result) {
    // gmp_printf("mod = %Zd\n", modulus);
    mpz_t temp, r, res;
    mpz_init(temp);
    mpz_init(r);
    mpz_init(res);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            mpz_set(temp, matrix[j][i]);
            // gmp_printf("i = %i, j = %i, matrix[j][i] = %Zd\n", i, j, temp);
            mpz_t x, y, g;
            mpz_init(x);
            mpz_init(y);
            mpz_init(g);
            xgcd(temp, modulus, g, x, y);
            int cm = mpz_cmp_ui(g, 1);
            // gmp_printf("g = %Zd\n", g);
            mpz_clears(x, y, g, NULL);
            if (cm != 0) {
                continue;
            }

            mpz_inv(temp, modulus, r);
            // gmp_printf("r = %Zd\n", r);
            for (int l = 0; l < width; l++) {
                mpz_mul(res, matrix[j][l], r);
                mpz_mod(matrix[j][l], res, modulus);
            }
            // gmp_printf("sol arr j %Zd\n", solution_array[j]);
            mpz_mul(res, solution_array[j], r); 
            // gmp_printf("res %Zd\n", res);
            mpz_mod(solution_array[j], res, modulus);
            // gmp_printf("sol arr j %Zd\n", solution_array[j]);
            for (int l = 0; l < height; l++) {
                if (l == j) {
                    continue;
                }
                mpz_set(r, matrix[l][i]);
                if (mpz_cmp_ui(r, 0) != 0) {
                    for (int k = 0; k < width; k++) {
                        mpz_mul(temp, matrix[j][k], r);
                        mpz_sub(res, matrix[l][k], temp);
                        mpz_mod(matrix[l][k], res, modulus);
                    }
                    mpz_mul(temp, solution_array[j], r);
                    mpz_sub(solution_array[l], solution_array[l], temp);
                    mpz_mod(solution_array[l], solution_array[l], modulus);
                }
            }
            //log ---------------------------------------------
            // printf("Matrix:\n");
            // for (int it = 0; it < height; it++) {
            //     for (int jt = 0; jt < width; jt++) {
            //         gmp_printf("%Zd ", matrix[it][jt]);
            //     }
            //     gmp_printf(" = %Zd", solution_array[it]);
            //     printf("\n");
            // }
            // printf("\n");
            
            break;
        }
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (mpz_cmp_ui(matrix[i][j], 1) == 0) {
                mpz_set(result[j], solution_array[i]);
                break;
            }
        }
    }
    
    mpz_clears(temp, r, res, NULL);
}


void generateFactorBase(mpz_t n, mpz_t *result, size_t *result_len){ //works
    double e = exp(1);
    double n_log = mpz_get_d(n);
    double b = round(exp(0.5 * pow(log(n_log) * log(log(n_log)), 0.5)) * 3.38);

    *result_len = 0;
    for (int i = 2; i < b; i++) {
        if (isPrime(i)) {
            mpz_init_set_ui(result[*result_len], i);
            (*result_len)++;
        }
    }
}

int binary_search(mpz_t *array, int len, mpz_t value) {
    int left = 0;
    int right = len - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;

        if (mpz_cmp(array[mid], value) == 0) {
            return mid;
        } else if (mpz_cmp(array[mid], value) < 0) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return -1;
}

void clearFactorStruct(FactorPowerPair *factors, int numFactors) {
    for (int i = 0; i < numFactors; i++) {
        mpz_clear(factors[i].factor);
        mpz_clear(factors[i].power);
    }
    free(factors);
}

typedef struct {
    int id;
    int start;
    int end;
    mpz_t **matrix;
    mpz_t *solutions_array;
    size_t factor_base_len;
    mpz_t a;
    mpz_t p;
    mpz_t n;
    mpz_t *factor_base;
    size_t matrix_len;
} ThreadDataGen;



void* generateSubMatrix(void *arg) {
    ThreadDataGen *data = (ThreadDataGen*) arg;

    // printf("Thread %d started\n", data->id);
    
    size_t factor_base_len = data->factor_base_len;
    mpz_t a, p, n;
    mpz_init_set(a, data->a);
    mpz_init_set(p, data->p);
    mpz_init_set(n, data->n);
    mpz_t *factor_base = data->factor_base;

    mpz_t **matrix = data->matrix;
    mpz_t *solutions_array = data->solutions_array;

    gmp_randstate_t state;
    gmp_randinit_default(state);    
    gmp_randseed_ui(state, time(NULL)+(data->id * 1000));
    
    size_t matrix_len = data->matrix_len;
    

    int mCount = data->start;
    int end = data->end;
    //generating equations

    //inits
    mpz_t k, aK;
    mpz_init(k);
    mpz_init(aK);
    while (1) {

        //get random k
        mpz_urandomm(k, state, n);
        // gmp_printf("k: %Zd\n", k);

        //get aK = pow(a,k,p)
        mpz_power(a, k, p, aK);


        //log ---------------------------------------------
        // gmp_printf("aK: %Zd\n", aK);
        //end of log---------------------------------------
        //factorize aK
        int numFactors;
        FactorPowerPair* factors = factorize(aK, &numFactors);

        //log ---------------------------------------------
        // printf("Factors:\n");
        // for (int i = 0; i < numFactors; i++) {
        //     gmp_printf("%Zd^%Zd ", factors[i].factor, factors[i].power);
        // }
        // printf("\n");
        //end of log---------------------------------------
        //basic checks

        //log ---------------------------------------------
        //end of log---------------------------------------

        mpz_t maxF;
        mpz_init(maxF);
        getMax(factors, numFactors, maxF);
        //log ---------------------------------------------
        // printf("Max factor: ");
        // gmp_printf("%Zd\n", maxF);
        // gmp_printf("%Zd\n", factor_base[factor_base_len - 1]);
        //end of log---------------------------------------
        
        if (mpz_cmp(maxF, factor_base[factor_base_len - 1]) > 0) {
            clearFactorStruct(factors, numFactors);
            // break;
            continue;
        }
        mpz_clear(maxF);

        //log ---------------------------------------------
        // printf("check passed\n");
        // // printf("Factors:\n");
        // for (int i = 0; i < numFactors; i++) {
        //     gmp_printf("%Zd^%Zd ", factors[i].factor, factors[i].power);
        // }
        // printf("\n");
        //end of log---------------------------------------

        //check if flat
        int all_in_base = 1;
        for (int i = 0; i < numFactors; i++){
            
            int in_base = 1;
            int indx = binary_search(factor_base, factor_base_len, factors[i].factor);
            if (indx == -1) {
                in_base = 0;
            }

            if (!in_base) {
                all_in_base = 0;
                break;
            }
        }
        if (all_in_base) {
            for (int i = 0; i < factor_base_len; i++) {
                mpz_set_ui(matrix[mCount][i], 0);
            }

            for (int i = 0; i < numFactors; i++) {
                int indx = binary_search(factor_base, factor_base_len, factors[i].factor);
                mpz_set(matrix[mCount][indx], factors[i].power);
            }
            mpz_set(solutions_array[mCount], k);
            mCount++;
        } else {
            continue;
            // clearFactorStruct(factors, numFactors);
            // break;
        }



        clearFactorStruct(factors, numFactors);
        if (mCount == end) {
            break;
        }
    }
    mpz_clear(aK);
    mpz_clear(k);
    gmp_randclear(state);   
    // printf("Thread %d finished\n", data->id);
    return NULL;
}

typedef struct {
    int id;
    mpz_t a;
    mpz_t b;
    mpz_t p;
    mpz_t n;
    mpz_t* factor_base;
    mpz_t* solution;
    size_t factor_base_len;
    int *solution_found;
    mpz_t* result;
    pthread_mutex_t *mutex;
} ThreadDataSolve;


void* findSolution(void *arg) {
    ThreadDataSolve *data = (ThreadDataSolve*) arg;
    // printf("Thread %d started\n", data->id);


    mpz_t a, b, p, n;

    mpz_init_set(a, data->a);
    mpz_init_set(b, data->b);
    mpz_init_set(p, data->p);
    mpz_init_set(n, data->n);

    mpz_t* factor_base = data->factor_base;
    mpz_t* solution = data->solution;
    size_t factor_base_len = data->factor_base_len;


    int solution_is_found = 0;
    mpz_t tmp1,tmp2;
    mpz_init_set_ui(tmp1, 1);
    mpz_init_set_ui(tmp2, 1);
    gmp_randstate_t state;
    gmp_randinit_default(state);    
    gmp_randseed_ui(state, time(NULL)* (data->id * 1000));
    mpz_t k;
    mpz_init(k);
    mpz_urandomm(k, state, n);
    gmp_printf("Thread %d started with random number: %Zd\n", data->id, k);
    uint64_t counterTMP = 0;
    while (1) {
        pthread_mutex_lock(data->mutex);
        if (*data->solution_found) {
            pthread_mutex_unlock(data->mutex);
            break;
        }
        pthread_mutex_unlock(data->mutex);

        counterTMP++;
        if (counterTMP % 100000 == 0) printf("C%i:%lu\r", data->id, counterTMP);
        // printf("C:%i\n", counterTMP);

        mpz_urandomm(k, state, n);
        mpz_power(a, k, p, tmp1);
        mpz_mul(tmp2, b, tmp1);
        mpz_mod(tmp2, tmp2, p);
        int numFactors;
        FactorPowerPair* factors = factorize(tmp2, &numFactors);
        //log ---------------------------------------------
        // printf("Factors:\n");
        // for (int i = 0; i < numFactors; i++) {
        //     gmp_printf("%Zd^%Zd ", factors[i].factor, factors[i].power);
            
        // }
        // printf("\n");
        // fflush(stdout);
        //end of log---------------------------------------
        mpz_t maxF;
        mpz_init(maxF);
        getMax(factors, numFactors, maxF);
        //log ---------------------------------------------
        // printf("Max factor: ");
        // gmp_printf("%Zd\n", maxF);
        // gmp_printf("%Zd\n", factor_base[factor_base_len - 1]);
        // fflush(stdout);
        //end of log---------------------------------------
        
        if (mpz_cmp(maxF, factor_base[factor_base_len - 1]) > 0) {
            
            clearFactorStruct(factors, numFactors);
            mpz_clear(maxF);
            // printf("cleared struct\n");
            continue;
        }
        mpz_clear(maxF);

        int all_in_base = 1;
        for (int i = 0; i < numFactors; i++){
            
            int in_base = 1;
            int indx = binary_search(factor_base, factor_base_len, factors[i].factor);
            if (indx == -1) {
                in_base = 0;
            }

            if (!in_base) {
                all_in_base = 0;
                break;
            }
        }
        if (!all_in_base) {
            clearFactorStruct(factors, numFactors);
            continue;
        }


        mpz_t rlog;
        mpz_init_set_ui(rlog, 0);
        for (int i = 0; i < numFactors; i++) {
            int indx = binary_search(factor_base, factor_base_len, factors[i].factor);
            mpz_mul(tmp1, factors[i].power, solution[indx]);
            mpz_add(rlog, rlog, tmp1);
        }
        mpz_sub(rlog, rlog, k);
        mpz_mod(rlog, rlog, n);
        clearFactorStruct(factors, numFactors);
        mpz_power(a, rlog, p, tmp2);



        //log ---------------------------------------------
        // printf("tmpr\n");
        // gmp_printf("a = %Zd; b = %Zd; p = %Zd; x = %Zd; r = %Zd\n",a, b, p, rlog, tmp2);
        //end of log---------------------------------------

        if (mpz_cmp(tmp2, b) == 0) {
            solution_is_found = 1;
        }

    

        if (solution_is_found) {
            pthread_mutex_lock(data->mutex);
            *data->solution_found = 1;
            mpz_set(data->result[0], rlog);
            gmp_printf("Thread %d found solution %Zd\n", data->id, data->result);
            mpz_clear(rlog);
            pthread_mutex_unlock(data->mutex);
            break;
        }
    }

    mpz_clears(tmp1, tmp2, NULL);
    mpz_clear(k);
    gmp_randclear(state);  
    mpz_clears(a, b, p, n, NULL);
    // printf("Thread %d ended\n", data->id);
    return NULL;
}


void mpzIndexCalculus(mpz_t a, mpz_t b, mpz_t p, mpz_t result) {
    //create factor base

    int numThreads = 4;

    mpz_t n;
    mpz_init_set(n, p);
    mpz_sub_ui(n, n, 1);
    
    size_t max_fb_size = 1000; //limit of factor base size
    mpz_t *factor_base = malloc(max_fb_size * sizeof(mpz_t));
    size_t factor_base_len = 0;
    generateFactorBase(p, factor_base, &factor_base_len);

    //log ---------------------------------------------
    // printf("Factor base: ");
    // for (size_t i = 0; i < factor_base_len; i++) {
    //     gmp_printf("%Zd ", factor_base[i]);
    // }
    // printf("\n");
    //end of log---------------------------------------
    //end of factor base creation OK

    //thread mod started

    size_t mC = 15 + round(factor_base_len/5);
    size_t matrix_len = factor_base_len + mC;
    mpz_t **matrix = malloc(matrix_len * sizeof(mpz_t*));
    for (int i = 0; i < matrix_len; i++) {
        matrix[i] = malloc(factor_base_len * sizeof(mpz_t));
        for (int j = 0; j < factor_base_len; j++) {
            mpz_init(matrix[i][j]);
        }
    }
    

    mpz_t *solutions_array = malloc(matrix_len * sizeof(mpz_t));
    for (int i = 0; i < matrix_len; i++) {
        mpz_init(solutions_array[i]);
    }


    ThreadDataGen thread_data[numThreads];
    pthread_t threads[numThreads];

    int submatrix_len = matrix_len / numThreads;

    printf("Generating equations\n");
    for (int i = 0; i < numThreads; i++) {
        thread_data[i].id = i;
        thread_data[i].start = i * submatrix_len;
        thread_data[i].end = (i + 1) * submatrix_len;
        thread_data[i].matrix = matrix;
        thread_data[i].solutions_array = solutions_array;
        thread_data[i].factor_base_len = factor_base_len;
        mpz_init_set(thread_data[i].a, a);
        mpz_init_set(thread_data[i].p, p);
        mpz_init_set(thread_data[i].n, n);
        thread_data[i].factor_base = factor_base;
        thread_data[i].matrix_len = matrix_len;
        pthread_create(&threads[i], NULL, generateSubMatrix, &thread_data[i]);
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }
    //thread mod ended

    //log ---------------------------------------------

    printf("equations generated\n");
    // printf("Matrix:\n");
    // for (int i = 0; i < matrix_len; i++) {
    //     for (int j = 0; j < factor_base_len; j++) {
    //         gmp_printf("%Zd ", matrix[i][j]);
    //     }
    //     printf("\n");
    // }
    // Print the solutions array
    // printf("Solutions array:\n");
    // for (int i = 0; i < matrix_len; i++) {
    //     gmp_printf("%Zd ", solutions_array[i]);
    // }
    // printf("\n");
    //end of log---------------------------------------
    mpz_t *solution = malloc(factor_base_len * sizeof(mpz_t));
    for (int i = 0; i < factor_base_len; i++) {
        mpz_init(solution[i]);
    }
    gaussian_elim(matrix, solutions_array, n, matrix_len, factor_base_len, solution);

    printf("equations solved\n");
    //log ---------------------------------------------
    // printf("solution:\n");
    // for (int i = 0; i < factor_base_len; i++) {
    //     gmp_printf("%Zd ", solution[i]);
    // }
    // printf("\n");
    //end of log---------------------------------------

    // return;
    printf("starting solution search\n");
    ThreadDataSolve thread_data_s[numThreads];

    int solution_found = 0;
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);

    mpz_t* thread_result = malloc(sizeof(mpz_t));
    mpz_init(thread_result[0]);

    for (int i = 0; i < numThreads; i++) {
        thread_data_s[i].id = i;
        thread_data_s[i].result = thread_result; 
        thread_data_s[i].solution_found = &solution_found;
        thread_data_s[i].mutex = &mutex;
        mpz_init_set(thread_data_s[i].a, a);
        mpz_init_set(thread_data_s[i].b, b);
        mpz_init_set(thread_data_s[i].p, p);
        mpz_init_set(thread_data_s[i].n, n);
        thread_data_s[i].factor_base = factor_base;
        thread_data_s[i].solution = solution;
        thread_data_s[i].factor_base_len = factor_base_len;
        pthread_create(&threads[i], NULL, findSolution, &thread_data_s[i]);
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    mpz_set(result, thread_result[0]);
    mpz_clear(thread_result[0]);
    free(thread_result);


    printf("solution found\n");
    // matrix is now complete

    pthread_mutex_destroy(&mutex);


    for (int i = 0; i < factor_base_len; i++) {
        mpz_clear(solution[i]);
    }
    free(solution);


    //start clearing
    //clear matrix and solution space
    for (int i = 0; i < matrix_len; i++) {
        for (int j = 0; j < factor_base_len; j++) {
            mpz_clear(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);

    // Clear and free the solutions array
    for (int i = 0; i < matrix_len; i++) {
        mpz_clear(solutions_array[i]);
    }
    free(solutions_array);

    //clear factor base
    for (size_t i = 0; i < factor_base_len; i++) {
        mpz_clear(factor_base[i]);
    }
    free(factor_base); 
}


uint64_t indexCalculus(uint64_t a, uint64_t b, uint64_t p) {
    mpz_t mpz_a, mpz_b, mpz_p, result;
    mpz_init_set_ui(mpz_a, a);
    mpz_init_set_ui(mpz_b, b);
    mpz_init_set_ui(mpz_p, p);
    mpz_init_set_ui(result, 1);

    mpzIndexCalculus(mpz_a, mpz_b, mpz_p, result);  
    // mpz_t at, mod;
    // mpz_init_set_ui(at, 1);
    // mpz_sub_ui(at,at,2);
    // gmp_printf("%Zd\n", at);
    // mpz_init_set_ui(mod, 1337);
    // mpz_mod(at, at, mod);
    // gmp_printf("%Zd\n", at);
    // mpz_t *tmp = malloc(10 * sizeof(mpz_t));
    // for (int i = 0; i < 10; i++){
    //     mpz_init_set_ui(tmp[i], i);
    // }
    // mpz_t tmp2;
    // mpz_init_set_ui(tmp2,43);
    // for (int i = 0; i < 10; i++) {
    //     gmp_printf("%i: %Zd\n", i, tmp[i]);
    // }
    // int indx = binary_search(tmp, 10, tmp2);
    // printf("index = %i\n",indx);

    uint64_t res = mpz_get_ui(result);
    mpz_clear(mpz_a);
    mpz_clear(mpz_b);
    mpz_clear(mpz_p);
    mpz_clear(result);
    return res;
}