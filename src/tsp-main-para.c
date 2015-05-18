#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <unistd.h>

#include "tsp-types.h"
#include "tsp-job.h"
#include "tsp-genmap.h"
#include "tsp-print.h"
#include "tsp-tsp.h"
#include "tsp-lp.h"
#include "tsp-hkbound.h"
#include <pthread.h>

/* clock_gettime definition for OSX compatibility */
#ifdef __MACH__

#include <sys/time.h>

#define CLOCK_REALTIME 0

int clock_gettime(int clk_id, struct timespec* t) {
    struct timeval now;
    int rv = gettimeofday(&now, NULL);
    if (rv) return rv;
    t->tv_sec  = now.tv_sec;
    t->tv_nsec = now.tv_usec * 1000;
    return 0;
}

#endif

/* macro de mesure de temps, retourne une valeur en nanosecondes */
#define TIME_DIFF(t1, t2) \
  ((t2.tv_sec - t1.tv_sec) * 1000000000ll + (long long int) (t2.tv_nsec - t1.tv_nsec))


/* tableau des distances */
tsp_distance_matrix_t tsp_distance ={};

/** Paramètres **/

/* nombre de villes */
int nb_towns=10;
/* graine */
long int myseed= 0;
/* nombre de threads */
int nb_threads=1;

/* affichage SVG */
bool affiche_sol= false;
bool affiche_progress=false;
bool quiet=false;

static void generate_tsp_jobs (struct tsp_queue *q, int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len, int depth)
{
    if (len >= minimum) {
        (*cuts)++ ;
        return;
    }
    
    if (hops == depth) {
        /* On enregistre du travail à faire plus tard... */
      add_job (q, path, hops, len, vpres);
    } else {
        int me = path [hops - 1];        
        for (int i = 0; i < nb_towns; i++) {
	  if (!present (i, hops, path, vpres)) {
                path[hops] = i;
		vpres |= (1<<i);
                int dist = tsp_distance[me][i];
                generate_tsp_jobs (q, hops + 1, len + dist, vpres, path, cuts, sol, sol_len, depth);
		vpres &= (~(1<<i));
            }
        }
    }
}

static void usage(const char *name) {
  fprintf (stderr, "Usage: %s [-s] <ncities> <seed> <nthreads>\n", name);
  exit (-1);
}
struct param_job{
    tsp_path_t path;
    tsp_path_t solution;
    uint64_t vpres; // permet de reconstruire le chemin partiel <> tableau de bits
    tsp_path_t best_sol;   // meilleure solution connue
    int sol_len;  // longueur de la solution partiel
    long long int cuts; //
    struct tsp_queue q; // queue


};
void do_parallel(struct param_job* paramJob){
    while (!empty_queue (&(paramJob->q))) {
        int hops = 0, len = 0;

        pthread_mutex_lock(&mutex_jobs);
        get_job (&q, solution, &hops, &len, &(paramJob->vpres);
        pthread_mutex_unlock(&mutex_jobs);


        // le noeud est moins bon que la solution courante
        if (minimum < INT_MAX
            && (nb_towns - hops) > 10
            && ( (lower_bound_using_hk(paramJob->solution, hops, len, paramJob->vpres)) >= minimum
                 || (lower_bound_using_lp(paramJob->solution, hops, len, paramJob->vpres)) >= minimum)
                )

            continue;

        tsp (hops, len, paramJob->vpres, paramJob->solution, &(paramJob->cuts), paramJob->best_sol, &(paramJob->sol_len));
    }
}
int main (int argc, char **argv)
{
    unsigned long long perf;
    tsp_path_t path;
    uint64_t vpres=0; // permet de reconstruire le chemin partiel <> tableau de bits
    tsp_path_t best_sol;   // meilleure solution connue
    int sol_len;  // longueur de la solution partiel
    long long int cuts = 0; //
    struct tsp_queue q; // queue
    struct timespec t1, t2;




    /* lire les arguments */
    int opt;
    while ((opt = getopt(argc, argv, "spq")) != -1) {
      switch (opt) {
      case 's':
	affiche_sol = true;
	break;
      case 'p':
	affiche_progress = true;
	break;
      case 'q':
	quiet = true;
	break;
      default:
	usage(argv[0]);
	break;
      }
    }

    if (optind != argc-3)
      usage(argv[0]);

    nb_towns = atoi(argv[optind]);
    myseed = atol(argv[optind+1]);
    nb_threads = atoi(argv[optind+2]);
    assert(nb_towns > 0);
    assert(nb_threads > 0);
   
    minimum = INT_MAX;
      
    /* generer la carte et la matrice de distance */
    if (! quiet)
      fprintf (stderr, "ncities = %3d\n", nb_towns);
    genmap ();

    init_queue (&q);

    clock_gettime (CLOCK_REALTIME, &t1);

    memset (path, -1, MAX_TOWNS * sizeof (int));
    path[0] = 0;
    vpres=1;

    /* mettre les travaux dans la file d'attente */
    generate_tsp_jobs (&q, 1, 0, vpres, path, &cuts, best_sol, & sol_len, 3);
    no_more_jobs (&q);


    /* mutex  */
    pthread_mutex_t mutex_jobs;//,mutex_minimum, mutex_variables_globales;
    pthread_mutex_t mutex_sol_len;
    pthread_mutex_t mutex_cuts;
    pthread_mutex_t mutex_best_sol;
    pthread_mutex_t mutex_solution;

    pthread_mutex_init(&mutex_jobs,NULL);
    pthread_mutex_init(&mutex_sol_len,NULL);
    pthread_mutex_init(&mutex_best_sol,NULL);
    pthread_mutex_init(&mutex_cuts,NULL);
    pthread_mutex_init(&mutex_solution,NULL);

    pthread_t threads[nb_threads];
    int ret_threads[nb_threads];
    struct param_job paramJob;
    paramJob.vpres=0;
    paramJob.cuts =0;

    /* calculer chacun des travaux */
    tsp_path_t solution;
    memset (solution, -1, MAX_TOWNS * sizeof (int));
    solution[0] = 0;
    (*paramJob.solution) = solution;

    // Création des threads
    for (int i = 0; i < nb_threads; ++i) {
        ret_threads[i] = pthread_create(&threads[i],NULL,do_parallel,&paramJob);
    }


    while (!empty_queue (&q)) {
        int hops = 0, len = 0;

        pthread_mutex_lock(&mutex_jobs);
        get_job (&q, solution, &hops, &len, &vpres);
        pthread_mutex_unlock(&mutex_jobs);


        // le noeud est moins bon que la solution courante
        if (minimum < INT_MAX
            && (nb_towns - hops) > 10
            && ( (lower_bound_using_hk(solution, hops, len, vpres)) >= minimum
                 || (lower_bound_using_lp(solution, hops, len, vpres)) >= minimum)
                )

            continue;

        tsp (hops, len, vpres, solution, &cuts, best_sol, &sol_len);
    }
// to arrete parell
    clock_gettime (CLOCK_REALTIME, &t2);

    if (affiche_sol)
        print_solution_svg (best_sol, sol_len);

    perf = TIME_DIFF (t1,t2);
    printf("<!-- # = %d seed = %ld len = %d threads = %d time = %lld.%03lld ms ( %lld coupures ) -->\n",
           nb_towns, myseed, sol_len, nb_threads,
           perf/1000000ll, perf%1000000ll, cuts);

    return 0 ;
}

