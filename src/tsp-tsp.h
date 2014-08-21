#ifndef TSP_TSP_H
#define TSP_TSP_H
#include <stdint.h>

/* dernier minimum trouvé */
extern int minimum;

int present (int city, int hops, tsp_path_t path, uint64_t vpres);
void tsp (int hops, int len, uint64_t vpres, tsp_path_t path, long long int *cuts, tsp_path_t sol, int *sol_len);

#endif
