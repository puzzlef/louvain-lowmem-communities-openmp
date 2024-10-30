#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <iterator>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include "inc/main.hxx"

using namespace std;




#pragma region CONFIGURATION
#ifndef TYPE
/** Type of edge weights. */
#define TYPE float
#endif
#ifndef MAX_THREADS
/** Maximum number of threads to use. */
#define MAX_THREADS 64
#endif
#ifndef REPEAT_METHOD
/** Number of times to repeat each method. */
#define REPEAT_METHOD 5
#endif
#pragma endregion




#pragma region METHODS
#pragma region HELPERS
/**
 * Obtain the modularity of community structure on a graph.
 * @param x original graph
 * @param a rak result
 * @param M sum of edge weights
 * @returns modularity
 */
template <class G, class K>
inline double getModularity(const G& x, const LouvainResult<K>& a, double M) {
  auto fc = [&](auto u) { return a.membership[u]; };
  return modularityBy(x, fc, M, 1.0);
}
#pragma endregion




#pragma region PERFORM EXPERIMENT
/**
 * Perform the experiment.
 * @param x original graph
 */
template <class G>
void runExperiment(const G& x) {
  int repeat = REPEAT_METHOD;
  double   M = edgeWeightOmp(x)/2;
  // Follow a specific result logging format, which can be easily parsed later.
  auto flog = [&](const auto& ans, const char *technique, size_t numSlots=0) {
    printf(
      "{%03d threads} -> "
      "{%09.1fms, %09.1fms mark, %09.1fms init, %09.1fms first, %09.1fms move, %09.1fms aggr, %.3e slots, %04d iters, %04d passes, %01.9f modularity} %s\n",
      MAX_THREADS,
      ans.time, ans.markingTime, ans.initializationTime, ans.firstPassTime, ans.localMoveTime, ans.aggregationTime,
      double(numSlots), ans.iterations, ans.passes, getModularity(x, ans, M), technique
    );
  };
  // Find static Louvain.
  {
    auto b0 = louvainStaticOmp(x, {repeat});
    flog(b0, "louvainStaticOmp");
  }
  // Get community memberships on original graph (low memory).
  {
    auto b1 = louvainLowmemStaticOmp<false>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajority", 0);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 4>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 4);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 8>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 8);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 16>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 16);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 32>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 32);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 64>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 64);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 128>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 128);
  }
  {
    auto b1 = louvainLowmemStaticOmp<true, 256>(x, {repeat});
    flog(b1, "louvainLowmemStaticOmpMajorities", 256);
  }
}


/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns zero on success, non-zero on failure
 */
int main(int argc, char **argv) {
  using K = uint32_t;
  using V = TYPE;
  install_sigsegv();
  char *file     = argv[1];
  bool symmetric = argc>2? stoi(argv[2]) : false;
  bool weighted  = argc>3? stoi(argv[3]) : false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  DiGraph<K, None, V> x;
  readMtxOmpW(x, file, weighted); LOG(""); println(x);
  if (!symmetric) { symmetrizeOmpU(x); LOG(""); print(x); printf(" (symmetrize)\n"); }
  runExperiment(x);
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
