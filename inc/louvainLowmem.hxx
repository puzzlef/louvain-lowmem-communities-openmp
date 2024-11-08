#pragma once
#include <utility>
#include <array>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "_main.hxx"
#include "Graph.hxx"
#include "properties.hxx"
#include "csr.hxx"
#include "louvain.hxx"

using std::array;
using std::vector;
using std::make_pair;
using std::swap;
using std::max;




#pragma region METHODS
#pragma region HASHTABLES
/**
 * Allocate a number of hashtables.
 * @param mcs majority communities vertex u is linked to (updated)
 * @param mws total edge weight from vertex u to community C (updated)
 */
template <class K, class V, size_t SLOTS>
inline void louvainLowmemAllocateHashtablesW(vector<array<K, SLOTS>*>& mcs, vector<array<V, SLOTS>*>& mws) {
  size_t N = mcs.size();
  for (size_t i=0; i<N; ++i) {
    mcs[i] = new array<K, SLOTS>();
    mws[i] = new array<V, SLOTS>();
  }
}


/**
 * Free a number of hashtables.
 * @param mcs majority communities vertex u is linked to (updated)
 * @param mws total edge weight from vertex u to community C (updated)
 */
template <class K, class V, size_t SLOTS>
inline void louvainLowmemFreeHashtablesW(vector<array<K, SLOTS>*>& mcs, vector<array<V, SLOTS>*>& mws) {
  size_t N = mcs.size();
  for (size_t i=0; i<N; ++i) {
    delete mcs[i];
    delete mws[i];
  }
}
#pragma endregion




#pragma region CHANGE COMMUNITY
/**
 * Scan an edge community connected to a vertex.
 * @param mcs majority communities vertex u is linked to (updated)
 * @param mws total edge weight from vertex u to community C (updated)
 * @param u given vertex
 * @param v outgoing edge vertex
 * @param w outgoing edge weight
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class K, class V, size_t SLOTS>
inline void louvainLowmemScanCommunityW(array<K, SLOTS>& mcs, array<V, SLOTS>& mws, K u, K v, V w, const vector<K>& vcom) {
  if (!SELF && u==v) return;
  K c = vcom[v];
  // Add edge weight to community.
  for (int i=0; i<SLOTS; ++i)
    mws[i] += mcs[i]==c? w : V();
  // Check if community is already in the list.
  int has = 0;
  for (int i=0; i<SLOTS; ++i)
    has |= mcs[i]==c? -1 : 0;
  if (has) return;
  // Find empty slot.
  int f = -1;
  for (int i=0; i<SLOTS; ++i)
    if (mws[i]==V()) f = i;
  // Add community to list.
  if (f>=0) {
    mcs[f] = c;
    mws[f] = w;
  }
  // Subtract edge weight from non-matching communities.
  else {
    for (int i=0; i<SLOTS; ++i)
      mws[i] = max(mws[i] - w, V());
  }
}


/**
 * Scan communities connected to a vertex.
 * @param mcs majority communities vertex u is linked to (updated)
 * @param mws total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 */
template <bool SELF=false, class G, class K, class V, size_t SLOTS>
inline void louvainLowmemScanCommunitiesW(array<K, SLOTS>& mcs, array<V, SLOTS>& mws, const G& x, K u, const vector<K>& vcom) {
  x.forEachEdge(u, [&](auto v, auto w) { louvainLowmemScanCommunityW<SELF>(mcs, mws, u, v, V(w), vcom); });
}


/**
 * Scan communities connected to a vertex.
 * @param x original graph
 * @param u given vertex
 * @param vcom community each vertex belongs to
 * @returns [majority community, total edge weight to community]
 */
template <bool SELF=false, class G, class K>
inline auto louvainLowmemScanCommunitiesMajorityW(const G& x, K u, const vector<K>& vcom) {
  using V = typename G::edge_value_type;
  K mc = K();
  V mw = V();
  x.forEachEdge(u, [&](auto v, auto w) {
    if (!SELF && u==v) return;
    K c = vcom[v];
    if (c==mc)     mw += w;
    else if (mw>w) mw -= w;
    else { mc = c; mw  = w; }
  });
  return make_pair(mc, mw);
}


/**
 * Clear communities scan data.
 * @param mws communities vertex u is linked to (updated)
 */
template <class V, size_t SLOTS>
inline void louvainLowmemClearScanW(array<V, SLOTS>& mws) {
  for (int i=0; i<SLOTS; ++i)
    mws[i] = V();
}


/**
 * Choose connected community with best delta modularity.
 * @param mws total edge weight from vertex u to community C (updated)
 * @param x original graph
 * @param u given vertex
 * @param mcs majority communities vertex u is linked to
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns [best community, delta modularity]
 */
template <bool SELF=false, class G, class K, class V, class W, size_t SLOTS>
inline auto louvainLowmemChooseCommunityW(array<V, SLOTS>& mws, const G& x, K u, const array<K, SLOTS>& mcs, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R) {
  K d  = vcom[u];
  V dw = V();
  louvainLowmemClearScanW(mws);
  // Compute total edge weight to communities.
  x.forEachEdge(u, [&](auto v, auto w) {
    if (!SELF && u==v) return;
    K c = vcom[v];
    if (c==d) dw += w;
    for (int i=0; i<SLOTS; ++i)
      mws[i] += mcs[i]==c? w : V();
  });
  // Choose community with best delta modularity.
  K cmax = K();
  W emax = W();
  for (int i=0; i<SLOTS; ++i) {
    K c = mcs[i];
    V w = mws[i];
    if (!w) continue;
    if (!SELF && c==d) continue;
    W e = deltaModularity(w, dw, vtot[u], ctot[c], ctot[d], M, R);
    if (e>emax) { emax = e; cmax = c; }
  }
  return make_pair(cmax, emax);
}


/**
 * Compute delta modularity to majority community.
 * @param x original graph
 * @param u given vertex
 * @param c majority community vertex u is linked to
 * @param vcom community each vertex belongs to
 * @param vtot total edge weight of each vertex
 * @param ctot total edge weight of each community
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @returns delta modularity to majority community
 * @note You need to ensure that c!=d
 */
template <bool SELF=false, class G, class K, class W>
inline W louvainLowmemDeltaModularityMajority(const G& x, K u, K c, const vector<K>& vcom, const vector<W>& vtot, const vector<W>& ctot, double M, double R) {
  using V = typename G::edge_value_type;
  K d  = vcom[u];
  V dw = V();
  V cw = V();
  x.forEachEdge(u, [&](auto v, auto w) {
    if (!SELF && u==v) return;
    K b = vcom[v];
    if (b==d) dw += w;
    if (b==c) cw += w;
  });
  return deltaModularity(cw, dw, vtot[u], ctot[c], ctot[d], M, R);
}
#pragma endregion




#pragma region LOCAL-MOVING PHASE
/**
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param mcs majority communities vertex u is linked to (temporary buffer, updated)
 * @param mws total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @param fa is vertex allowed to be updated?
 * @returns iterations performed (0 if converged already)
 */
template <bool MULTI=true, class G, class K, class V, class W, class B, size_t SLOTS, class FC, class FA>
inline int louvainLowmemMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<array<K, SLOTS>*>& mcs, vector<array<V, SLOTS>*>& mws, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc, FA fa) {
  size_t S = x.span();
  int l = 0;
  W  el = W();
  for (; l<L;) {
    el = W();
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:el)
    for (K u=0; u<S; ++u) {
      int t = omp_get_thread_num();
      if (!x.hasVertex(u)) continue;
      if (!fa(u) || !vaff[u]) continue;
      // Perform local-move using multiple majority communities (Misra-Gries sketch).
      while (MULTI) {
        louvainLowmemClearScanW(*mws[t]);
        louvainLowmemScanCommunitiesW(*mcs[t], *mws[t], x, u, vcom);
        auto [c, e] = louvainLowmemChooseCommunityW(*mws[t], x, u, *mcs[t], vcom, vtot, ctot, M, R);
        if (e<=0) break;
        louvainChangeCommunityOmpW(vcom, ctot, x, u, c, vtot);
        x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); });
        el += e;  // l1-norm
        break;
      }
      // Perform local-move using a single majority community (Boyer-Moore majority vote).
      while (!MULTI) {
        auto [c, w] = louvainLowmemScanCommunitiesMajorityW(x, u, vcom);
        if   (c==vcom[u]) break;
        auto  e = louvainLowmemDeltaModularityMajority(x, u, c, vcom, vtot, ctot, M, R);
        if (e<=0) break;
        louvainChangeCommunityOmpW(vcom, ctot, x, u, c, vtot);
        x.forEachEdgeKey(u, [&](auto v) { vaff[v] = B(1); });
        el += e;  // l1-norm
        break;
      }
      vaff[u] = B();
    }
    if (fc(el, l++)) break;
  }
  return l>1 || el? l : 0;
}


/**
 * Louvain algorithm's local moving phase.
 * @param vcom community each vertex belongs to (initial, updated)
 * @param ctot total edge weight of each community (precalculated, updated)
 * @param vaff is vertex affected flag (updated)
 * @param mcs majority communities vertex u is linked to (temporary buffer, updated)
 * @param mws total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vtot total edge weight of each vertex
 * @param M total weight of "undirected" graph (1/2 of directed graph)
 * @param R resolution (0, 1]
 * @param L max iterations
 * @param fc has local moving phase converged?
 * @returns iterations performed (0 if converged already)
 */
template <bool MULTI=true, class G, class K, class V, class W, class B, size_t SLOTS, class FC>
inline int louvainLowmemMoveOmpW(vector<K>& vcom, vector<W>& ctot, vector<B>& vaff, vector<array<K, SLOTS>*>& mcs, vector<array<V, SLOTS>*>& mws, const G& x, const vector<W>& vtot, double M, double R, int L, FC fc) {
  auto fa = [](auto u) { return true; };
  return louvainLowmemMoveOmpW<MULTI>(vcom, ctot, vaff, mcs, mws, x, vtot, M, R, L, fc, fa);
}
#pragma endregion




#pragma region AGGREGATION PHASE
/**
 * Aggregate outgoing edges of each community.
 * @param ydeg degree of each community (updated)
 * @param yedg vertex ids of outgoing edges of each community (updated)
 * @param ywei weights of outgoing edges of each community (updated)
 * @param mcs majority communities vertex u is linked to (temporary buffer, updated)
 * @param mws total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param coff offsets for vertices belonging to each community
 * @param cedg vertices belonging to each community
 * @param yoff offsets for vertices belonging to each community
 */
template <int CHUNK_SIZE=2048, class G, class K, class V, class W, size_t SLOTS>
inline void louvainLowmemAggregateEdgesOmpW(vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<array<K, SLOTS>*>& mcs, vector<array<V, SLOTS>*>& mws, const G& x, const vector<K>& vcom, const vector<K>& coff, const vector<K>& cedg, const vector<size_t>& yoff) {
  size_t C = coff.size() - 1;
  fillValueOmpU(ydeg, K());
  #pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
  for (K c=0; c<C; ++c) {
    int t = omp_get_thread_num();
    K   n = csrDegree(coff, c);
    if (n==0) continue;
    louvainLowmemClearScanW(*mws[t]);
    csrForEachEdgeKey(coff, cedg, c, [&](auto u) {
      louvainLowmemScanCommunitiesW<true>(*mcs[t], *mws[t], x, u, vcom);
    });
    for (int i=0; i<SLOTS; ++i) {
      K d = (*mcs[t])[i];
      V w = (*mws[t])[i];
      if (!w) continue;
      csrAddEdgeU(ydeg, yedg, ywei, yoff, c, d, W(w));
    }
  }
}


/**
 * Louvain algorithm's community aggregation phase.
 * @param yoff offsets for vertices belonging to each community (updated)
 * @param ydeg degree of each community (updated)
 * @param yedg vertex ids of outgoing edges of each community (updated)
 * @param ywei weights of outgoing edges of each community (updated)
 * @param bufs buffer for exclusive scan of size |threads| (scratch)
 * @param mcs majority communities vertex u is linked to (temporary buffer, updated)
 * @param mws total edge weight from vertex u to community C (temporary buffer, updated)
 * @param x original graph
 * @param vcom community each vertex belongs to
 * @param coff offsets for vertices belonging to each community
 * @param cedg vertices belonging to each community
 */
template <int CHUNK_SIZE=2048, class G, class K, class V, class W, size_t SLOTS>
inline void louvainLowmemAggregateOmpW(vector<size_t>& yoff, vector<K>& ydeg, vector<K>& yedg, vector<W>& ywei, vector<size_t>& bufs, vector<array<K, SLOTS>*>& mcs, vector<array<V, SLOTS>*>& mws, const G& x, const vector<K>& vcom, vector<K>& coff, vector<K>& cedg) {
  size_t C = coff.size() - 1;
  louvainCommunityTotalDegreeOmpW(yoff, x, vcom);
  yoff[C] = exclusiveScanOmpW(yoff.data(), bufs.data(), yoff.data(), C);
  louvainLowmemAggregateEdgesOmpW<CHUNK_SIZE>(ydeg, yedg, ywei, mcs, mws, x, vcom, coff, cedg, yoff);
}
#pragma endregion




#pragma region ENVIRONMENT SETUP
/**
 * Setup and perform the Louvain algorithm.
 * @param x original graph
 * @param o louvain options
 * @param fi initializing community membership and total vertex/community weights (vcom, vtot, ctot)
 * @param fm marking affected vertices (vaff, mcs, mws, vcom, vtot, ctot)
 * @param fa is vertex allowed to be updated? (u)
 * @returns louvain result
 */
template <bool MULTI=true, size_t SLOTS=8, class G, class FI, class FM, class FA>
inline auto louvainLowmemInvokeOmp(const G& x, const LouvainOptions& o, FI fi, FM fm, FA fa) {
  using  K = typename G::key_type;
  using  V = typename G::edge_value_type;
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = char;
  // Options.
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  // Get graph properties.
  size_t X = x.size();
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
  // Allocate buffers.
  int    T = omp_get_max_threads();
  vector<B> vaff(S);            // Affected vertex flag (any pass)
  vector<K> ucom(S), vcom(S);   // Community membership (first pass, current pass)
  vector<W> utot(S), vtot(S);   // Total vertex weights (first pass, current pass)
  vector<W> ctot(S);            // Total community weights (any pass)
  vector<K> bufk(T);            // Buffer for exclusive scan
  vector<size_t> bufs(T);       // Buffer for exclusive scan
  vector<array<K, SLOTS>*> mcs(T);  // Hashtable keys
  vector<array<V, SLOTS>*> mws(T);  // Hashtable values
  louvainLowmemAllocateHashtablesW(mcs, mws);
  size_t Z = max(size_t(o.aggregationTolerance * X), X);
  size_t Y = max(size_t(o.aggregationTolerance * Z), Z);
  DiGraphCsr<K, None, None, K> cv(S, S);  // CSR for community vertices
  DiGraphCsr<K, None, W> y(S, Y);         // CSR for aggregated graph (input);  y(S, X)
  DiGraphCsr<K, None, W> z(S, Z);         // CSR for aggregated graph (output); z(S, X)
  // Perform Louvain algorithm.
  float tm = 0, ti = 0, tp = 0, tl = 0, ta = 0;  // Time spent in different phases
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el<=E; };
    // Reset buffers, in case of multiple runs.
    fillValueOmpU(vaff, B());
    fillValueOmpU(ucom, K());
    fillValueOmpU(vcom, K());
    fillValueOmpU(utot, W());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    cv.respan(S);
    y .respan(S);
    z .respan(S);
    // Time the algorithm.
    mark([&]() {
      // Initialize community membership and total vertex/community weights.
      ti += measureDuration([&]() { fi(ucom, utot, ctot); });
      // Mark affected vertices.
      tm += measureDuration([&]() { fm(vaff, mcs, mws, ucom, utot, ctot); });
      // Start timing first pass.
      auto t0 = timeNow(), t1 = t0;
      // Start local-moving, refinement, aggregation phases.
      // NOTE: In first pass, the input graph is a DiGraph.
      // NOTE: For subsequent passes, the input graph is a DiGraphCsr (optimization).
      for (l=0, p=0; M>0 && P>0;) {
        if (p==1) t1 = timeNow();
        bool isFirst = p==0;
        int m = 0;
        tl += measureDuration([&]() {
          if (isFirst) m += louvainLowmemMoveOmpW<MULTI>(ucom, ctot, vaff, mcs, mws, x, utot, M, R, L, fc, fa);
          else         m += louvainLowmemMoveOmpW<MULTI>(vcom, ctot, vaff, mcs, mws, y, vtot, M, R, L, fc);
        });
        l += max(m, 1); ++p;
        if (m<=1 || p>=P) break;
        size_t GN = isFirst? x.order() : y.order();
        size_t CN = 0;
        if (isFirst) CN = louvainCommunityExistsOmpW(cv.degrees, x, ucom);
        else         CN = louvainCommunityExistsOmpW(cv.degrees, y, vcom);
        if (double(CN)/GN >= o.aggregationTolerance) break;
        if (isFirst) louvainRenumberCommunitiesOmpW(ucom, cv.degrees, bufk, x);
        else         louvainRenumberCommunitiesOmpW(vcom, cv.degrees, bufk, y);
        if (isFirst) {}
        else         louvainLookupCommunitiesOmpU(ucom, vcom);
        ta += measureDuration([&]() {
          cv.respan(CN); z.respan(CN);
          if (isFirst) louvainCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, x, ucom);
          else         louvainCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, y, vcom);
          if (isFirst) louvainLowmemAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, mcs, mws, x, ucom, cv.offsets, cv.edgeKeys);
          else         louvainLowmemAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, mcs, mws, y, vcom, cv.offsets, cv.edgeKeys);
        });
        swap(y, z);
        // fillValueOmpU(vcom.data(), CN, K());
        // fillValueOmpU(ctot.data(), CN, W());
        fillValueOmpU(vtot.data(), CN, W());
        fillValueOmpU(vaff.data(), CN, B(1));
        louvainVertexWeightsOmpW(vtot, y);
        louvainInitializeOmpW(vcom, ctot, y, vtot);
        E /= o.toleranceDrop;
      }
      if (p<=1) {}
      else      louvainLookupCommunitiesOmpU(ucom, vcom);
      if (p<=1) t1 = timeNow();
      tp += duration(t0, t1);
    });
  }, o.repeat);
  louvainLowmemFreeHashtablesW(mcs, mws);
  return LouvainResult<K>(ucom, utot, ctot, l, p, t, tm/o.repeat, ti/o.repeat, tp/o.repeat, tl/o.repeat, ta/o.repeat);
}
#pragma endregion




#pragma region STATIC APPROACH
/**
 * Obtain the community membership of each vertex with Static Louvain.
 * @param x original graph
 * @param o louvain options
 * @returns louvain result
 */
template <bool MULTI=true, size_t SLOTS=8, class G>
inline auto louvainLowmemStaticOmp(const G& x, const LouvainOptions& o={}) {
  using B = char;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    louvainVertexWeightsOmpW(vtot, x);
    louvainInitializeOmpW(vcom, ctot, x, vtot);
  };
  auto fm = [ ](auto& vaff, auto& mcs, auto& mws, const auto& vcom, const auto& vtot, const auto& ctot) {
    fillValueOmpU(vaff, B(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainLowmemInvokeOmp<MULTI, SLOTS>(x, o, fi, fm, fa);
}
#pragma endregion
#pragma endregion
