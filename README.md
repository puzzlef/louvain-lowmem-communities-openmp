Design of memory-efficient OpenMP-based Parallel [Louvain] algorithm for [community detection].

Community detection in graphs identifies groups of nodes that are more densely connected within the groups than between them. While many studies focus on enhancing detection performance, managing memory becomes critical for processing large graphs on shared-memory systems. Recently, we developed efficient implementations of the Louvain, Leiden, and Label Propagation Algorithms (LPA) for community detection, though these methods incur high memory costs due to collision-free per-thread hashtables. To mitigate this, we introduce memory-efficient alternatives based on weighted Misra-Gries (MG) sketches, which replace the per-thread hashtables and significantly reduce memory usage in Louvain, Leiden, and LPA implementations. This approach achieves only a minor quality reduction (up to 1%) with moderate runtime penalties. We believe these slightly slower but memory-efficient methods are well-suited for parallel processing and could offer better performance than current memory-intensive techniques on systems with numerous threads.

Below we plot the time taken by [Default Louvain], weighted Boyer-Moore (BM) based Louvain, and weighted Misra-Gries (MG) based Louvain with `k=8` slots, on 13 different graphs. MG8 Louvain is, on average, `1.48x` slower than Default Louvain, but requires significantly less memory (i.e., just `64 bytes` per thread).

[![](https://i.imgur.com/derv6de.png)][sheets-o1]

Next, we plot the modularity of communities identified by Default Louvain, BM Leiden, and MG8 Louvain. MG8 Louvain, on average, achieves a modularity within `1%` of Default Louvain - and is thus a viable memory-efficient alternative to the Default Louvain. Higher values of `k` can be chosen for even better modularity, if needed.

[![](https://i.imgur.com/EiM7Hdn.png)][sheets-o1]

Refer to our technical report for more details: \
[Memory-Efficient Community Detection on Large Graphs Using Weighted Sketches][report].

<br>

> [!NOTE]
> You can just copy `main.sh` to your system and run it. \
> For the code, refer to `main.cxx`.


[Louvain]: https://en.wikipedia.org/wiki/Louvain_method
[community detection]: https://en.wikipedia.org/wiki/Community_structure
[Default Louvain]: https://github.com/puzzlef/louvain-communities-openmp
[sheets-o1]: https://docs.google.com/spreadsheets/d/1fZhh2oPDB2nBKUbBowU-VxpeJ0cymWJHLrtkmP7aG5Y/edit?usp=sharing
[report]: https://arxiv.org/abs/2411.02268

<br>
<br>


### Code structure

The code structure is as follows:

```bash
- inc/_algorithm.hxx: Algorithm utility functions
- inc/_bitset.hxx: Bitset manipulation functions
- inc/_cmath.hxx: Math functions
- inc/_ctypes.hxx: Data type utility functions
- inc/_cuda.hxx: CUDA utility functions
- inc/_debug.hxx: Debugging macros (LOG, ASSERT, ...)
- inc/_iostream.hxx: Input/output stream functions
- inc/_iterator.hxx: Iterator utility functions
- inc/_main.hxx: Main program header
- inc/_mpi.hxx: MPI (Message Passing Interface) utility functions
- inc/_openmp.hxx: OpenMP utility functions
- inc/_queue.hxx: Queue utility functions
- inc/_random.hxx: Random number generation functions
- inc/_string.hxx: String utility functions
- inc/_utility.hxx: Runtime measurement functions
- inc/_vector.hxx: Vector utility functions
- inc/batch.hxx: Batch update generation functions
- inc/bfs.hxx: Breadth-first search algorithms
- inc/csr.hxx: Compressed Sparse Row (CSR) data structure functions
- inc/dfs.hxx: Depth-first search algorithms
- inc/duplicate.hxx: Graph duplicating functions
- inc/Graph.hxx: Graph data structure functions
- inc/louvain.hxx: Louvain community detection algorithm functions
- inc/louvainLowmem.hxx: Memory-efficient Louvain community detection algorithm functions
- inc/main.hxx: Main header
- inc/mtx.hxx: Graph file reading functions
- inc/properties.hxx: Graph Property functions
- inc/selfLoop.hxx: Graph Self-looping functions
- inc/symmetrize.hxx: Graph Symmetrization functions
- inc/transpose.hxx: Graph transpose functions
- inc/update.hxx: Update functions
- main.cxx: Experimentation code
- process.js: Node.js script for processing output logs
```

Note that each branch in this repository contains code for a specific experiment. The `main` branch contains code for the final experiment. If the intention of a branch in unclear, or if you have comments on our technical report, feel free to open an issue.

<br>
<br>


## References

- [Fast unfolding of communities in large networks; Vincent D. Blondel et al. (2008)](https://arxiv.org/abs/0803.0476)
- [Community Detection on the GPU; Md. Naim et al. (2017)](https://arxiv.org/abs/1305.2006)
- [Scalable Static and Dynamic Community Detection Using Grappolo; Mahantesh Halappanavar et al. (2017)](https://ieeexplore.ieee.org/document/8091047)
- [From Louvain to Leiden: guaranteeing well-connected communities; V.A. Traag et al. (2019)](https://www.nature.com/articles/s41598-019-41695-z)
- [CS224W: Machine Learning with Graphs | Louvain Algorithm; Jure Leskovec (2021)](https://www.youtube.com/watch?v=0zuiLBOIcsw)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)

<br>
<br>


[![](https://img.youtube.com/vi/M6npDdVGue4/maxresdefault.jpg)](https://www.youtube.com/watch?v=M6npDdVGue4)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
![](https://ga-beacon.deno.dev/G-KD28SG54JQ:hbAybl6nQFOtmVxW4if3xw/github.com/puzzlef/louvain-lowmem-communities-openmp)

[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
