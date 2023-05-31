pcg3: 
- instrumented PCG code to record the data to be analyzed later (this is the code used in the paper)

pcg4: 
- instrumented PCG code to also compute first and second temporal gradient (this is the new idea to be explore in this project)

Step 2:
- matrices: compute different norms of matrix A, and incomplete cholesky fractorization for PCG solvers (saved in matrices folder);

Step 3:
- solving: inject error in pcg3 and record information
- matrices: call solving for different matrices
- Step 2 must be executed before Step 3

Step5:
- analysis: computes overhead and slowdown of our algorithm and random algorithm with different protection fractions. 
- Step 3 must be executed before Step 5

Step 6:
- impact: plot Fig 1 in the paper
- correlate: plot Fig 2 in the paper; must run Step 2 & 3 first
- overhead: plot Fig 6, 7, 8, in the paper; must run Step 2 & 3 & 5 first
- plotall: plot Fig 9 in the paper; must run everything first

