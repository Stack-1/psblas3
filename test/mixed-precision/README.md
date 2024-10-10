# Mixed Precision Conjugate Gradient Algorithm
This directory holds the informations about a new possible functionality to implement in the library. The main goal is always to solve Ax = b, using a mixed precision implementation of CG algorithm.
## Algorithm description
The main idea is to use both single precision and double precision matrix to save on computation speed for a subset of Conjugate Gradient operations. 
Steps:
- Read matrix A from file in double precision
- Read matrix A from file in single precision 
- Read vector b from file in double precision
- Read vector b from file in single precision
- Apply Conjugate Gradient
- Save results on file
## Error Analysys
## Getting Started
As a first thing you have to provide a file to specify the parameters to use in a single experiment:
The file should be composed as follows:
1) Number of parameters
2) Matrix file name (The file should be copied in the runs/ directory) [MM - ]
3) RHS vector [NONE - ]
4) Krilov methos [CG - ]
5) cxd
6) The format in which the matrix will be stored in memory [CSR - COO]
7) The type of partition to use to divide the data [BLOCK - -] 
8) ...
9) ...
10) ...
11) ...
12) ...
13) ...
14) ...
15) ...
16) ...
17) ...

To try this code you can use a specs.txt file containing the standard parameters, the only other thing you have to do is to download the matrix and the rhs vector, use the following links:
Matrix link: https://math.nist.gov/pub/MatrixMarket2/SPARSKIT/fidap/fidap019.mtx.gz
RHS vector link: https://math.nist.gov/pub/MatrixMarket2/SPARSKIT/fidap/fidap019_rhs1.mtx.gz
