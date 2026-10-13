# miga
A matlab implementation of isogeometric analysis based on sum factorization

## Installation

1. Clone the repository and initialize the submodules[^1] by using the shell command
```bash
 git clone --recurse-submodules git@github.com:AndreaBressan/miga.git REPOSITORY_PATH
```
2. Start matlab and add the repository to its path: in the command window type
```matlab
 addpath(genpath(REPOSITORY_PATH));
```
3. Optionally compile the provided mesh file to speed-up sparse matrix forming[^2]: in matlab command window type
```matlab
 cd REPOSITORY_PATH
 mex_build_all
```

[^1]: `miga` uses [chebfun](https://github.com/chebfun/chebfun) to compute quadrature weights and nodes.
[^2]: The non-zeros computed by the sum-factorization algorithm are positioned in a known order so that sorting them requires order ($O(\textrm{nnz})$) operations instead of order ($O(\textrm{nnz} \log(\textrm{nnz}))$) required by calling `sparse(row,col,val,n,m)`.
