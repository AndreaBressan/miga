(#) MIGA — Installation and quick start

MIGA stands "for Matlab IsoGeometric Analysis" and it is an implementation of isogeometric analysis implementing using state of the art integration and solving algorithms.

## Prerequisites

- MATLAB (R2020b or newer required). Some compiled components require the interleaved-complex API available since R2018a
- the `git` version system is required to collaborate to the library development
- a working C/C++ toolchain configured for MATLAB's `mex` is required to achieve better performances by using the provided c++ implementation

## Obtaining the library

Download the library from `https://github.com/AndreaBressan/miga.git` using either your browser or the `git` program

```bash
git clone https://github.com/AndreaBressan/miga.git
cd miga
```

## Add MIGA to the MATLAB path

Next add the root of folder of miga to the MATLAB path using either the graphical interface or the following commands (where '/full/path/to/miga' depends on your installation path) 

```matlab
addpath(genpath('/full/path/to/miga/miga'))
% optionally make the change persistent
% savepath
```
## Optionally compile mex files

Set up the mex compiler in MATLAB, i.e., choose a compiler and then compile the mex files

```matlab
mex -setup
mex_build_all
```

At the moment the mex files are 
- `miga/tensor/computation/tSparsityCPP/tSparsityToSparse.cpp` that speeds up the matrix assembly.

In case of failures check the source files as they contain additional documentation.

## Quick usage examples

Once the path is set you can use the library from MATLAB. Try running the examples in the PDEs folder.

## Run the built-in tests

The repository contains a test runner in `miga/utils/runMigaTests.m`. From MATLAB run:

```
runMigaTests
```

This will execute the project's unit / integration tests and report results. 

## Troubleshooting

- mex: "No supported compiler found"
	- Run `mex -setup` in MATLAB and install a supported compiler on your OS (on Linux, `g++` usually works).

- mex compile errors about C++ standard or missing symbols
	- Try adding the flags suggested in the `.cpp` file: `mex -R2018a CXXFLAGS='$CXXFLAGS -std=c++20' filename.cpp`.
	- Verify your `g++` version supports C++20 or adjust the standard flag to a supported version.

- Runtime errors related to interleaved complex arrays or assertions about `MX_HAS_INTERLEAVED_COMPLEX`
	- Use MATLAB R2018a or newer and add `-R2018a` to the mex command.
