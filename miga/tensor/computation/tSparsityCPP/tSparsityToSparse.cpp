/*=================================================================
* tSparsityToSparse
* This program populates a sparse matrix from the non-zeros and the array of
* univariate sparsities.
*
* Compile with mex -R2018a CXXFLAGS='$CXXFLAGS -std=c++20' tSparsityToSparse.cpp.
*
* This file contains dimension specific implementations for dimensions 1 to 7
* and a dimension generic implementation for dimensions >=8. 
* The function is chosen at run-time.
*
*
* This is a MEX-file for MATLAB.
* Copyright 2025 Andrea Bressan, IMATI-CNR
* All rights reserved.
* 
*=================================================================*/

#include "mex.h"
#include "matrix.h"

#include <math.h>
#include <span>
#include <vector>
#include <chrono>
#include <numeric>
#include <thread>
#include <algorithm>

constexpr  size_t minThreadWorkSize=1024 * 8;

/*
*   THREAD PARALLEL IMPLEMENTATION
*/
template <typename ...T>
void splitWork(size_t totalWorkSize, size_t numThreads, size_t minWorkSize, T...args);
size_t getNumThreads(int nrhs, const mxArray *prhs[]);

/*
*   TIME KEEPING INTERFACE
*   the time log is a vector of n-pairs (time instant, label)
*   the time report transforms it in a Matlab struct with fields
*   - time: the (n-1) vector of durations, i.e. differences of instants
*   - label: the last (n-1) labels
*/
typedef std::chrono::time_point<std::chrono::high_resolution_clock> timePoint;
typedef std::pair<timePoint, std::string> logEntry;
void addPartialTime(std::vector<logEntry> &time_log, std::string name);
mxArray* constructTimeReport(const std::vector<logEntry>& time_log);

/*
*   MXARRAY TO STD:SPAN CONVERSION (so that the reordering routines are not linked to Matlab)
*/
template <typename T>
std::span<T> getSpan(mxArray *arr);

/*
*   HELPER FUNCTIONS
*/
int checkInputs(int nrhs, const mxArray *prhs[]);
template <typename T>
void copyDataIntoSparse(
    std::span<T> nnzDest,
    std::span<mwIndex> rowDest,
    const std::span<T> nnzSrc,
    const std::vector<std::span<mxDouble>> &rowSrc,
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumRow,
    const std::span<mwIndex> cumCol,
    const std::vector<std::vector<mwIndex>> &colStart,
    const std::span<mwIndex> colStartTot,
    size_t numThreads);
void computeColStart(
    const std::vector<std::vector<mwIndex>> &colStart, 
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumCol,
    std::span<mwIndex> dest,
    size_t numThreads);

/*
* MAIN MEX FUNCTION
*
* in arguments are 
* -sparsityData: a cell array of size 3 x dim containing the following
* univariate sparsity data
*   -on the first  row the row-positions of the non-zeros
*   -on the second row the col-positions of the non-zeros
*   -on the third  row the number of non-zeros per column
* -shapeData: a dim x 2 double array containing in each row the univariate
*   dimensions
* -nnzData: a double array containing the values of the non-zeros
*/
void mexFunction(
    int nlhs,       mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    // init multi-threading and time keeping
    size_t numThreads=getNumThreads(nrhs, prhs);
    std::vector<logEntry> time_log;
    addPartialTime(time_log, "start");
    // Check inputs and outputs
    int dim=checkInputs(nrhs,prhs);
    if (nlhs > 2) mexErrMsgIdAndTxt("tSparsityToSparse:nrhs", "Wrong number of output arguments. Check the documentation.");
    if (nlhs == 0) return;
    addPartialTime(time_log, "checked input");
    // prepare the computation data
    mxArray * sparsityData = (mxArray *)prhs[0];
    mxArray * shapeData    = (mxArray *)prhs[1];
    mxArray * nnzData      = (mxArray *)prhs[2];
    mxArray *&sparseMat    = plhs[0];
    mwSize numNNZ = mxGetNumberOfElements(nnzData);
    // collect univariate data in std::spans
    std::vector<std::span<mxDouble>> rows(dim), cols(dim), nnzPerCol(dim);
    std::vector<mwIndex> cumNNZ(dim + 1);
    std::vector<mwIndex> cumRow(dim + 1);
    std::vector<mwIndex> cumCol(dim + 1);
    std::vector<std::vector<mwIndex>> colStart(dim);
    cumRow[0] = 1;
    cumCol[0] = 1;
    cumNNZ[0] = 1;

    auto shape=getSpan<mxDouble>(shapeData);

    for (int i = 0; i < dim; ++i)
    {
        rows[i]      = getSpan<mxDouble>(mxGetCell(sparsityData, 3*i));
        cols[i]      = getSpan<mxDouble>(mxGetCell(sparsityData, 3*i+1));
        nnzPerCol[i] = getSpan<mxDouble>(mxGetCell(sparsityData, 3*i+2));
        cumRow[i + 1] = cumRow[i] * mwIndex(shape[i]);
        cumCol[i + 1] = cumCol[i] * mwIndex(shape[dim+i]);
        cumNNZ[i + 1] = cumNNZ[i] * rows[i].size();
        colStart[i] = std::vector<mwIndex>(shape[dim+i] + 1);
        colStart[i][0] = 0;
        std::partial_sum(nnzPerCol[i].begin(), nnzPerCol[i].end(), colStart[i].begin() + 1);
    }
    bool    isComplex = mxIsComplex(nnzData);
    mwIndex numRows   = cumRow[dim];
    mwIndex numCols   = cumCol[dim];
    addPartialTime(time_log, "extract univariate data");

    // construct output
    sparseMat    = mxCreateSparse(numRows, numCols, numNNZ, isComplex ? mxCOMPLEX : mxREAL);
    mwIndex *IrB = mxGetIr(sparseMat);
    mwIndex *JcB = mxGetJc(sparseMat);
    auto colStartTot = std::span<mwIndex>(JcB, JcB + cumCol.back() + 1);
    std::span<mwIndex> rowDest(IrB, IrB + numNNZ);
    addPartialTime (time_log, "construct output matrix");
    computeColStart(colStart, cumNNZ, cumCol, colStartTot, numThreads);
    addPartialTime (time_log, "compute the column start positions");
    if (!isComplex)
    {
        std::span<mxDouble> nnzSrc =  getSpan<mxDouble>(nnzData);
        std::span<mxDouble> nnzDest = getSpan<mxDouble>(sparseMat);
        copyDataIntoSparse<mxDouble>(nnzDest, rowDest, nnzSrc,
                                     rows, cumNNZ, cumRow, cumCol, colStart, colStartTot,
                                     numThreads);
    }
    else
    {
        static_assert(MX_HAS_INTERLEAVED_COMPLEX, "tSparsityToSparse requires interleaved complex and does not work with old Matlab versions");
        std::span<mxComplexDouble> nnzSrc  = getSpan<mxComplexDouble>(nnzData);
        std::span<mxComplexDouble> nnzDest = getSpan<mxComplexDouble>(sparseMat);
        copyDataIntoSparse<mxComplexDouble>(nnzDest, rowDest, nnzSrc,
                                            rows, cumNNZ, cumRow, cumCol, colStart, colStartTot,
                                            numThreads);
    }
    addPartialTime(time_log, "copy non zeros");
    if (nlhs == 2) plhs[1]=constructTimeReport(time_log);
}

template <typename T, int dim>
void copyDataIntoSparseDimSpecific(
    std::span<T> nnzDest,
    std::span<mwIndex> rowDest,
    const std::span<T> nnzSrc,
    const std::vector<std::span<mxDouble>> &rowSrc,
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumRow,
    const std::span<mwIndex> cumCol,
    const std::vector<std::vector<mwIndex>> &colStart,
    const std::span<mwIndex> colStartTot,
    size_t start,
    size_t end)
{
    mwIndex col, colMI[dim], colRow[dim];
    mwIndex row, rowMI[dim];
    mwIndex srcPos = 0;
    mwIndex colEnd = 0;
    colRow[0] = 1;

    // init data for the middle of a column
    col = std::upper_bound(colStartTot.begin(), colStartTot.end(), start) - colStartTot.begin() - 1;
    row = start - colStartTot[col];
    srcPos = 0;
    for (int i = 0; i < dim; ++i)
    {
        colMI[i] = (col % cumCol[i + 1]) / cumCol[i];
        colRow[i] = (colStart[i][colMI[i] + 1] - colStart[i][colMI[i]]);
        rowMI[i] = row % colRow[i];
        row /= colRow[i];
        srcPos += (colStart[i][colMI[i]] + rowMI[i]) * cumNNZ[i];
    }
    colEnd = std::min(colStartTot[col + 1], end);
    while (start < end)
    {
        // copy non zeros of current column and compute their row position
        for (; start < colEnd; ++start)
        {
            nnzDest[start] = nnzSrc[srcPos];
            rowDest[start] = 0;
            int carryon = 1;
            for (int i = 0; i < dim; ++i)
            {
                rowDest[start] += (rowSrc[i][colStart[i][colMI[i]] + rowMI[i]] - 1) * cumRow[i];
                rowMI[i] = (rowMI[i] + carryon) % colRow[i];
                srcPos += carryon * cumNNZ[i]; // rowMI[i] has been incremented
                carryon = carryon * (rowMI[i] == 0 ? 1 : 0);
                srcPos -= carryon * colRow[i] * cumNNZ[i]; // rowMI[i] reached 0
            }
        }
        // set up data for next column
        ++col;
        srcPos = 0;
        for (int i = 0; i < dim; ++i)
        {
            colMI[i] = (col % cumCol[i + 1]) / cumCol[i];
            colRow[i] = (colStart[i][colMI[i] + 1] - colStart[i][colMI[i]]);
            rowMI[i] = 0;
            srcPos += colStart[i][colMI[i]] * cumNNZ[i];
        }
        colEnd = std::min(colStartTot[col + 1], end);
    }
}

template <typename T>
void copyDataIntoSparseDimGeneric(
    std::span<T> nnzDest,
    std::span<mwIndex> rowDest,
    const std::span<T> nnzSrc,
    const std::vector<std::span<mxDouble>> &rowSrc,
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumRow,
    const std::span<mwIndex> cumCol,
    const std::vector<std::vector<mwIndex>> &colStart,
    const std::span<mwIndex> colStartTot,
    size_t start,
    size_t end)
{
    int dim = rowSrc.size();
    std::vector<mwIndex> colMI(dim), colRow(dim), rowMI(dim);
    mwIndex col;
    mwIndex row;
    mwIndex srcPos = 0;
    mwIndex colEnd = 0;
    colRow[0] = 1;

    // init data for the middle of a column
    col = std::upper_bound(colStartTot.begin(), colStartTot.end(), start) - colStartTot.begin() - 1;
    row = start - colStartTot[col];
    srcPos = 0;
    for (int i = 0; i < dim; ++i)
    {
        colMI[i] = (col % cumCol[i + 1]) / cumCol[i];
        colRow[i] = (colStart[i][colMI[i] + 1] - colStart[i][colMI[i]]);
        rowMI[i] = row % colRow[i];
        row /= colRow[i];
        srcPos += (colStart[i][colMI[i]] + rowMI[i]) * cumNNZ[i];
    }
    colEnd = std::min(colStartTot[col + 1], end);
    while (start < end)
    {
        // copy non zeros of current column and compute their row position
        for (; start < colEnd; ++start)
        {
            nnzDest[start] = nnzSrc[srcPos];
            rowDest[start] = 0;
            int carryon = 1;
            for (int i = 0; i < dim; ++i)
            {
                rowDest[start] += (rowSrc[i][colStart[i][colMI[i]] + rowMI[i]] - 1) * cumRow[i];
                rowMI[i] = (rowMI[i] + carryon) % colRow[i];
                srcPos += carryon * cumNNZ[i]; // rowMI[i] has been incremented
                carryon = carryon * (rowMI[i] == 0 ? 1 : 0);
                srcPos -= carryon * colRow[i] * cumNNZ[i]; // rowMI[i] reached 0
            }
        }
        // set up data for next column
        ++col;
        srcPos = 0;
        for (int i = 0; i < dim; ++i)
        {
            colMI[i] = (col % cumCol[i + 1]) / cumCol[i];
            colRow[i] = (colStart[i][colMI[i] + 1] - colStart[i][colMI[i]]);
            rowMI[i] = 0;
            srcPos += colStart[i][colMI[i]] * cumNNZ[i];
        }
        colEnd = std::min(colStartTot[col + 1], end);
    }
}

template <typename T>
void copyDataIntoSparse(
    std::span<T> nnzDest,
    std::span<mwIndex> rowDest,
    const std::span<T> nnzSrc,
    const std::vector<std::span<mxDouble>> &rowSrc,
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumRow,
    const std::span<mwIndex> cumCol,
    const std::vector<std::vector<mwIndex>> &colStart,
    const std::span<mwIndex> colStartTot,
    size_t numThreads)
{
    int dim = rowSrc.size();
    switch (dim)
    {
    case 0:
        break;
    case 1:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 1>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    case 2:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 2>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    case 3:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 3>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    case 4:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 4>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    case 5:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 5>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    case 6:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 6>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    case 7:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimSpecific<T, 7>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
        break;
    default:
        splitWork(
            cumNNZ.back(), numThreads, minThreadWorkSize,
            copyDataIntoSparseDimGeneric<T>,
            nnzDest, rowDest, nnzSrc, rowSrc,
            cumNNZ, cumRow, cumCol, colStart, colStartTot
        );
    }
}

// Compute col starting points
void computeColStartST(
    const std::vector<std::vector<mwIndex>> &colStart, 
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumCol,
    std::span<mwIndex> dest,
    size_t start, size_t end) 
{
switch (colStart.size())
{
    case 1:
        for (; start<end; ++start)
            dest[start]=colStart[0][start];
        break;
    case 2:
        for (; start<end; ++start)
        {
            size_t c0= ( start /cumCol[0] ) % cumCol[1];
            size_t c1= ( start /cumCol[1] ) % cumCol[2];
            dest[start]=
                colStart[0][c0]*(colStart[1][c1+1]-colStart[1][c1])
                +colStart[1][c1]*cumNNZ[1];
        }
        break;
    case 3:
        for (; start<end; ++start)
        {
            size_t c0= (start % cumCol[1]) / cumCol[0];
            size_t c1= (start % cumCol[2]) / cumCol[1];
            size_t c2= (start % cumCol[3]) / cumCol[2];
            dest[start]=
                (
                colStart[0][c0]*(colStart[1][c1+1]-colStart[1][c1])
                +colStart[1][c1]*cumNNZ[1]
                    )*(colStart[2][c2+1]-colStart[2][c2])
                +colStart[2][c2]*cumNNZ[2]
                ;
        }
        break;
    case 4:
        for (; start<end; ++start)
        {
            size_t c0= (start % cumCol[1] ) / cumCol[0];
            size_t c1= (start % cumCol[2] ) / cumCol[1];
            size_t c2= (start % cumCol[3] ) / cumCol[2];
            size_t c3= (start % cumCol[4] ) / cumCol[3];

            dest[start]=
                ((
                colStart[0][c0]*(colStart[1][c1+1]-colStart[1][c1])
                +colStart[1][c1]*cumNNZ[1]
                    )*(colStart[2][c2+1]-colStart[2][c2])
                +colStart[2][c2]*cumNNZ[2]
                    )*(colStart[3][c3+1]-colStart[3][c3])
                +colStart[3][c3]*cumNNZ[3]
                ;
        }
        break;
    case 5:
        for (; start<end; ++start)
        {
            size_t c0= (start % cumCol[1] ) / cumCol[0];
            size_t c1= (start % cumCol[2] ) / cumCol[1];
            size_t c2= (start % cumCol[3] ) / cumCol[2];
            size_t c3= (start % cumCol[4] ) / cumCol[3];
            size_t c4= (start % cumCol[5] ) / cumCol[4];

            dest[start]=
                (((
                colStart[0][c0]*(colStart[1][c1+1]-colStart[1][c1])
                +colStart[1][c1]*cumNNZ[1]
                +colStart[1][c1]
                    )*(colStart[2][c2+1]-colStart[2][c2])
                +colStart[2][c2]*cumNNZ[2]
                    )*(colStart[3][c3+1]-colStart[3][c3])
                +colStart[3][c3]*cumNNZ[3]
                    )*(colStart[4][c4+1]-colStart[4][c4])
                +colStart[4][c4]*cumNNZ[4]
                ;
        }
        break;
    case 6:
        for (; start<end; ++start)
        {
            size_t c0= (start % cumCol[1] ) / cumCol[0];
            size_t c1= (start % cumCol[2] ) / cumCol[1];
            size_t c2= (start % cumCol[3] ) / cumCol[2];
            size_t c3= (start % cumCol[4] ) / cumCol[3];
            size_t c4= (start % cumCol[5] ) / cumCol[4];
            size_t c5= (start % cumCol[6] ) / cumCol[5];

            dest[start]=
                ((((
                colStart[0][c0]*(colStart[1][c1+1]-colStart[1][c1])
                +colStart[1][c1]*cumNNZ[1]
                +colStart[1][c1]
                    )*(colStart[2][c2+1]-colStart[2][c2])
                +colStart[2][c2]*cumNNZ[2]
                    )*(colStart[3][c3+1]-colStart[3][c3])
                +colStart[3][c3]*cumNNZ[3]
                    )*(colStart[4][c4+1]-colStart[4][c4])
                +colStart[4][c4]*cumNNZ[4]
                    )*(colStart[5][c5+1]-colStart[5][c5])
                +colStart[5][c5]*cumNNZ[5]
                ;
        }    
        break;
    case 7:
        for (; start<end; ++start)
        {           
            size_t c0= (start % cumCol[1] ) / cumCol[0];
            size_t c1= (start % cumCol[2] ) / cumCol[1];
            size_t c2= (start % cumCol[3] ) / cumCol[2];
            size_t c3= (start % cumCol[4] ) / cumCol[3];
            size_t c4= (start % cumCol[5] ) / cumCol[4];
            size_t c5= (start % cumCol[6] ) / cumCol[5];
            size_t c6= (start % cumCol[7] ) / cumCol[6];
            
            dest[start]=
                (((((
                colStart[0][c0]*(colStart[1][c1+1]-colStart[1][c1])
                +colStart[1][c1]*cumNNZ[1]
                +colStart[1][c1]
                    )*(colStart[2][c2+1]-colStart[2][c2])
                +colStart[2][c2]*cumNNZ[2]
                    )*(colStart[3][c3+1]-colStart[3][c3])
                +colStart[3][c3]*cumNNZ[3]
                    )*(colStart[4][c4+1]-colStart[4][c4])
                +colStart[4][c4]*cumNNZ[4]
                    )*(colStart[5][c5+1]-colStart[5][c5])
                +colStart[5][c5]*cumNNZ[5]
                    )*(colStart[6][c6+1]-colStart[6][c6])
                +colStart[6][c6]*cumNNZ[6];
        }
        break;
    default:
        size_t dim=cumCol.size()-1;
        std::vector<mwIndex> c(dim);
        for (; start<end; ++start)
        {
            for (size_t i=0;i<dim;++i)
                c[i] = (start /cumCol[i]) % cumCol[i+1];
            dest[start]=colStart[0][c[0]];   
            for (size_t i=1;i<dim;++i)
            {
                dest[start] *= (colStart[i][c[i]+1]-colStart[i][c[i]]);
                dest[start] += colStart[i][c[i]]*cumNNZ[i];
            }
        }
    }
}


void computeColStart(
    const std::vector<std::vector<mwIndex>> &colStart, 
    const std::span<mwIndex> cumNNZ,
    const std::span<mwIndex> cumCol,
    std::span<mwIndex> dest,
    size_t numThreads)
{
    splitWork(cumCol.back(), numThreads, minThreadWorkSize,
         computeColStartST, colStart,cumNNZ,cumCol, dest );
    dest.back()=cumNNZ.back();
}

template <typename ...T>
void splitWork(size_t totalWorkSize, size_t numThreads, size_t minWorkSize, T...args)
{
    size_t numBl=numThreads;
    if (minWorkSize!=0) numBl=std::max(std::min(numThreads, totalWorkSize/minWorkSize),size_t(1));
    size_t chunk=totalWorkSize/numBl;
    size_t excess=totalWorkSize%numBl;
    std::vector<std::thread> worker;
    size_t b=0; size_t e; size_t j;
    for (j=0;j<numBl;++j) {
        e=b+chunk+(j<excess?1:0);
        worker.emplace_back(args...,b,e);
        b=e;
    }
    for (auto& thread : worker) {
        thread.join();
    }
}

size_t getNumThreads(int nrhs, const mxArray *prhs[])
{
    if (nrhs==4) 
    {
        long int tmp=0;
        if( !mxIsScalar(prhs[3]) )
            mexErrMsgIdAndTxt("tSparsityToSparse:numThreads", "The third argument must be scalar and specifies the number of threads.");   
        
        if      (mxIsDouble(prhs[3]) ) tmp=(long int)(* mxGetDoubles(prhs[3]) );     
        else if (mxIsInt8  (prhs[3]) ) tmp=(long int)(* mxGetInt8s  (prhs[3]) );   
        else if (mxIsInt16 (prhs[3]) ) tmp=(long int)(* mxGetInt16s (prhs[3]) );      
        else if (mxIsInt32 (prhs[3]) ) tmp=(long int)(* mxGetInt32s (prhs[3]) );      
        else if (mxIsInt64 (prhs[3]) ) tmp=(long int)(* mxGetInt64s (prhs[3]) );      
        else if (mxIsUint8 (prhs[3]) ) tmp=(long int)(* mxGetUint8s (prhs[3]) );     
        else if (mxIsUint16(prhs[3]) ) tmp=(long int)(* mxGetUint16s(prhs[3]) );     
        else if (mxIsUint32(prhs[3]) ) tmp=(long int)(* mxGetUint32s(prhs[3]) );     
        else if (mxIsUint64(prhs[3]) ) tmp=(long int)(* mxGetUint64s(prhs[3]) );
        else mexErrMsgIdAndTxt("tSparsityToSparse:numThreads", "num thread can only be of type double, IntX or UInt.");

        if (tmp>0 ) return size_t(tmp); else return std::thread::hardware_concurrency();
    }
    else return std::thread::hardware_concurrency();
}

void addPartialTime(std::vector<logEntry> &time_log, std::string name)
{
    time_log.emplace_back(std::pair(std::chrono::high_resolution_clock::now(), name));
}

mxArray* constructTimeReport(const std::vector<logEntry>& time_log)
{ 
        const char* field_names[]={"time", "label"};
        auto num_fields= (sizeof(field_names) / sizeof(*field_names));
        mwSize dims[2] = {1,  time_log.size()-1};
        mxArray* lhs = mxCreateStructArray(2, dims, num_fields, field_names);

        int time_field =  mxGetFieldNumber(lhs,field_names[0]);
        int label_field = mxGetFieldNumber(lhs,field_names[1]);
        for (size_t i = 1; i < time_log.size(); ++i)
        {
            mxArray* field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
            auto data=getSpan<mxDouble>(field_value);
            data[0]=std::chrono::duration<double>(time_log[i].first - time_log[i - 1].first).count();
            mxSetFieldByNumber(lhs, i, time_field,  field_value);
            mxSetFieldByNumber(lhs, i, label_field, mxCreateString(time_log[i].second.c_str()));
        }
    return lhs;
}

template <typename T>
std::span<T> getSpan(mxArray *arr)
{
    T *beg = (T *) mxGetData(arr);
    mwIndex size = mxGetNumberOfElements(arr);
    if (beg == NULL) mexErrMsgIdAndTxt("TSparsityToSparse:cannotAccessData", "got a null pointer, so I aborted");
    return std::span<T>(beg, beg + size);
}

int checkInputs(int nrhs, const mxArray *prhs[])
{
    if ( !(nrhs == 3 || nrhs == 4) ) 
        mexErrMsgIdAndTxt("tSparsityToSparse:nrhs", "Wrong number of input arguments. Check the documentation.");

    if ( !mxIsCell(prhs[0]) )
        mexErrMsgIdAndTxt("TSparsityToSparse:sparsityData", "The sparsity data must be a cell array.");

    mwSize        nDimSparsityData = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *dimsSparsityData = mxGetDimensions(prhs[0]);

    if ( nDimSparsityData!=2 )
        mexErrMsgIdAndTxt("TSparsityToSparse:sparsityData", "The sparsity data must be a two dimensional cell array with 3 rows.");

    for (mwSize i=0; i<dimsSparsityData[0]*dimsSparsityData[1];++i)
    {
        mxArray *cellElement = mxGetCell(prhs[0], i);
        if ( !mxIsDouble(cellElement) || mxIsSparse(cellElement) || mxIsComplex(cellElement) )
            mexErrMsgIdAndTxt("TSparsityToSparse:sparsityData", "Each cell in the sparsity data should be a full double array.");
    }

    if ( !(mxIsDouble(prhs[1]) || mxIsLogical(prhs[1])) || mxIsSparse(prhs[1]))
        mexErrMsgIdAndTxt("TSparsityToSparse:nonZeroType", "Non Zeros must be a full double or logical array.");

    return dimsSparsityData[1];
}