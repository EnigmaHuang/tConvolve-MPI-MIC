/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///
/// This program was modified so as to use it in the contest.
/// The last modification was on April 2, 2015.
///

#ifndef BENCHMARK_H
#define BENCHMARK_H

#define BLOCK_LOW(id,p,n) ( (id)*(n) / (p) )
#define BLOCK_HIGH(id,p,n) ( BLOCK_LOW((id)+1,p,n) - 1 )
#define BLOCK_SIZE(id,p,n) ( BLOCK_HIGH(id,p,n) - BLOCK_LOW(id,p,n) + 1 )

// System includes
#include <mpi.h>
#include <complex>
#include <stdio.h>

// Typedefs
typedef double Coord;
typedef float Real;

class Benchmark {
    public:
        Benchmark();

        int randomInt();
        void init();
        void runGrid();

        void gridKernel(const int curr_samplesSize, const int gSize, int *iu,
                        int *iv, int *cOffset, double *data_r, double *data_i);

        void initC(const Coord* freq, const Coord cellSize, const int wSize,
                   int& support, int& overSample, Coord& wCellSize);

        void initCOffset(const Coord* u, const Coord* v, Coord* w, const Coord* freq,
                         const Coord cellSize, const Coord wCellSize, const int wSize,
                         const int gSize, const int support, const int overSample);

        int getSupport();

        int getsSize();

        void printGrid();

// Change these if necessary to adjust run time
        int nSamples_a; // Number of data samples
        int nSamples; // Number of data samples
        int wSize; // Number of lookup planes in w projection
        int nChan; // Number of spectral channels

// Don't change any of these numbers unless you know what you are doing!
        int gSize; // Size of output grid in pixels
        Coord cellSize; // Cellsize of output grid in wavelengths
        int baseline; // Maximum baseline in meters

        int myid;
        int np;
        MPI_Comm comm;


    private:
        Coord* u;
        Coord* v;
        Coord* w;

        // Change 'vector<Sample> samples' to original types
        int samplesSize;
        Coord* samples_data_r;      // Real part of sample data
        Coord* samples_data_i;      // Imaginary part of sample data
        int* samples_iu;
        int* samples_iv;
        int* samples_cOffset;

        // Change "vector<Value> outdata" to original types
        Coord* outdata_r;
        Coord* outdata_i;

        // Change"vector<Value> C" to original types
        Coord* C_r;
        Coord* C_i;

        // Change "Value" to original types
        Coord* grid_r;
        Coord* grid_i;
        Coord* grid0_r;
        Coord* grid0_i;

        int m_support;
        int overSample;
        int sSize;

        Coord wCellSize;

        // For random number generator
        unsigned long next;

        Coord *rd1;
        Coord *rd2;
        Coord *rd3;
        Coord *rd4;

};
#endif
