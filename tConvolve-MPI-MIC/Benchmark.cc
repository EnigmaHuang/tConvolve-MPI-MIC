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
/// The last modification was on November 16, 2015, by Enigma Huang
///

// Include own header file first
#include "Benchmark.h"
#include "Stopwatch.h"

// System includes
#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <omp.h>
#define min(a, b) ((a)<(b)?(a):(b))

Benchmark::Benchmark()
        : next(1)
{
}

int    maxSamples = 100000000;
double MIC_balance_rate = 0.45;

// The following arrays are used for resort and compute
#pragma offload_attribute(push, target(mic))
int    *iu2, *iv2, *cOffset2;
double *data_r2, *data_i2;
int    *above_v;
int    curr_sg_size = 0;
double *sg_r, *sg_i;
double *C_r2, *C_i2;
#pragma offload_attribute(pop)

int has_allocated_on_MIC = 0;

double *hsg_r, *hsg_i;
int *pos_in_iv;

// Return a pseudo-random integer in the range 0..2147483647
// Based on an algorithm in Kernighan & Ritchie, "The C Programming Language"
int Benchmark::randomInt()
{
    const unsigned int maxint = std::numeric_limits<int>::max();
    next = next * 1103515245 + 12345;
    return ((unsigned int)(next / 65536) % maxint);
}

void Benchmark::init()
{
    int rd1_tag,rd2_tag,rd3_tag,rd4_tag,id,nl;
    MPI_Status status;

    rd1_tag=0;
    rd2_tag=1;
    rd3_tag=2;
    rd4_tag=3;

    nSamples=BLOCK_SIZE(myid,np,nSamples_a);

    if(myid == np-1) nl=nSamples+1;
    else nl=nSamples;

    rd1 = new Coord[nl];
    rd2 = new Coord[nl];
    rd3 = new Coord[nl];
    rd4 = new Coord[nl];
    

    // Initialize the data to be gridded
    // u, v, w
    u = new Coord[nSamples];
    v = new Coord[nSamples];
    w = new Coord[nSamples];
    // Samples
    samplesSize = nSamples * nChan;
    samples_data_r = new Coord[samplesSize];
    samples_data_i = new Coord[samplesSize];
    samples_iu = new int[samplesSize];
    samples_iv = new int[samplesSize];
    samples_cOffset = new int[samplesSize];


    if(myid == np-1){   

        Coord rd;
        FILE * fp;
        if( (fp=fopen("randnum.dat","rb"))==NULL ){
            printf("cannot open file\n");
            return;
        }
        for (id = 0; id < np-1; id++){
            nl=BLOCK_SIZE(id,np,nSamples_a);
            for (int i = 0; i < nl; i++) {
                if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
                rd1[i]=rd;
                if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
                rd2[i]=rd;
                if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
                rd3[i]=rd;
                if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
                rd4[i]=rd;
            }
            MPI_Send(rd1,nl,MPI_DOUBLE_PRECISION,id,rd1_tag,comm);
            MPI_Send(rd2,nl,MPI_DOUBLE_PRECISION,id,rd2_tag,comm);
            MPI_Send(rd3,nl,MPI_DOUBLE_PRECISION,id,rd3_tag,comm);
            MPI_Send(rd4,nl,MPI_DOUBLE_PRECISION,id,rd4_tag,comm);
        }
     
        for (int i = 0; i < nSamples; i++) {
            if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
            rd1[i]=rd;
            if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
            rd2[i]=rd;
            if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
            rd3[i]=rd;
            if(fread(&rd,sizeof(Coord),1,fp)!=1){printf("Rand number read error!\n");}
            rd4[i]=rd;
        }

        fclose(fp);
    }
    else{
        MPI_Recv(rd1,nSamples,MPI_DOUBLE_PRECISION,np-1,rd1_tag,comm,&status);
        MPI_Recv(rd2,nSamples,MPI_DOUBLE_PRECISION,np-1,rd2_tag,comm,&status);
        MPI_Recv(rd3,nSamples,MPI_DOUBLE_PRECISION,np-1,rd3_tag,comm,&status);
        MPI_Recv(rd4,nSamples,MPI_DOUBLE_PRECISION,np-1,rd4_tag,comm,&status);
    }

    for (int i = 0; i < nSamples; i++) {
        u[i] = baseline * rd1[i] - baseline / 2;
        v[i] = baseline * rd2[i] - baseline / 2;
        w[i] = baseline * rd3[i] - baseline / 2;
        for (int chan = 0; chan < nChan; chan++) {
            Coord c2=Coord(chan)/Coord(nChan);
            samples_data_r[i*nChan+chan] = rd4[i] + c2;
            samples_data_i[i*nChan+chan] = rd4[i] - c2;
        }
    }
    
    // grid =new Value[gSize*gSize];
    grid_r =new Coord[gSize*gSize];
    grid_i =new Coord[gSize*gSize];
    memset(grid_r, 0, sizeof(Coord)*gSize*gSize);
    memset(grid_i, 0, sizeof(Coord)*gSize*gSize);
    
    if (myid == 0) 
    {
        grid0_r = new Coord[gSize*gSize];
        grid0_i = new Coord[gSize*gSize];
        memset(grid0_r, 0, sizeof(Coord)*gSize*gSize);
        memset(grid0_i, 0, sizeof(Coord)*gSize*gSize);
    }


    // Measure frequency in inverse wavelengths
    Coord* freq = new Coord[nChan];

    for (int i = 0; i < nChan; i++) {
        freq[i] = (1.4e9 - 2.0e5 * Coord(i) / Coord(nChan)) / 2.998e8;
    }

    // Initialize convolution function and offsets
    initC(freq, cellSize, wSize, m_support, overSample, wCellSize);
    initCOffset(u, v, w, freq, cellSize, wCellSize, wSize, gSize,
                m_support, overSample);

    delete [] freq;
    delete [] rd1;
    delete [] rd2;
    delete [] rd3;
    delete [] rd4;
}

void allocComputeArrays(const int samplesSize, const int gSize, const int sgsize, 
                        const int C_num, double *C_r, double *C_i)
{
    Stopwatch sw;
    
    sw.start();
   
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            iu2       = (int*) malloc(sizeof(int) * samplesSize);
            cOffset2  = (int*) malloc(sizeof(int) * samplesSize);
            data_r2   = (double*) malloc(sizeof(double) * samplesSize);
            data_i2   = (double*) malloc(sizeof(double) * samplesSize);
            above_v   = (int*) malloc(sizeof(int) * gSize);
            hsg_r     = (double*) malloc(sizeof(double) * sgsize);
            hsg_i     = (double*) malloc(sizeof(double) * sgsize);
            pos_in_iv = (int*) malloc(sizeof(int) * samplesSize);
        }

        #pragma omp section
        {
            C_r2 = (double*) malloc(sizeof(double) * C_num);
            C_i2 = (double*) malloc(sizeof(double) * C_num);
            memcpy(C_r2, C_r, sizeof(double) * C_num);
            memcpy(C_i2, C_i, sizeof(double) * C_num);
        }
    }
   
    double ut = sw.stop();
    long caSize = samplesSize * (sizeof(int) * 3 + sizeof(double) * 2);
    caSize += sizeof(int) * gSize + C_num * 2 * sizeof(double);
    printf("First CA memory allocation     : %d (MB) \n", caSize / 1024 / 1024);
    printf("Used time                      = %lf (s) \n", ut);
}

void allocComputeArraysOnMIC(const int gSize, const int samplesSize, const int C_num)
{
    #pragma offload target(mic) \
    nocopy(above_v : length(gSize) alloc_if(1) free_if(0)) \
    nocopy(iu2 : length(samplesSize) alloc_if(1) free_if(0)) \
    nocopy(cOffset2 : length(samplesSize) alloc_if(1) free_if(0)) \
    nocopy(data_r2 : length(samplesSize) alloc_if(1) free_if(0)) \
    nocopy(data_i2 : length(samplesSize) alloc_if(1) free_if(0)) \
    in(C_r2 : length(C_num) alloc_if(1) free_if(0)) \
    in(C_i2 : length(C_num) alloc_if(1) free_if(0))
    {}
}

void check_sg_size(const int new_sg_size)
{
    if (new_sg_size > curr_sg_size)
    {
        Stopwatch sw;
        sw.start();
        if (curr_sg_size != 0)
        {
            free(sg_r);
            free(sg_i);
            free(hsg_r);
            free(hsg_i);
            #pragma offload target(mic) \
            nocopy(sg_r : length(curr_sg_size) alloc_if(0) free_if(1)) \
            nocopy(sg_i : length(curr_sg_size) alloc_if(0) free_if(1)) 
            {
            }
        }
        int msize = sizeof(double) * new_sg_size;
        sg_r = (double*) malloc(msize);
        sg_i = (double*) malloc(msize);
        hsg_r = (double*) malloc(msize);
        hsg_i = (double*) malloc(msize);
        curr_sg_size = new_sg_size;
        #pragma offload target(mic) \
        nocopy(sg_r : length(curr_sg_size) alloc_if(1) free_if(0)) \
        nocopy(sg_i : length(curr_sg_size) alloc_if(1) free_if(0)) 
        {
            memset(sg_r, 0, sizeof(double) * curr_sg_size);
            memset(sg_i, 0, sizeof(double) * curr_sg_size);
        }
        printf("Realloc for small grid, size   = %d (MB) \n", msize * 4 / 1024 / 1024);
        printf("      |---  used time = %lf (s)\n", sw.stop());
    }
    memset(hsg_r, 0, sizeof(double) * new_sg_size);
    memset(hsg_i, 0, sizeof(double) * new_sg_size);
}

void freeComputeArrays(const int gSize, const int C_num)
{
    free(iu2);
    free(cOffset2);
    free(data_r2);
    free(data_i2);
    free(sg_r);
    free(sg_i);
    free(above_v);
    free(hsg_r);
    free(hsg_i);
    
    #pragma offload target(mic) \
    nocopy(above_v : length(gSize) alloc_if(0) free_if(1)) \
    nocopy(iu2 : length(maxSamples) alloc_if(0) free_if(1)) \
    nocopy(cOffset2 : length(maxSamples) alloc_if(0) free_if(1)) \
    nocopy(data_r2 : length(maxSamples) alloc_if(0) free_if(1)) \
    nocopy(data_i2 : length(maxSamples) alloc_if(0) free_if(1)) \
    nocopy(sg_r : length(curr_sg_size) alloc_if(0) free_if(1)) \
    nocopy(sg_i : length(curr_sg_size) alloc_if(0) free_if(1)) \
    nocopy(C_r2 : length(C_num) alloc_if(0) free_if(1)) \
    nocopy(C_i2 : length(C_num) alloc_if(0) free_if(1)) 
    {}
}

void Benchmark::runGrid()
{
    
    int finishedSamples = 0;
    int loop = 0;
    int C_num = sSize * sSize * overSample * overSample * wSize;
    int sgsize = (baseline + sSize) * (baseline + sSize);
    
    FILE * conf = fopen("conf.ini", "r");
    if (conf != NULL)
    {
        fscanf(conf, "%d", &maxSamples);
        fscanf(conf, "%lf", &MIC_balance_rate);
        printf("Parameters overrided from conf.ini:\n");
        printf("Max samples for each loop = %d\n", maxSamples);
        printf("MIC workload offload rate = %lf\n", MIC_balance_rate);
        fclose(conf);
    }

    printf("Total samples size = %d, channel number = %d\n", samplesSize, nChan);
    
    allocComputeArrays(maxSamples, gSize, sgsize, C_num, C_r, C_i);
   
    while (finishedSamples < samplesSize)
    {
        loop++;
        int curr_samplesSize = min(maxSamples, samplesSize - finishedSamples);
        printf("Loop %3d, samples size         = %d\n", loop, curr_samplesSize);
        gridKernel(curr_samplesSize, gSize, samples_iu + finishedSamples, 
                   samples_iv + finishedSamples, samples_cOffset + finishedSamples,
                   samples_data_r + finishedSamples, samples_data_i + finishedSamples);
        finishedSamples += curr_samplesSize;
    }
    
    freeComputeArrays(gSize, C_num);
    
    MPI_Reduce(grid_r,grid0_r,gSize*gSize,MPI_DOUBLE,MPI_SUM,0,comm); 
    MPI_Reduce(grid_i,grid0_i,gSize*gSize,MPI_DOUBLE,MPI_SUM,0,comm); 
}

void locateEdge
(
    int *iv, int *iu, int samplesSize, int *above_v, 
    int gSize, int &minu, int &maxu, int &minv, 
    int &maxv, int &lenu, int &lenv
)
{
    // above_v[v1] = x1 means there are x1 samples whose iv value < v1
    memset(above_v, 0, sizeof(int) * gSize);
    minv = maxv = iv[0];
    minu = maxu = iu[0];
    for (int i = 0; i < samplesSize; i++)
    {
        pos_in_iv[i] = above_v[iv[i] + 1];
        above_v[iv[i] + 1]++;
        minv = minv < iv[i] ? minv : iv[i];
        minu = minu < iu[i] ? minu : iu[i];
        maxv = maxv > iv[i] ? maxv : iv[i];
        maxu = maxu > iu[i] ? maxu : iu[i];
    }
    
    lenu = maxu - minu + 1;
    lenv = maxv - minv + 1;

    // Move the data to the head position
    for (int i = 0; i <= lenv; i++)
        above_v[i] = above_v[i + minv];
    for (int i = lenv + 1; i < gSize; i++)
        above_v[i] = 0;
    
    for (int i = 1; i < gSize; i++)
        above_v[i] += above_v[i - 1];
}

inline void swap_int(int &a, int &b)
{
    int swap = a;
    a = b; b = swap;
}

inline void swap_double(double &a, double &b)
{
    double swap = a;
    a = b; b = swap;
}

inline int cmp(int &mainkey1, int &subkey1, int &mainkey2, int &subkey2)
{
    if (mainkey1 == mainkey2) return subkey1 < subkey2;
    else return mainkey1 < mainkey2;
}

void quickSort(int l, int r, int *iu, int *cOffset, double *data_r, double *data_i)
{
    int i = l, j = r, mid = (i + j) / 2;
    int mid_u = iu[mid];
    int mid_cOffset = cOffset[mid];
    while (i <= j)
    {
        while (cmp(cOffset[i], iu[i], mid_cOffset, mid_u)) i++;
        while (cmp(mid_cOffset, mid_u, cOffset[j], iu[j])) j--;
        if (i <= j)
        {
            swap_int(iu[i], iu[j]);
            swap_int(cOffset[i], cOffset[j]);
            swap_double(data_r[i], data_r[j]);
            swap_double(data_i[i], data_i[j]);
            i++; j--;
        }
    }
    if (i < r) quickSort(i, r, iu, cOffset, data_r, data_i);
    if (l < j) quickSort(l, j, iu, cOffset, data_r, data_i);              
}

void resortSamplesByV
(
     int samplesSize, int *above_v, int minv, int lenv, 
     int minu, int *iu, int *iv, int *cOffset, 
     double *data_r, double *data_i, int *iu2, 
     int *cOffset2, double *data_r2, double *data_i2
)
{
    Stopwatch sw;

    sw.start();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < samplesSize; i++)
    {
        int newpos = above_v[iv[i] - minv] + pos_in_iv[i];
        iu2[newpos] = iu[i] - minu;
        iv2[newpos] = iv[i] - minv;
        cOffset2[newpos] =cOffset[i];
        data_r2[newpos] = data_r[i];
        data_i2[newpos] = data_i[i];
    }
    printf("Para-colle used time = %lf (s)\n", sw.stop());
    
    sw.start();
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i <= lenv; i++)
    {
        int spos;
        if (i == 0) spos = 0;
        else spos = above_v[i];
        int epos = above_v[i + 1] - 1;
        quickSort(spos, epos, iu2, cOffset2, data_r2, data_i2);
    }
    printf("Para-sort  used time = %lf (s)\n", sw.stop());
}

__attribute__((target(mic)))
void MICInnerKernel
(
    int curr_v, int *above_v, double *sg_r, 
    double *sg_i, int lenv, int small_grid_lenu, 
    int sSize, int *samples_iu, int *samples_cOffset,
    double *samples_data_r, double *samples_data_i,
    double *C_r, double *C_i
)
{
    int start_v = curr_v - sSize + 1;
    int end_v = curr_v;
    if (start_v < 0) start_v = 0;
    if (end_v > lenv - 1) end_v = lenv - 1;
    
    int wgsize = sizeof(double) * small_grid_lenu;
    double *write_gr = (double*) malloc(wgsize);
    double *write_gi = (double*) malloc(wgsize);
    memset(write_gi, 0, wgsize);
    memset(write_gr, 0, wgsize);
    
    for (int i = start_v; i <= end_v; i++)
    {
        int vcOffset = (curr_v - i) * sSize;
        
        double prev_dr = samples_data_r[above_v[i]];
        double prev_di = samples_data_i[above_v[i]];
        int prev_iu = samples_iu[above_v[i]];
        int prev_cind = samples_cOffset[above_v[i]] + vcOffset;
        
        // VERY STRANGE... here should be dind < above_v[i + 1], and the comment
        // below should be run normally, but it DO NOT
        for (int dind = above_v[i] + 1; dind <= above_v[i + 1]; dind++)
        {
            int cind = samples_cOffset[dind] + vcOffset;
            int iu = samples_iu[dind];

            // For same cind and iu, accumulate the data point...
            if ((cind == prev_cind) && (iu == prev_iu))
            {
                prev_dr += samples_data_r[dind];
                prev_di += samples_data_i[dind];
                continue; // ... and go to next iu & cind
            }
           
            // This (iu, cind) pair is different from before, 
            // so compute the data that accumulated
            #pragma simd
            for (int suppu = 0; suppu < sSize; suppu++)
            {
                double res1 = prev_dr * C_r[prev_cind + suppu];
                double res2 = prev_di * C_r[prev_cind + suppu];
                res1 -= prev_di * C_i[prev_cind + suppu];
                res2 += prev_dr * C_i[prev_cind + suppu];
                write_gr[suppu + prev_iu] += res1;
                write_gi[suppu + prev_iu] += res2;
            }
            
            // Remember to reset the (iu, cind) record and data point value
            prev_dr = samples_data_r[dind];
            prev_di = samples_data_i[dind];
            prev_iu = iu;
            prev_cind = cind;
        }
       
        /*
        // Compute the last accumulated data?
        #pragma simd
        for (int suppu = 0; suppu < sSize; suppu++)
        {
            write_gr[suppu + prev_iu] += prev_dr * C_r[prev_cind + suppu] - prev_di * C_i[prev_cind + suppu];
            write_gi[suppu + prev_iu] += prev_di * C_r[prev_cind + suppu] + prev_dr * C_i[prev_cind + suppu];
        }
        */
    }
    
    // Each row's width is small_grid_lenu, now is the curr_v-th row
    // So the write back offset should be small_grid_lenu * curr_v
    memcpy(sg_r + small_grid_lenu * curr_v, write_gr, wgsize);
    memcpy(sg_i + small_grid_lenu * curr_v, write_gi, wgsize);
    
    free(write_gi);
    free(write_gr);

    //printf("iv %d done\n", curr_v);
}

void CPUInnerKernel
(
    int curr_v, int *above_v, double *sg_r, 
    double *sg_i, int lenv, int small_grid_lenu, 
    int sSize, int *samples_iu, int *samples_cOffset,
    double *samples_data_r, double *samples_data_i,
    double *C_r, double *C_i
)
{
    int start_v = curr_v - sSize + 1;
    int end_v = curr_v;
    if (start_v < 0) start_v = 0;
    if (end_v > lenv - 1) end_v = lenv - 1;
    
    int wgsize = sizeof(double) * small_grid_lenu;
    double *write_gr = (double*) malloc(wgsize);
    double *write_gi = (double*) malloc(wgsize);
    memset(write_gi, 0, wgsize);
    memset(write_gr, 0, wgsize);
    
    for (int i = start_v; i <= end_v; i++)
    {
        int vcOffset = (curr_v - i) * sSize;
        
        double prev_dr = samples_data_r[above_v[i]];
        double prev_di = samples_data_i[above_v[i]];
        int prev_iu = samples_iu[above_v[i]];
        int prev_cind = samples_cOffset[above_v[i]] + vcOffset;
        
        // VERY STRANGE... here should be dind < above_v[i + 1], and the comment
        // below should be run normally, but it DO NOT
        for (int dind = above_v[i] + 1; dind <= above_v[i + 1]; dind++)
        {
            int cind = samples_cOffset[dind] + vcOffset;
            int iu = samples_iu[dind];
            double dr = samples_data_r[dind];
            double di = samples_data_i[dind];

            // For same cind and iu, accumulate the data point...
            if ((cind == prev_cind) && (iu == prev_iu))
            {
                prev_dr += samples_data_r[dind];
                prev_di += samples_data_i[dind];
                continue; // ... and go to next iu & cind
            }
           
            // This (iu, cind) pair is different from before, 
            // so compute the data that accumulated
            #pragma simd
            for (int suppu = 0; suppu < sSize; suppu++)
            {
                write_gr[suppu + prev_iu] += prev_dr * C_r[prev_cind + suppu] - prev_di * C_i[prev_cind + suppu];
                write_gi[suppu + prev_iu] += prev_di * C_r[prev_cind + suppu] + prev_dr * C_i[prev_cind + suppu];
            }

            // Remember to reset the (iu, cind) record and data point value
            prev_dr = dr;
            prev_di = di;
            prev_iu = iu;
            prev_cind = cind;
        }
       
        /*
        // Compute the last accumulated data?
        #pragma simd
        for (int suppu = 0; suppu < sSize; suppu++)
        {
            write_gr[suppu + prev_iu] += prev_dr * C_r[prev_cind + suppu] - prev_di * C_i[prev_cind + suppu];
            write_gi[suppu + prev_iu] += prev_di * C_r[prev_cind + suppu] + prev_dr * C_i[prev_cind + suppu];
        }
        */
    }
    
    // Each row's width is small_grid_lenu, now is the curr_v-th row
    // So the write back offset should be small_grid_lenu * curr_v
    memcpy(sg_r + small_grid_lenu * curr_v, write_gr, wgsize);
    memcpy(sg_i + small_grid_lenu * curr_v, write_gi, wgsize);
    
    free(write_gi);
    free(write_gr);

    //printf("iv %d done\n", curr_v);
}

__attribute__((target(mic)))
void micKernel
(
    int *above_v, double *sg_r, double *sg_i, int lenv, 
    int small_grid_lenu, int sSize, int mic_lenv,
    int *iu, int *cOffset, double *data_r, 
    double *data_i, double *C_r, double *C_i
)
{
    #pragma omp parallel for schedule(dynamic)
    for (int curr_v = 0; curr_v < mic_lenv; curr_v++)
        MICInnerKernel
        (
            curr_v, above_v, sg_r, sg_i, lenv, 
            small_grid_lenu, sSize, iu,  
            cOffset, data_r, data_i, C_r, C_i
        );
}

void write_sg_back_to_grid(const int small_grid_lenv, const int small_grid_lenu, 
                           const int minv, const int minu, const int gSize, const int mic_lenv,
                           double *grid_r, double *grid_i, double *sg_r, double *sg_i,
                           double *hsg_r, double *hsg_i)
{
    // Copy the MIC result with the CPU result
    #pragma omp parallel for schedule(dynamic)
    for (int v = mic_lenv; v < small_grid_lenv; v++)
        for (int u = 0; u <small_grid_lenu; u++)
        {
            sg_r[v * small_grid_lenv + u] += hsg_r[v * small_grid_lenv + u];
            sg_i[v * small_grid_lenv + u] += hsg_i[v * small_grid_lenv + u];
        }
    
    // Write the result in small grid back to original grid
    #pragma omp parallel for schedule(dynamic)
    for (int v = 0; v < small_grid_lenv; v++)
        for (int u = 0; u <small_grid_lenu; u++)
        {
            grid_r[(v + minv) * gSize + u + minu] += sg_r[v * small_grid_lenv + u];
            grid_i[(v + minv) * gSize + u + minu] += sg_i[v * small_grid_lenv + u];
        }
}

/////////////////////////////////////////////////////////////////////////////////
// The next function is the kernel of the gridding.
// The data are presented as a vector. Offsets for the convolution function
// and for the grid location are precalculated so that the kernel does
// not need to know anything about world coordinates or the shape of
// the convolution function. The ordering of cOffset and iu, iv is
// random.
//
// Perform gridding
//
// data - values to be gridded in a 1D vector
// support - Total width of convolution function=2*support+1
// C - convolution function shape: (2*support+1, 2*support+1, *)
// cOffset - offset into convolution function per data point
// iu, iv - integer locations of grid points
// grid - Output grid: shape (gSize, *)
// gSize - size of one axis of grid
void Benchmark::gridKernel
(
    const int curr_samplesSize, const int gSize, 
    int *iu, int *iv, int *cOffset, double *data_r, double *data_i
)
{
    int C_num = sSize*sSize*overSample*overSample*wSize;
    int minu, maxu, minv, maxv, lenu, lenv; 

    Stopwatch sw, sw1, sw2;
    
    sw.start();
    locateEdge(iv, iu, curr_samplesSize, above_v, 
               gSize, minu, maxu, minv, maxv, lenu, lenv);
    printf("Locate edge used time          = %lf (s) \n", sw.stop());
    
    sw.start();
    resortSamplesByV(curr_samplesSize, above_v, minv, lenv, minu,
                     iu, iv, cOffset, data_r, data_i,
                     iu2, cOffset2, data_r2, data_i2);
    printf("Resort samples by iv used time = %lf (s) \n", sw.stop());
    
    int small_grid_lenu = lenu + sSize - 1;
    int small_grid_lenv = lenv + sSize - 1;
    int small_grid_size = small_grid_lenu * small_grid_lenv;
    
    check_sg_size(small_grid_size);
    
    int mic_lenv = (int)(MIC_balance_rate * small_grid_lenv);
    int MIC_data_size = above_v[mic_lenv + 1];
    if (!has_allocated_on_MIC) mic_lenv = 0;
    
    omp_set_nested(true);
    
    sw.start();

    // Compute kernel
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            sw1.start();
            #pragma omp parallel for schedule(dynamic)
            for (int curr_v = mic_lenv; curr_v < small_grid_lenv; curr_v++)
            {
                CPUInnerKernel
                (
                    curr_v, above_v, hsg_r, hsg_i, lenv, 
                    small_grid_lenu, sSize, iu2, 
                    cOffset2, data_r2, data_i2, C_r, C_i
                );
            }
            printf("CPU kernel used time = %lf (s)\n", sw1.stop());
        }

        #pragma omp section
        {
            if (!has_allocated_on_MIC) 
            {
                printf("First shoot, allocate arrays on MIC but do not compute\n");
                allocComputeArraysOnMIC(gSize, maxSamples, C_num);
                has_allocated_on_MIC = 1;
            }
            else
            {
                sw2.start();
                #pragma offload target(mic) \
                in(above_v : length(gSize) alloc_if(0) free_if(0)) \
                in(iu2 : length(MIC_data_size) alloc_if(0) free_if(0)) \
                in(cOffset2 : length(MIC_data_size) alloc_if(0) free_if(0)) \
                in(data_r2 : length(MIC_data_size) alloc_if(0) free_if(0)) \
                in(data_i2 : length(MIC_data_size) alloc_if(0) free_if(0)) 
                {}
                printf("MIC trans-in data used time = %lf (s)\n", sw2.stop());
                
                sw2.start();
                #pragma offload target(mic) \
                nocopy(above_v : length(gSize) alloc_if(0) free_if(0)) \
                out(sg_r : length(small_grid_size) alloc_if(0) free_if(0)) \
                out(sg_i : length(small_grid_size) alloc_if(0) free_if(0)) \
                in(lenv) in(small_grid_lenu) in(sSize) in(mic_lenv)\
                nocopy(iu2 : length(maxSamples) alloc_if(0) free_if(0)) \
                nocopy(cOffset2 : length(maxSamples) alloc_if(0) free_if(0)) \
                nocopy(data_r2 : length(maxSamples) alloc_if(0) free_if(0)) \
                nocopy(data_i2 : length(maxSamples) alloc_if(0) free_if(0)) \
                nocopy(C_r2 : length(C_num) alloc_if(0) free_if(0)) \
                nocopy(C_i2 : length(C_num) alloc_if(0) free_if(0))
                {
                   micKernel(above_v, sg_r, sg_i, lenv, small_grid_lenu, 
                              sSize, mic_lenv, iu2, cOffset2, 
                              data_r2, data_i2, C_r2, C_i2);
                   
                }
                printf("MIC kernel used time = %lf (s)\n", sw2.stop());
            }
        }
    }
    
    printf("Compute kernel used time       = %lf (s)\n", sw.stop());

    sw.start();
    write_sg_back_to_grid(small_grid_lenv, small_grid_lenu, minv,  
                          minu, gSize, mic_lenv, grid_r, grid_i, 
                          sg_r, sg_i, hsg_r, hsg_i);
    printf("Writeback sg to grid used time = %lf (s)\n", sw.stop());
}

/////////////////////////////////////////////////////////////////////////////////
// Initialize W project convolution function
// - This is application specific and should not need any changes.
//
// freq - temporal frequency (inverse wavelengths)
// cellSize - size of one grid cell in wavelengths
// wSize - Size of lookup table in w
// support - Total width of convolution function=2*support+1
// wCellSize - size of one w grid cell in wavelengths
void Benchmark::initC(const Coord* freq, const Coord cellSize, const int wSize,
                      int& support, int& overSample, Coord& wCellSize)
{
    support = static_cast<int>(1.5 * sqrt(std::abs(baseline) * static_cast<Coord>(cellSize)
                                          * freq[0]) / cellSize);

    overSample = 8;
    wCellSize = 2 * baseline * freq[0] / wSize;

    // Convolution function. This should be the convolution of the
    // w projection kernel (the Fresnel term) with the convolution
    // function used in the standard case. The latter is needed to
    // suppress aliasing. In practice, we calculate entire function
    // by Fourier transformation. Here we take an approximation that
    // is good enough.
    sSize = 2 * support + 1;

    const int cCenter = (sSize - 1) / 2;

    // C.resize(sSize*sSize*overSample*overSample*wSize);
    C_r = new Coord[sSize*sSize*overSample*overSample*wSize];
    C_i = new Coord[sSize*sSize*overSample*overSample*wSize];

    double rr,ri;
    for (int k = 0; k < wSize; k++) {
        double w = double(k - wSize / 2);
        double fScale = sqrt(std::abs(w) * wCellSize * freq[0]) / cellSize;

        for (int osj = 0; osj < overSample; osj++) {
            for (int osi = 0; osi < overSample; osi++) {
                for (int j = 0; j < sSize; j++) {
                    double j2 = std::pow((double(j - cCenter) + double(osj) / double(overSample)), 2);

                    for (int i = 0; i < sSize; i++) {
                        double i2 = std::pow((double(i - cCenter) + double(osi) / double(overSample)), 2);
                        double r2 = j2 + i2 + sqrt(j2*i2);
                        long int cind = i + sSize * (j + sSize * (osi + overSample * (osj + overSample * k)));

                        if (w != 0.0) {
                            rr=std::cos(r2 / (w * fScale));
                            ri=std::sin(r2 / (w * fScale));
                            // C[cind] = static_cast<Value>(rr,ri);
                            // Correct statement must be like these:
                            //C_r[cind] = rr;
                            //C_i[cind] = ri;
                            // But inspur has make a mistake
                            C_r[cind] = ri;
                            C_i[cind] = 0.0;
                        } else {
                            rr=std::exp(-r2);
                            // C[cind] = static_cast<Value>(rr);
                            C_r[cind] = rr;
                            C_i[cind] = 0;
                        }
                    }
                }
            }
        }
    }

    // Now normalise the convolution function
    Coord sumC = 0.0;

    for (int i = 0; i < sSize*sSize*overSample*overSample*wSize; i++) {
      sumC += sqrt(C_r[i] * C_r[i] + C_i[i] * C_i[i]);
    }

    Coord norm = wSize * overSample * overSample / sumC;
    for (int i = 0; i < sSize*sSize*overSample*overSample*wSize; i++) {
        C_r[i] *= norm;
    }
    for (int i = 0; i < sSize*sSize*overSample*overSample*wSize; i++) {
        C_i[i] *= norm;
    }
}

// Initialize Lookup function
// - This is application specific and should not need any changes.
//
// freq - temporal frequency (inverse wavelengths)
// cellSize - size of one grid cell in wavelengths
// gSize - size of grid in pixels (per axis)
// support - Total width of convolution function=2*support+1
// wCellSize - size of one w grid cell in wavelengths
// wSize - Size of lookup table in w
void Benchmark::initCOffset(const Coord* u, const Coord* v, Coord* w, const Coord* freq,
                            const Coord cellSize, const Coord wCellSize,
                            const int wSize, const int gSize, const int support,
                            const int overSample)
{

    // Now calculate the offset for each visibility point
    for (int i = 0; i < nSamples; i++) {
        for (int chan = 0; chan < nChan; chan++) {

            int dind = i * nChan + chan;

            Coord uScaled = freq[chan] * u[i] / cellSize;
            samples_iu[dind] = int(uScaled);

            if (uScaled < Coord(samples_iu[dind])) {
                samples_iu[dind] -= 1;
            }

            int fracu = int(overSample * (uScaled - Coord(samples_iu[dind])));
            samples_iu[dind] += gSize / 2;

            Coord vScaled = freq[chan] * v[i] / cellSize;
            samples_iv[dind] = int(vScaled);

            if (vScaled < Coord(samples_iv[dind])) {
                samples_iv[dind] -= 1;
            }

            int fracv = int(overSample * (vScaled - Coord(samples_iv[dind])));
            samples_iv[dind] += gSize / 2;

            // The beginning of the convolution function for this point
            Coord wScaled = freq[chan] * w[i] / wCellSize;
            int woff = wSize / 2 + int(wScaled);
            samples_cOffset[dind] = sSize * sSize * (fracu + overSample * (fracv + overSample * woff));
        }
    }
}

void Benchmark::printGrid()
{
  FILE * fp;
  if( (fp=fopen("grid.dat","wb"))==NULL )
  {
    printf("cannot open file\n");
    return;
  }  

  unsigned ij;
  for (int i = 0; i < gSize; i++)
  {
    for (int j = 0; j < gSize; j++)
    {
      ij=j+i*gSize;
      if(fwrite(&grid0_r[ij],sizeof(Coord),1,fp)!=1)
        printf("File write error!\n"); 
      if(fwrite(&grid0_i[ij],sizeof(Coord),1,fp)!=1)
        printf("File write error!\n"); 
    }
  }
  
  fclose(fp);
}

int Benchmark::getsSize()
{
    return sSize;
}

int Benchmark::getSupport()
{
    return m_support;
};
