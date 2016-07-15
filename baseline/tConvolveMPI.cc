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

// Include own header file first
#include <mpi.h>
#include "tConvolveMPI.h"
#include <iostream>

// Local includes
#include "Benchmark.h"
#include "Stopwatch.h"
#include <cstdlib>

// Main testing routine
int main(int argc, char *argv[])
{
    int rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Setup the benchmark class
    Benchmark bmark;
    bmark.myid=rank;
    bmark.np=numtasks;
    bmark.comm=MPI_COMM_WORLD;

  int inputi[5];
  Coord cellSize;

  if(rank == 0){
    FILE * fp;
    if( (fp=fopen("input.dat","r"))==NULL )
    {
      printf("cannot open input file\n");
      return 1;
    }
    fscanf(fp,"nSamples=%d\n",&bmark.nSamples_a);
    fscanf(fp,"wSize=%d\n",&bmark.wSize);
    fscanf(fp,"nChan=%d\n",&bmark.nChan);
    fscanf(fp,"gSize=%d\n",&bmark.gSize);
    fscanf(fp,"baseline=%d\n",&bmark.baseline);
    fscanf(fp,"cellSize=%lf\n",&bmark.cellSize);
    fclose(fp);

    inputi[0]=bmark.nSamples_a;
    inputi[1]=bmark.wSize;
    inputi[2]=bmark.nChan;
    inputi[3]=bmark.gSize;
    inputi[4]=bmark.baseline;
    cellSize=bmark.cellSize;
  }
    
    MPI_Bcast(inputi,5,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&cellSize,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);

  if(rank != 0)
  {
    bmark.nSamples_a=inputi[0];
    bmark.wSize=inputi[1];
    bmark.nChan=inputi[2];
    bmark.gSize=inputi[3];
    bmark.baseline=inputi[4];
    bmark.cellSize=cellSize;
  }
    
    bmark.init();

    // Determine how much work will be done 
    const int sSize = bmark.getsSize();
    const double griddings = (double(bmark.nSamples_a * bmark.nChan) * double((sSize) * (sSize))); 

    Stopwatch sw;
    sw.start();
    bmark.runGrid();
    double time_l = sw.stop();
    double time;
    MPI_Reduce(&time_l,&time,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD);
 

  if(rank == 0){
    FILE * flog;
    // Report on timings (master reports only)
    std::cout << "    Time " << time << " (s) " << std::endl;
    std::cout << "    Gridding rate   " << (griddings / 1000000) / time << " (million grid points per second)" << std::endl;

    if( (flog=fopen("log.dat","w"))==NULL )
    {
      printf("cannot open log file\n");
      return 1;
    }
    fprintf(flog,"Time %f (s)\n",time);
    fprintf(flog,"Gridding rate %f (million grid points per second)\n",(griddings/1000000)/time);
    fclose(flog);

    // Output the grid array.
    bmark.printGrid();
  }

    MPI_Finalize();

    return 0;
}
