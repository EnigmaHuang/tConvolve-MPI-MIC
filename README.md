# ASKAP tConvolve(Gridding) Benchmark on Xeon Phi

This is an optimized program for ASKPA tConvolve(Gridding) Benchmark. The original algorithm and source code is here : [ATNF/askap-benchmarks](https://github.com/ATNF/askap-benchmarks/tree/master/tConvolveMPI). The "baseline" code is from ASC15 final contest, and it is slightly different from [ATNF/askap-benchmarks](https://github.com/ATNF/askap-benchmarks/tree/master/tConvolveMPI) version. This "baseline" code is in folder "baseline".



tConvole is a convolutional resampling algorithm, which is used in radio astronomy data processing. The program is configured to reflect the computing needs of the Australian Square Kilometer Array Pathfinder (ASKAP) central processor. 



This optimized program uses Message Passing Interface (MPI) to work on multi-nodes, and uses OpenMP to work on CPU & Xeon Phi coprocessor. This is a **persional re-implementation** and its algorithm is different from the [e-Prize version](http://www.sysu.edu.cn/2012/en/news/new05/22235.htm) in ASC15 (but have the same performance).



To generate different random data, compile *random_data_generator.cpp* and run the program. The program will ask you to enter 4 parameters, then the program will generate *input.dat* and *randnum.dat* . tConvolve program will read these two files. 



To compile tConvolve-MPI-MIC program, you need an ICC 2014 or later version and a Xeon Phi coprocessor (I used the earlist KNL cards). If you don't have a MIC card, *tConvolve-MPI-CPU* folder provides a pure CPU version. Just type `make`  to compile this program. To run the program, the command should looks like `mpirun -np <num_of_nodes> -machinefile <node_list_file> ./tConvolve` , because in a node the program will use OpenMP for parallel(so you just need 1 process in each node). Make sure *input.dat*, *randnum.dat* and *tConvolve* are in the same folder. 



*conf.ini* allows you to override two parameters required by this optimized program: max samples for each computation loop, in the first line; MIC workload balance rate, in the second line. I don't suggest you to change the second line. But you should change the first line according to the memory your computing nodes have. If you have 64GB memory, the first line can be 120000000, bigger is better. 



To verify the result, you can compile *verify.cpp* and use this program to check: put the reference result(in file name *grid_std.dat*) and your result(in file name *grid.dat*) output in the same folder, and run this program.  The *"Relative error of the L1 norm"* should be smaller than 1e-13.



