# Test Folder

To run the test, write on the command line
```bash
sh ./run_test.sh
```
Ensure to have a folder called `data` in this location.

Thanks to the `run_test.sh` script, it is possible to:

+ Creating the executable and the object file related to `main.cpp` file.
+ Run the program with 1,2 and 4 processors, so run the problem sequentially and parallel, and with $n=2^k$ size of the grid, with k=4,5,6,7,8.
+ For each of the combinations of number of process and size of the grid n, collect L2 norm errors and time of computation in the `data` folder.
+ Plot the values colllected by generated two png files ( `error_plot.png` and `time_plot.png`)


On your command line, type
```bash
make clean
```
to remove .png and .dat files

## Content of the folder
+ A data folder, to save the files that contains errors in L2 norm and time of computing, as the size of the grid and number of cores increases (make sure to create it).
+ [plot_errors.gnuplot](https://github.com/irene-fagnani/challenge3/blob/main/test/plot_errors.gnuplot) and [plot_times.gnuplot](https://github.com/irene-fagnani/challenge3/blob/main/test/plot_times.gnuplot), two scripts useful to create the two graphs as .png files.
+ [run_test.sh](https://github.com/irene-fagnani/challenge3/blob/main/test/run_test.sh) script, to run the test with 1,2 and 4 cores and with a  grid of lenght $2^k$, where k=4,5,6,7,8.
+ A [Makefile](https://github.com/irene-fagnani/challenge3/blob/main/test/Makefile) to delete the .dat files in the data folder and the .png graphs.

## Note
In `error_plot.png`, the graph for 1,2 and 4 cores overlaps.