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

## Note
In `error_plot.png`, the graph for 2 and 4 cores overlaps.