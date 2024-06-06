set terminal pngcairo enhanced size 1000,800
set output 'time_plot.png'
set format y "%3.1e"
set title "Computational time as function of matrix dimension"
set xlabel "Grid size n"
set ylabel "Computational time norm error"
set logscale y 10
set grid

set style line 1 lc rgb "red" lw 2 pt 7
set style line 2 lc rgb "blue" lw 2 pt 7
set style line 3 lc rgb "green" lw 2 pt 7
plot "data/times.dat" using ($1 == 1 ? $2 : 1/0):3 with linespoints linestyle 1 title "Core number 1", \
     '' using ($1 == 2 ? $2 : 1/0):3 with linespoints linestyle 2 title "Core number 2", \
     '' using ($1 == 4 ? $2 : 1/0):3 with linespoints linestyle 3 title "Core number 4"