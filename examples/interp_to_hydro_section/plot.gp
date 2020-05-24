#!/usr/bin/gnuplot

set xrange [0:6000]
set yrange [0:16]
plot 'data_in.dat'  t 'model grid' w p pt 4 ps 3, \
     'data_out.dat' t 'target grid' w p pt 2 ps 3
pause 100
