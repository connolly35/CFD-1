g++ code.cpp
./a.out > toplot
gnuplot -e "plot 'toplot' using 1:2 with lines; pause 10"
#gnuplot -e "set title 'Representations';set xlabel 'x';set ylabel 'y';set zeroaxis;set mxtics 5;set mytics 5;set xtics 0,0.5,7;set ytics -1.0,0.25,1.0;set xrange[0:7];set yrange [-1.5:1.5];plot 'toplot' using 1:2 with lines title 'sin(x)'; pause 5"

#set border N
#N = 1+2+4+8
#1 bottom
#2 left
#4 top
#8 right
# don't forget to use nomirror on xtics and ytics when removing soem borders

#set terminal x11 
#set terminal png picsize 512 512
#set output "cubic.png"
#set size ratio -1.0
