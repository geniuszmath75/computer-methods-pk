set title "Euler bezposredni - graniczny T=5 N=155"
set xlabel "t"
set ylabel "y(t)"
set grid

analytic(x) = 1 + (1+x)**90 * exp(-100*x)

plot analytic(x) with lines title "analityczne", \
     "solution.txt" using 1:3 with points pt 7 title "Euler bezposredni"
pause -1