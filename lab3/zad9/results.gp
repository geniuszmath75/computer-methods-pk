set xlabel 'x'
set ylabel 'U(x)'
set grid
set style data line

plot 'wyniki_x_ux.txt' using 1:2 title 'analytic', \
'' using 1:3 title '3 points', \
'' using 1:4 title 'shoots'