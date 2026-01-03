set xlabel 'log_{10}(h)'
set ylabel 'log_{10}(|MaxError|)'
set grid
set style data line

plot 'wyniki_h_errors.txt' using 1:2 title '3 points', \
'' using 1:3 title 'shoots', \