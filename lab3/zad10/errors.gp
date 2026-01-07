set title "Maksymalne bledy bezwzgledne vs krok czasowy"
set xlabel "log10(dt)"
set ylabel "log10(|MaxError|)"
set grid
set key left top

plot "errors.txt" using 1:2 with linespoints pt 7 title "Euler bezposredni", \
     "errors.txt" using 1:3 with linespoints pt 5 title "Euler posredni", \
     "errors.txt" using 1:4 with linespoints pt 9 title "Metoda trapezow"

pause -1
