# Set output image file name and format
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'perf_iter.png'

# Set chart title and axis labels
set title 'Latences de résolution de méthodes Richardson'
set xlabel 'taille de matrice'
set ylabel 'latence (ns)'

# Plot both datasets
plot 'ALPHA' with lines title 'alpha', \
     'JAC' with lines title 'Jacobi', \
     'GS' with lines title 'Gauss-Seidel'