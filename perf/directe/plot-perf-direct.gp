# Set output image file name and format
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'perf_direct.png'

# Set chart title and axis labels
set title 'Latences de résolution de 2 méthodes'
set xlabel 'taille de matrice'
set ylabel 'latence (ns)'

# Plot both datasets
plot 'SV' with lines title 'dgbsv', \
     'TRF' with lines title 'dgbtrf+dgbtrs'