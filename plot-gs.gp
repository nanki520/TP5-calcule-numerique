set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'Richardson_gs.png'

set title 'Historique de convergence'
set xlabel 'Nombre Itérations'
set ylabel 'Norme du résidu relatif'

set xrange [0:130]

plot 'RESVEC.dat' with line title 'Gauss-Seidel'



