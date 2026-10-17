reset
set terminal pngcairo enhanced size 1400,900 font "Arial,14"
set output "S_all.png"

set xlabel "Итерация k"
set ylabel "S(theta_k)"
set grid
set logscale y
set key top right font "Arial,12"
set title "Сходимость при разных theta_0"

plot \
    "S_theta0_0.000000.txt" using 1:2 with linespoints pt 7 lw 2 lc rgb "blue"   title "theta_0 = 0", \
    "S_theta0_0.300000.txt" using 1:2 with linespoints pt 7 lw 2 lc rgb "red"    title "theta_0 = 0.3", \
    "S_theta0_0.700000.txt" using 1:2 with linespoints pt 7 lw 2 lc rgb "green"  title "theta_0 = 0.7", \
    "S_theta0_1.500000.txt" using 1:2 with linespoints pt 7 lw 2 lc rgb "orange" title "theta_0 = 1.5", \
    "S_theta0_3.000000.txt" using 1:2 with linespoints pt 7 lw 2 lc rgb "purple" title "theta_0 = 3.0", \
    "S_theta0_5.000000.txt" using 1:2 with linespoints pt 7 lw 2 lc rgb "brown"  title "theta_0 = 5.0"

set output
