reset
set terminal pngcairo enhanced size 1400,900 font "Arial,14"
set output "scene.png"

set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set grid
set title "Поворот эллипсоида"
set key top left font "Arial,12"

set view 60, 30, 1, 1
set ticslevel 0

splot \
    "x.txt"     using 1:2:3 with points pt 7 ps 0.5 lc rgb "blue"  title "Исходные x_i", \
    "y.txt"     using 1:2:3 with points pt 7 ps 0.5 lc rgb "red"   title "Измерения y_i", \
    "y_fit.txt" using 1:2:3 with points pt 7 ps 0.5 lc rgb "green" title "Модель R(theta*) x_i"

set output
