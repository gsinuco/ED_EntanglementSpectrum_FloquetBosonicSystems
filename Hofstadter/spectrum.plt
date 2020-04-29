reset
set view 90,90
unset xtics
#set xrange [250:500]
set zrange [-1:1]
set ztics ("-1" -1, "0" 0, "1" 1)
set ticslevel 0
unset key
set border
set ylabel "k_y" offset 0,-2
set label "{/Symbol e} ({/Symbol p})" at 3,-25
set yrange [0:128]
set ytics offset 0,-1
set ytics ("-1" 0 , "0" 64, "1" 128)
set terminal postscrip eps color enhanced "Helvetica,22"
set output "spectrum.eps"
splot "datos.dat" matrix w p ps 0.5 pt 5, "datosPeriodic.dat" matrix w p ps 0.5 pt 5
unset output