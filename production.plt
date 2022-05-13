set key off
set terminal pdfcairo size 11in,8in 
set xtics (1E-9,1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E0,1E1,1E2,1E3)
#set ytics (5E5,1E6,1.5E6,2E6,2.5E6,3E6)
set format y "%.2tx10^{%L}"
set format x "10^{%L}"
set logscale x
set xlabel "Production Rate [m^{-3} s^{-1} eV^{-1}]" offset 0,-2  font ",28"
set ylabel "Altitude [m]" offset -8,0 font ",28"
set title "Production Rate vs Altitude at 20 eV Photon Energy" font ",28"
set tmargin 10
set rmargin 10
set lmargin 30
set bmargin 20
set tics scale 1
set xtics font ", 22"
set ytics font ", 22"
set yrange [*:*]
set xrange [1E-7:*]
set output 'production.pdf'
plot 'production.txt' lt 7 lw 4 lc 1 w lines

