set key off
set terminal pdfcairo size 8in,11in 
#set ytics (2E5,3E5, 5E5, 1E6, 3E6)
#set xtics (1E-26,1E-25,1E-24,1E-23,1E-22,1E-21,1E-20,1E-19,1E-18,1E-17,1E-16,1E-15,1E-14,1E-13,1E-12,1E-11,1E-10,1E-9,1E-8)
#set logscale x
set logscale y
set format y "10^{%L}"
set xlabel "Photon Energy [eV]" offset 0,-4  font ",25"
set ylabel "Photon Intensity [m^{-2} s^{-1} sr^{-1} eV^{-1}]" offset -6,0 font ",28"
set title "Photon Intensity vs Photon Energy" font ",25"
set tmargin 10
set rmargin 10
set lmargin 20
set bmargin 20
set tics scale 1
set xtics font ", 22"
set ytics font ", 22"
set yrange [*:*]
set xrange [*:*]
set output 'intensity.pdf'
plot 'intensity.txt' lt 7 lw 4 lc 1 w lines

