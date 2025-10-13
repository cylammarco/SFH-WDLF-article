set terminal pngcairo dashed enhanced color size 640,640 font "DejaVuSans"

set style line 1 lt 1 lc rgb "red" lw 1
set style line 2 lt 3 lc rgb "black" lw 1 pt 5 ps 0.5

set xlabel ""
set xrange [5:18.5]
set xtics 2 out nomirror scale 0.5
set format x ''

set ylabel "Log {/Symbol F} [N Mag^{-1} pc^{-3}]"
set yrange [-6.41535251634884:-2.354584287560365]
set ytics 1 out nomirror scale 0.5
set mytics 2

set key off

set style line 20 lc rgb '#ddccdd' lt 1 lw 1.5
set style line 21 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 20, ls 21

f(s,n) = (s-n>0) ? (s-n) : 9E-18

set output 'nr_wdlf_bootstrap.png'

set lmargin 13
set bar 0.0

set multiplot

# Upper panel
set size 1.0,0.666
set origin 0.0, 0.333

plot 'MonteCarlo_wdlf_mean.txt' u 1:(log10($3)) w l ls 1 notitle,\
     'obsWDLF' u 1:(log10($3)):(log10(f($3,$4))):(log10($3+$4)) w yerrorbars ls 2 notitle

# Lower panel
set size 1.0,0.333
set origin 0.0, 0.0
set tmargin 0
set xlabel "M_{bol} [mag]"
set format x '%g'
set ylabel "Δ{/Symbol F} [N Mag^{-1} pc^{-3}]" offset 1,0
set yrange [-0.00013:0.00013]
set ytics 0.0001 out nomirror

plot '< paste obsWDLF MonteCarlo_wdlf_mean.txt' u 1:($3-$7):($3-$7-$4):($3-$7+$4) w yerrorbars ls 2 notitle

