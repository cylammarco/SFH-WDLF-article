set terminal pngcairo dashed enhanced color size 640,480 font "DejaVuSans"

# Shaded region style
set style line 1 lt 1 lc rgb "gray40" lw 1
# Line style
set style line 2 lt 1 lc rgb "black" lw 2

set xlabel "Lookback time [Gyr]"
set xrange [0:14.5]
set xtics 2 out nomirror
set xtics format "%g"

set ylabel "Star Formation Rate [N Gyr^{-1} pc^{-3}]"
set yrange [0:0.00327]

set ytics 0.001 out nomirror
set mytics 2

set key off

# Background grid lines
set style line 20 lc rgb '#ddccdd' lt 1 lw 1.5 
set style line 21 lc rgb '#ddccdd' lt 1 lw 0.5
set grid xtics mxtics ytics mytics back ls 20, ls 21

set output 'nr_sfr_bootstrap.png'

plot 'MonteCarlo_sfhPlotData.txt' u ($1/1E9):(($3+$5)*1E9):(($3-$5)*1E9) w filledcurves fillstyle transparent solid 0.5 ls 1 notitle,\
     'MonteCarlo_sfhPlotData.txt' u ($1/1E9):($3*1E9) w l ls 2 notitle

