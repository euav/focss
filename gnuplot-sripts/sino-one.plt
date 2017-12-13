set term pngcairo truecolor font 'Fira Sans Reqular,16' fontscale 1.8 size 2560,1600
set output 'results/plots/sino-one.png' 
set datafile separator ','

set xtics 2 
set xrange[-10:10]
set ytics 1
set yrange[10:25]
set grid lw 2

set key top left invert
set xlabel 'Average power [dBm]'
set ylabel 'Q^2-factor [dB]'

set multiplot layout 1, 2

set title 'Fig. 1. Statistics of 2^1^4 = 16384 symbols and 0.5 dBm step.'
plot\
'dbp-one-1.csv' u 1:3 w lines lw 4 lc 'gray' title ' DBP 1:10',\
'dbp-one-1.csv' u 1:4 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:5 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:6 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:7 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:8 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:9 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:10 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:11 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:12 w lines lw 4 lc 'grey' notitle,\
'dbp-one-1.csv' u 1:2 w lines lw 4 lc 'black' title ' DBP Ideal',\
'sino-one.csv' u 1:2 w lines lw 4 lc 'skyblue' dt 2 title ' PSE',\
'sino-one.csv' u 1:3 w lines lw 4 lc 'orange' title ' SINO',\


set title 'Fig. 2. Same but with cubic splines smoothing.'
plot\
'dbp-one-1.csv' u 1:3 w lines lw 4 lc 'gray' title ' DBP 1:10' smooth acsplines,\
'dbp-one-1.csv' u 1:4 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:5 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:6 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:7 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:8 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:9 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:10 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:11 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:12 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-one-1.csv' u 1:2 w lines lw 4 lc 'black' title ' DBP Ideal' smooth acsplines,\
'sino-one.csv' u 1:2 w lines lw 4 lc 'skyblue' dt 2 title ' PSE' smooth acsplines,\
'sino-one.csv' u 1:3 w lines lw 4 lc 'orange' title ' SINO' smooth acsplines,\


unset multiplot
