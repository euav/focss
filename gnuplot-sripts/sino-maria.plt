set term pngcairo truecolor font 'Fira Sans Reqular,16' fontscale 1.8 size 2560,1600
set output 'sino-score.png' 
set datafile separator ','

set xtics 2 
set xrange[-10:10]
set ytics 1
set yrange[12:28]
set grid lw 2

set key top left invert
set xlabel 'Average power [dBm]'
set ylabel 'Q^2-factor [dB]'

set multiplot layout 1, 2

set title 'Fig. 1. Statistics of 4096 symbols and 0.5 dBm step.'
plot\
'dbp-kappa-4k.txt' u 1:2 w lines lw 4 lc 'gray' title ' DBP 1:10',\
'dbp-kappa-4k.txt' u 1:3 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:4 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:5 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:6 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:7 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:8 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:9 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:10 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:11 w lines lw 4 lc 'grey' notitle,\
'dbp-kappa-4k.txt' u 1:12 w lines lw 4 lc 'black' title ' DBP ideal',\
'sino_score.csv' u 1:3 w lines lw 4 lc 'blue' title ' PSE',\
'sino_score.csv' u 1:5 w lines lw 4 lc 'red' title ' SINO',\

set title 'Fig. 2. Same but with cubic splines smoothing.'
plot\
'dbp-kappa-4k.txt' u 1:2 w lines lw 4 lc 'gray' title ' DBP 1:10' smooth acsplines,\
'dbp-kappa-4k.txt' u 1:3 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:4 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:5 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:6 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:7 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:8 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:9 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:10 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:11 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'dbp-kappa-4k.txt' u 1:12 w lines lw 4 lc 'black' title ' DBP ideal' smooth acsplines,\
'sino_score.csv' u 1:3 w lines lw 4 lc 'blue' title ' PSE' smooth acsplines,\
'sino_score.csv' u 1:5 w lines lw 4 lc 'red' title ' SINO' smooth acsplines,\

unset multiplot
