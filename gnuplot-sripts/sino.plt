set term pngcairo truecolor font 'Fira Sans Reqular,16' fontscale 1.8 size 2560,1600
set output 'sino.png' 
set datafile separator ','

set xtics 2 
set xrange[-10:10]
set ytics 1
set yrange[12:25]
set grid lw 2

set key top left invert
set xlabel 'Average power [dBm]'
set ylabel 'Q^2-factor [dB]'

set multiplot layout 1, 2

set title 'Fig. 1. Statistics of 4096 symbols and 0.5 dBm step.'
plot\
'gold-final.txt' u 1:2 w lines lw 4 lc 'gray' title ' DBP 1:10',\
'gold-final.txt' u 1:3 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:4 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:5 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:6 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:7 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:8 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:9 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:10 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:11 w lines lw 4 lc 'grey' notitle,\
'gold-final.txt' u 1:12 w lines lw 4 lc 'black' dt 2 title ' DBP ideal',\
'sino11.txt' u 1:2 w lines lw 4 lc 'pink' title ' LMS',\
'sino11.txt' u 1:3 w lines lw 4 lc 'skyblue' title ' IPR',\
'sino11.txt' u 1:4 w lines lw 4 lc 'blue' dt 2 title ' NLPR',\
'sino11.txt' u 1:5 w lines lw 4 lc 'red' title ' SINO 13',\

set title 'Fig. 2. Same but with cubic splines smoothing.'
plot\
'gold-final.txt' u 1:2 w lines lw 4 lc 'gray' title ' DBP 1:10' smooth acsplines,\
'gold-final.txt' u 1:3 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:4 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:5 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:6 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:7 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:8 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:9 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:10 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:11 w lines lw 4 lc 'grey' notitle smooth acsplines,\
'gold-final.txt' u 1:12 w lines lw 4 lc 'black' dt 2 title ' DBP ideal' smooth acsplines,\
'sino11.txt' u 1:2 w lines lw 4 lc 'pink' title ' LMS' smooth acsplines,\
'sino11.txt' u 1:3 w lines lw 4 lc 'skyblue' title ' IPR' smooth acsplines,\
'sino11.txt' u 1:4 w lines lw 4 lc 'blue' dt 2 title ' NLPR' smooth acsplines,\
'sino11.txt' u 1:5 w lines lw 4 lc 'red' title ' SINO 13' smooth acsplines,\

 unset multiplot
