set term pngcairo truecolor font 'Fira Sans Reqular,16' fontscale 1.8 size 2560,1600
set output 'plots/two.png' 
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

set title 'Fig. 1. Statistics of 1024 symbols and 0.5 dBm step.'
plot\
'two-pol.csv' u 1:2 w lines lw 4 lc 'grey' title ' PSE x-pol',\
'two-pol.csv' u 1:3 w lines lw 4 lc 'black' dt 2 title ' PSE y-pol',\
'two-pol-half.csv' u 1:4 w lines lw 4 lc 'skyblue' title ' SINO x-pol',\
'two-pol-half.csv' u 1:5 w lines lw 4 lc 'blue' dt 2 title ' SINO y-pol',\
'two-pol.csv' u 1:4 w lines lw 4 lc 'pink' title ' SINO (full) x-pol',\
'two-pol.csv' u 1:5 w lines lw 4 lc 'red' dt 2 title ' SINO (full) y-pol',\

set title 'Fig. 2. Same but with cubic splines smoothing.'
plot\
'two-pol.csv' u 1:2 w lines lw 4 lc 'grey' title ' PSE x-pol' smooth acsplines,\
'two-pol.csv' u 1:3 w lines lw 4 lc 'black' dt 2 title ' PSE y-pol' smooth acsplines,\
'two-pol-half.csv' u 1:4 w lines lw 4 lc 'skyblue' title ' SINO x-pol' smooth acsplines,\
'two-pol-half.csv' u 1:5 w lines lw 4 lc 'blue' dt 2 title ' SINO y-pol' smooth acsplines,\
'two-pol.csv' u 1:4 w lines lw 4 lc 'pink' title ' SINO (full) x-pol' smooth acsplines,\
'two-pol.csv' u 1:5 w lines lw 4 lc 'red' dt 2 title ' SINO (full) y-pol' smooth acsplines,\

unset multiplot
