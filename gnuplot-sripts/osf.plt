set term pngcairo truecolor font 'Fira Sans Reqular,16' fontscale 1.8 size 2560,1600
set output 'oversampling.png' 
set datafile separator ','

set xtics 2 
set xrange[-10:10]
set ytics 1
set yrange[0:28]
set grid lw 2

set key top left invert
set xlabel 'Average power [dBm]'
set ylabel 'Q^2-factor [dB]'

set multiplot layout 1, 2

set title 'Fig. 1. Statistics of 4096 symbols and 0.5 dBm step.'
plot\
'oversampling.csv' u 1:2 w lines lw 4 lc 'grey' dt 2 title ' OSF 2',\
'oversampling.csv' u 1:3 w lines lw 4 lc 'grey' title ' OSF 16',\
'oversampling.csv' u 1:4 w lines lw 4 lc 'black' title ' OSF 128',\

set title 'Fig. 2. Same but with cubic splines smoothing.'
plot\
'oversampling.csv' u 1:2 w lines lw 4 lc 'grey' dt 2 title ' OSF 2' smooth acsplines,\
'oversampling.csv' u 1:3 w lines lw 4 lc 'grey' title ' OSF 16' smooth acsplines,\
'oversampling.csv' u 1:4 w lines lw 4 lc 'black' title ' OSF 128' smooth acsplines,\

unset multiplot
