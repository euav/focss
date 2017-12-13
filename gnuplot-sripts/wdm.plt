set term pngcairo truecolor font 'Fira Sans Reqular,16' fontscale 1.8 size 2560,1600
set output 'wdm-dbp-alt.png' 
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
'dbp-kappa-4k.txt' u 1:11 w lines lw 4 lc 'black' title ' No-WDM DBP 10',\
'ipr-sino-kappa-4k.txt' u 1:2 w lines lw 4 lc 'skyblue' title ' No-WDM IPR',\
'wdm-ipr.csv' u 1:2 w lines lw 4 lc 'blue' dt 2 title ' WDM IPR',\
'wdm-dbp-alt.csv' u 1:6 w lines lw 4 lc 'grey' dt 3 title ' false 5-channel WDM DBP 130',\
'wdm-dbp-alt.csv' u 1:5 w lines lw 4 lc 'grey' title ' 5-channel WDM DBP 130',\
'wdm-dbp-alt.csv' u 1:4 w lines lw 4 lc 'grey' dt 2 title ' alt 3-channel WDM DBP 130',\
'wdm-dbp-alt.csv' u 1:3 w lines lw 4 lc 'grey' title ' 3-channel WDM DBP 130',\
'wdm-dbp-alt.csv' u 1:2 w lines lw 4 lc 'grey' title '1-channel WDM DBP 130',\

set title 'Fig. 2. Same but with cubic splines smoothing.'
plot\
'dbp-kappa-4k.txt' u 1:11 w lines lw 4 lc 'black' title ' No-WDM DBP 10' smooth acsplines,\
'ipr-sino-kappa-4k.txt' u 1:2 w lines lw 4 lc 'skyblue' title ' No-WDM IPR' smooth acsplines,\
'wdm-ipr.csv' u 1:2 w lines lw 4 lc 'blue' dt 2 title ' IPR' smooth acsplines,\
'wdm-dbp-alt.csv' u 1:6 w lines lw 4 lc 'grey' dt 3 title ' false 5-channel WDM DBP 130' smooth acsplines,\
'wdm-dbp-alt.csv' u 1:5 w lines lw 4 lc 'grey' title ' 5-channel WDM DBP 130' smooth acsplines,\
'wdm-dbp-alt.csv' u 1:4 w lines lw 4 lc 'grey' dt 2 title ' alt 3-channel WDM DBP 130' smooth acsplines,\
'wdm-dbp-alt.csv' u 1:3 w lines lw 4 lc 'grey' title ' 3-channel WDM DBP 130' smooth acsplines,\
'wdm-dbp-alt.csv' u 1:2 w lines lw 4 lc 'grey' title '1-channel WDM DBP 130' smooth acsplines,\

unset multiplot
