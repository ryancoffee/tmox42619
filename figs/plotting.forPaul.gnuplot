#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 8    last modified 2019-12-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2019
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal qt 0 font "Sans,9"
# set output
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "% h" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics back
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key fixed right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 0.25
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view map scale 1
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels 5
set cntrparam levels auto
set cntrparam firstlinetype 0 unsorted
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data image
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set title "" 
set title  font "" textcolor lt -1 norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" textcolor lt -1 norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 950.000 : 1200.00 ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ 1000.00 : 2000.00 ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
file(p) = sprintf('/media/coffee/9C33-6BBD/ascii/port_%i.hist.dat',p)
GNUTERM = "qt"
file = "/media/coffee/9C33-6BBD/ascii/port_0.hist.dat"
x = 0.0
R = .4
unset xtics
unset ytics
unset colorbox
set cbrange [0:5000]
set yrange [1400:2000]
set xrange [950:1200]
X(th) = 0.5+R*sin((2*pi*(th))/16)-.11
Y(th) = 0.5+R*cos((2*pi*(th))/16)-.125
## Last datafile plotted: "/media/coffee/9C33-6BBD/ascii/port_14.hist.dat"
set term png size 1000,1000
set output 'figs/plotting.forEliane.png'
set multiplot
set size .21,.25
positions = '0 1 4 5 12 13 14 15'
retardations = '300V 75V 75V 300V 75V 0V 75V 300V'
do for [n=1:words(positions)] {
	set title word(retardations,n) offset 0,-1
	set origin X(int(word(positions,n))),Y(int(word(positions,n)))
	splot file(int(word(positions,n))) mat every 2:2 notitle
}
set colorbox
set label 1 'O-KLL' at screen .6,.59 front textcolor 'white'
set label 2 'N-KLL' at screen .6,.61 front textcolor 'white'
set label 3 'N-photo' at screen .6,.7 front textcolor 'white'
set label 11 'O-KLL' at screen .6,.33 front textcolor 'white'
set label 12 'N-KLL' at screen .6,.35 front textcolor 'white'
set label 13 'N-photo' at screen .6,.4 front textcolor 'white'
set label 21 'O-KLL' at screen .6,.115 front textcolor 'white'
set label 22 'N-KLL' at screen .6,.16 front textcolor 'white'
set ylabel 'log2(ToF) [arb]'
set xlabel 'vls [arb]' offset 0,2
set size .5,.35
set origin .25,0
set title '300V' offset 0,-1
splot file(5) mat notitle
set origin .25,.25
set title '0V' offset 0,-1
splot file(13) mat notitle
set origin .25,.5
set title '75V' offset 0,-1
set arrow 1 from screen .33,.78 to screen .59,.81 front nohead
set arrow 2 from screen .68,.78 to screen .71,.81 front nohead 
set arrow 3 from screen .19,.71 to screen .32,.53 front nohead
set arrow 4 from screen .19,.59 to screen .32,.33 front nohead
set arrow 5 from screen .8,.4 to screen .68,.28 front nohead
set arrow 6 from screen .8,.28 to screen .68,.07 front nohead
splot file(1) mat notitle
unset multiplot

#    EOF
