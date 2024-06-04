reset
set terminal postscript eps size 10,9 enhanced color font 'Times-Roman,55'
set pm3d map
#set pm3d interpolate 2,2
set output 'Figura1_exc.eps'
set xrange [0.0:10.0]
set yrange [0:10.0]
set xtics 5 offset 0.0,0.75
set ytics 5 offset 1.0
set ylabel "{g_{/[Times-Roman] ie} }  " offset 1.0,0.0 font 'Times-Italic,55'
set palette defined (1 "black", 1 "black", 2.0 "blue", 3 "red", 6 "yellow")
set palette defined (1 "black", 3 "red", 4.0 "orange", 5 "yellow")
set multiplot layout 2,2
###########
set origin 0.05,0.43
set size 0.5,0.6
unset xtics
set cbtics 0.5 offset -1.0,0.0
set cblabel "R1" offset -6.0,4.7 rotate by 0 left font 'Times-Italic,55'
set palette defined (1 "black", 3 "red", 4.0 "orange", 5 "yellow")
set origin 0.05,0.43
set size 0.5,0.6
set cbrange [0.0:0.5]
set cbtics 0.5 offset -1.0,0.0
set ytics 5.0 offset 1.0
set cblabel "R" offset -6.4,4.7 rotate by 0 left font 'Times-Italic,55'
splot "POK_CV_F11.dat" u 2:5:9 t"
#####################################3
set palette defined (0 "black", 1 "white", 2.0 "yellow")
set palette defined (0 "blue", 1 "white", 2.0 "red")
#set palette defined (0 "blue", 0.2 "white",0.7 "white",  1 "red")
set origin 0.5,0.43
set size 0.5,0.6
set cbrange [0:1.0]
unset ytics
unset xtics
unset ylabel
unset ylabel
#set xtics 100 offset 0.0,0.75
set cbtics 0.5 offset -1.0,0.0
set cblabel "CV" offset -7.0,4.7 rotate by 0 left font 'Times-Roman,55'
splot "POK_CV_F11.dat" u 2:5:8 t"
#############################
#set ytics  5.0 offset 1.0
set palette defined (1 "black", 3 "red", 4.0 "orange", 5 "yellow")
set palette defined (1 "blue", 4 "green", 4.5 "yellow", 6 "red")
set origin 0.05,-0.02
set size 0.5,0.6
set cbrange [1.5:2.5]
#unset xlabel
#unset ytics
#unset xtics
#unset ylabel
set ytics 5.0 offset 1.0
set ylabel "{g_{/[Times-Roman] ie}}  " offset 1.0,0.0 font 'Times-Italic,55'
set cblabel "F" offset -6.3,4.75 rotate by 0 left font 'Times-Italic,55'
set xtics 100 offset 0.0,0.5
set cbtics 0.5
set cbtics ("1.5" 1.5, ">2.5" 2.5 )
set xtics 5
set xlabel "g_{/[Times-Roman] ei} {/[Times-Roman](nS)} " offset 0.0,0.8 font 'Times-Italic,55'
splot "POK_CV_F11.dat" u 2:5:7 t"
#####################################3
set origin 0.5,-0.02
set size 0.5,0.6
unset ytics
unset ylabel
set palette defined ( 0 0 0 0, 0.1667 0 0 1, 0.5 0 1 0,\
     0.6667 1 0 0, 1 1 1 1 )
set palette defined (0 "white",1.5 "blue", 2 "black",3 "green", 3.2 "yellow", 5 "red")
set palette defined (0 "white",1.8 "blue", 2.5 "black", 3.9 "green", 4.3 "yellow", 4.7 "orange", 5.2 "red")
set cbrange [-20.0:15.0]
set cbtics 5.0 offset -1.0,0.0
set xtics 5.0
set cblabel "I_{/[Times-Roman] s}" offset -4.5,4.8 rotate by 0 left font 'Times-Italic,55'
splot "POK_CV_F11.dat" u 2:5:10 t"
######################################
unset multiplot
