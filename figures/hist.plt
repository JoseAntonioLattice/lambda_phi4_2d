reset
min = -2.0
max =  2.0

n = 50
width = 0.05#(max - min)/n

hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set style fill solid 0.5 # fill style

lmb = 0.05
msq = -0.11

set title sprintf("{/Symbol l} = %1.2f, m^2 = %2.3f", lmb, msq) font ",20"
set xlabel "magnetización" font ",16"
p "../data/data.dat" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "red" notitle
