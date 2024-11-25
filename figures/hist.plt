reset
min = -1.0
max =  1.0

n = 50
width = (max - min)/n

hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set style fill solid 0.5 # fill style

p "../data/data.dat" u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "red"
