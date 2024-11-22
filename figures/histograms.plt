reset

set terminal png

min = -200
max =  200

n = 50#floor(sqrt(10000))
width = (max - min)/n

hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set style fill solid 0.5 # fill style

set size sq
set xlabel "<{/Symbol F}>" font ",16"

do for [i=101:600]{
    set output sprintf("%.3d.png",i) 
p "../data/data.dat" every ::0::i u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "red" title sprintf("%d",i)
}
unset output
