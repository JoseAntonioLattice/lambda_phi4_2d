reset

set terminal png enhanced

min = -200
max =  200

n = 50#floor(sqrt(10000))
width = (max - min)/n

hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set size sq
set cbrange[-1:1]
set palette rgb 33,13,10
j = 0
do for [i=0:1000:10]{
    
    set output sprintf("images/all/%.3d.png",j)
    set multiplot layout 2,2 title sprintf("λ = 1, m^2 = 1, sweeps %d",(j+1)*10) font ",20"
 
    unset xtics
    unset ytics
    set xrange[-0.5:149.5]
    set yrange[-0.5:149.5]
    p "../data/thermalized_configurations.dat" i i matrix with image notitle
    unset xrange 
    unset yrange
    set xtics
    set ytics 
    p "../data/data.dat" every ::0::i u 0:($1/150**2) w l t "action/V"
    p "../data/data.dat" every ::0::i u 0:($2/150**2) w l t "magnetization/V"
    p "../data/data.dat" every ::0::i u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "red" notitle #sprintf("%d",)
    j = j + 1
unset multiplot
}
unset output
