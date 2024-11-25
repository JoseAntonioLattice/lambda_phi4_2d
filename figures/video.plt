reset

set terminal png enhanced

min = -1.0
max =  1.0

L = 32

n = 50#floor(sqrt(10000))
width = (max - min)/n

hist(x,width) =  width*(floor((x-min)/width)+0.5) + min
set boxwidth width
set size sq
set cbrange[-1:1]
set palette rgb 33,13,10
j = 0

do for [i=0:5000:10]{

    set output sprintf("images/all/%.3d.png",j)
    set multiplot layout 2,2 title sprintf("λ = 0.05, m^2 = -0.075, L = %d, measurements = %d",L,j*10+1) font ",18"

    unset xtics
    unset ytics
    unset xlabel
    set xrange[-0.5:L-0.5]
    set yrange[-0.5:L-0.5]
    p "../data/thermalized_configurations.dat" i i matrix with image notitle
    unset xrange
    unset yrange

    set ytics
    set xlabel 'measurements'
    set title "action/V"
    p "../data/data.dat" every ::0::i u 0:1 w l notitle
    set title "|magnetization/V|"
    p "../data/data.dat" every ::0::i u 0:3 w l notitle
    set xrange[-2:2]
    unset title
    set xtics
    set xlabel 'magnetization/V'
    p "../data/data.dat" every ::0::i u (hist($2,width)):(1.0) smooth freq w boxes lc rgb "red" notitle #sprintf("%d",)
    unset xrange
    j = j + 1
unset multiplot
}
unset output
