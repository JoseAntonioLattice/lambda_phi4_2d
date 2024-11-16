reset

set terminal png

set xrange [-0.5:149.5]
set yrange [-0.5:149.5]

set cbrange [-1:1]

unset xtics
unset ytics

set size sq
set palette rgb 33,13,10
do for [i=0:1000]{
  set output sprintf("images/%.4d.png",i)
  p "../data/configurations.dat" i i matrix with image
}
unset output
