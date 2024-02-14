#set term tikz
#set out "orbit.tex"

set xlabel '$x$'
set ylabel '$y$'

plot "orbit.dat" u 2:4 w l

unset out
unset term
