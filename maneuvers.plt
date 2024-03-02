# Plot shadowing orbit in position space, i.e. the points after correction
# maneuvers.

set datafile separator ','
set key autotitle columnhead	# use the first line as key titles
splot "maneuvers.csv" using 2:3:4 with points 
