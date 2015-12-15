meshfunction = ARG1
if (ARG2 eq "") ARG2 = 20
if (ARG3 eq "") ARG3 = 20

set dgrid3d ARG2, ARG3, 4
unset pm3d

set table "interpolated_data.dat"
splot meshfunction


unset table
unset dgrid3d

stats 'interpolated_data.dat' u 3 nooutput

splot 'interpolated_data.dat' with pm3d title 'mesh function', \
      'plt/edges.gnu' u 1:2:(STATS_min) with lines linetype rgb 'gray30' title 'mesh'
