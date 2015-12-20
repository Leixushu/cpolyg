meshfunction = ARG1
if (ARG2 eq "") { ARG2 = 20 }
if (ARG3 eq "") { ARG3 = 20 }

set dgrid3d ARG2, ARG3, 4
unset pm3d

set table "interpolated_data.gnu"
splot meshfunction

unset table
unset dgrid3d

set yrange [*:*]
stats 'interpolated_data.gnu' u 3 nooutput

splot 'interpolated_data.gnu' with pm3d title 'mesh function', \
      'edges.gnu' u 1:2:(STATS_min) with lines linetype rgb 'gray30' title 'mesh'
