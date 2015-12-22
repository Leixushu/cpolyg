meshfunction = ARG1
if (ARG2 eq "") { ARG2 = 20 }
if (ARG3 eq "") { ARG3 = 20 }
if (!exists("c")) { c = 0 }
component = 3 + c

set dgrid3d ARG2, ARG3, 4
unset pm3d

set table "interpolated_data.gnu"
splot meshfunction u 1:2:component

unset table
unset dgrid3d

set yrange [*:*]

if (!exists("meshz")) {
    stats 'interpolated_data.gnu' u 3 nooutput
    fn_mesh_z = STATS_min
} else {
    fn_mesh_z = meshz
}

splot 'interpolated_data.gnu' with pm3d title 'mesh function', \
      'edges.gnu' u 1:2:(fn_mesh_z) with lines linetype rgb 'gray30' title 'mesh'
