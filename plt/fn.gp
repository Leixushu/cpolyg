meshfunction = ARG1
xgrid = ARG2
ygrid = ARG3
component = ARG4
if (xgrid eq "") { xgrid = 20 }
if (ygrid eq "") { ygrid = 20 }
if (!exists("c")) { c = 0 }

component = c+3

set dgrid3d xgrid, ygrid, 4

if (exists("mesh")) {
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
} else {
    splot meshfunction u 1:2:component with pm3d palette title 'mesh function'
}

