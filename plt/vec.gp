meshfunction = ARG1
xgrid = ARG2
ygrid = ARG3
component = ARG4
if (xgrid eq "") { xgrid = 20 }
if (ygrid eq "") { ygrid = 20 }

set dgrid3d xgrid, ygrid, 4

set table "interpolated_data_x.gnu"
splot meshfunction u 1:2:4
unset table

set table "interpolated_data_y.gnu"
splot meshfunction u 1:2:5
unset table
unset dgrid3d

plot "< paste interpolated_data_x.gnu interpolated_data_x.gnu" \
    u 1:2:3:7:($3*$3+$7*$7) with vectors title 'vector field' linecolor palette
