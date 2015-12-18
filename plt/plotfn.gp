meshfunction = ARG1
xgrid = ARG2
ygrid = ARG3
if (xgrid eq "") { xgrid = 20 }
if (ygrid eq "") { ygrid = 20 }

set dgrid3d xgrid, ygrid, 4
splot meshfunction with pm3d title 'mesh function'
