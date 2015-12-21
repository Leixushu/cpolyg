meshfunction = ARG1
xgrid = ARG2
ygrid = ARG3
component = ARG4
if (xgrid eq "") { xgrid = 20 }
if (ygrid eq "") { ygrid = 20 }
if (!exists("c")) { c = 0 }

component = c+3

set dgrid3d xgrid, ygrid, 4
splot meshfunction u 1:2:component with pm3d palette title 'mesh function'
