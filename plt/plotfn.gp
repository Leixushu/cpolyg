meshfunction = ARG1
if (ARG2 eq "") { ARG2 = 20 }
if (ARG3 eq "") { ARG3 = 20 }

set dgrid3d ARG2, ARG3, 4
splot meshfunction with pm3d title 'mesh function'
