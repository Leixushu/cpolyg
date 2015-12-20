final = ARG1
xgrid = ARG2
ygrid = ARG3
step = ARG4

if (final eq "") { final = 0 }
if (xgrid eq "") { xgrid = 20 }
if (ygrid eq "") { ygrid = 20 }
if (step eq "") { step = 1 }
if (!exists("c")) { c = 0 }
if (!exists("start")) { start = 0 }

if (ARG1 ne "") {
    do for [i=start:final:step] {
        #print sprintf("Plot number %d", i);
        call "fn.gp" sprintf("u%d.gnu", i) (xgrid) (ygrid) (c);
        pause 0.02;
    }
}
