if (ARG2 eq "") { ARG2 = 20 }
if (ARG3 eq "") { ARG3 = 20 }
if (ARG4 eq "") { ARG4 = 1 }

if (ARG1 ne "")
{
    do for [i=0:ARG1:ARG4] {
        call "plotfn.gp" sprintf("plt/u%d.gnu", i) 30 30 ; pause 0.01
    }
}
