plot 'M.gnu' u 1:(-$2):1:($1+1):(-$2):(-$2-1):3 \
    with boxxyerrorbars \
    fill solid border linecolor 'black' linecolor palette title 'Mass matrix sparsity'
