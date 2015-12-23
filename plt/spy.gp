matrix_file = ARG1

plot matrix_file u 2:(-$1):2:($2+1):(-$1):(-$1-1):3 \
    with boxxyerrorbars \
    fill solid border linecolor 'black' linecolor palette title 'Mass matrix sparsity'
