startIdx = int(ARG1)
totalNum = int(ARG2)
iterStep = ARG3

if (iterStep == "") {
    iterStep = 1
} else {
    iterStep = int(iterStep)
}

do for [i=startIdx:totalNum:iterStep] {
    plot 'edges.gnu' every :::0::i with lines
    
    if (i < totalNum) {
        pause mouse any
    }
}
