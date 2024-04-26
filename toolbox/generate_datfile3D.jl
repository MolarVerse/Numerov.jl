using DelimitedFiles

function generate_datfile(inputfile, imag_inputfile, outputfile, state)

    infile = readdlm(inputfile)
    imag_infile = readdlm(imag_inputfile)

    outfile = open(outputfile, "w")

    state = parse(Int64, state)

    for i in 1:size(infile)[1]
        println(outfile, infile[i, 1], " ", infile[i, 2], " ", infile[i, 3], " ", infile[i, 5+state]^2 + imag_infile[i, 5+state]^2)
    end

    close(outfile)
end
    
function main()
    inputfile = ARGS[1]
    imag_inputfile = ARGS[2]
    outputfile = ARGS[3]
    if length(ARGS) == 4
        state = ARGS[4]
        generate_datfile(inputfile, imag_inputfile, outputfile, state)
    else
        generate_datfile(inputfile, imag_inputfile, outputfile, "0")
    end

end

main()