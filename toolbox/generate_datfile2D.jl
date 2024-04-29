using DelimitedFiles

function generate_datfile(inputfile, imag_inputfile, outputfile, npoints, state)

    infile = readdlm(inputfile)
    imag_infile = readdlm(imag_inputfile)

    outfile = open(outputfile, "w")

    state = parse(Int64, state)
    npoints = parse(Int64, npoints)

    for i in 1:size(infile)[1]
        println(outfile, infile[i, 1], " ", infile[i, 2], " ", infile[i, 4+state]^2 + imag_infile[i, 4+state]^2)
        if i % npoints == 0
            println(outfile)
        end
    end

    close(outfile)
end
    
function main()
    inputfile = ARGS[1]
    imag_inputfile = ARGS[2]
    outputfile = ARGS[3]
    npoints = ARGS[4]
    if length(ARGS) == 5
        state = ARGS[5]
        generate_datfile(inputfile, imag_inputfile, outputfile, npoints, state)
    else
        generate_datfile(inputfile, imag_inputfile, outputfile, npoints, "0")
    end

end

main()