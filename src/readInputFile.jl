function readInputFile(inputFileName::String)

    inputFile = open(inputFileName, "r")

    lines = readlines(inputFile)

    lineElements = split.(lines)

    for line in lineElements
        
        line[2] != "=" && @error "Parsing error in inputfile -- second entry in a line has to be a \"=\""

        if("potential-file" == line[1])
            continue
        end
        @error "Keyword $(line[1]) not defined!"
        break
    end
    
end