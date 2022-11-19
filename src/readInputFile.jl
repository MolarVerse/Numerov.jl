function readInputFile(inputFileName::String)

    inputFile = open(inputFileName, "r")

    lines = readlines(inputFile)
    
    # removing all comments starting with #
    lines = getindex.(split.(lines, "#"), 1)

    # removing all blank lines from input
    filter!(x -> !isempty(strip(x)), lines)

    lineElements = split.(lowercase.(lines))

    for line in lineElements

        keyFound = false
        
        #length(line) > 3   && (@error "There are to many entries in line $(line)"; exit())
        length(line) < 3   && (@error "There are to few entries in line $(line)"; exit())
        line[2]     != "=" && (@error "Parsing error in inputfile -- second entry in a line has to be a \"=\"" ; exit())

        for (key, _) in inputDictionary
            if key == line[1] && keyFound == false
                inputDictionary[key] = join(line[3:end], " ")
                keyFound = true
            elseif key == line[1]
                @error "You have defined the keyword $(line[1]) multiple times" ; exit()
            end
        end

        keyFound == false && (@error "Keyword $(line[1]) not defined!" ; exit())

    end

end