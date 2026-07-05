function readInputFile(inputFileName::String)

    ##################################################
    #                                                #
    # prepare input file for processing to read data #
    #                                                #
    ##################################################

    inputFile = open(inputFileName, "r")

    lines = readlines(inputFile)
    
    lines = getindex.(split.(lines, "#"), 1) # removing all comments starting with #
    
    filter!(x -> !isempty(strip(x)), lines) # removing all blank lines from input

    lineElements = split.(lines) # values are kept verbatim - only the keyword is lowercased below

    ##############################################################################
    #                                                                            #
    # read line per line and save values assigned to keywords in inputDictionary #
    #                                                                            #
    ##############################################################################

    for line in lineElements

        keyFound = false

        length(line) < 3   && throw(ArgumentError("There are too few entries in line $(line)"))
        line[2]     != "=" && throw(ArgumentError("Parsing error in inputfile -- second entry in a line has to be a \"=\""))

        keyword = lowercase(line[1])

        for (key, _) in inputDictionary
            if key == keyword && keyFound == false
                inputDictionary[key] = join(line[3:end], " ")
                keyFound = true
            elseif key == keyword
                throw(ArgumentError("You have defined the keyword $(keyword) multiple times"))
            end
        end

        keyFound == false && throw(ArgumentError("Keyword $(line[1]) not defined!"))

    end

end