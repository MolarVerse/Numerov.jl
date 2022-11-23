using DelimitedFiles
using LinearAlgebra

function bandStructure(args::Vector{String})
    
    if length(args) < 3
        @error ""
    end

    dimensions = parse.(Int64, args[2:end])

    data = readdlm(args[1]; skipstart=1)

    spacing    = zeros(2)
    spacing[1] = data[2,2] - data[1,2]
    spacing[2] = data[dimensions[2]+1,1] - data[dimensions[2],1]
    
    for i in 1:dimensions[2]
        print((i-1)*spacing[1], "  ")
        for j in length(dimensions)+1:length(data[1,:])
            print(data[i, j], " ")
        end
        println()
    end
 
    for i in 2:dimensions[1]
        print((i-1)*spacing[2]+(dimensions[1]-1)*spacing[1], "  ")
        for j in length(dimensions)+1:length(data[1,:])
            print(data[i*dimensions[2], j], " ")
        end
        println()
    end

    index = 1
    for i in dimensions[1]-1:-1:1
        print(index*norm(spacing) + (dimensions[2]-1)*spacing[2]+(dimensions[1]-1)*spacing[1], "  ")
        for j in length(dimensions)+1:length(data[1,:])
            print(data[i*dimensions[2]+i, j], " ")
        end
        println()
        index += 1
    end
end