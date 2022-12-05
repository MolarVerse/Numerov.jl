using DelimitedFiles
using Printf

function generate_cubefile(filename::String, x::Int, y::Int, z::Int; outputfilename="mycube.cube", column=4)
    data = readdlm(filename)

    step = data[2,3] - data[1,3]

    file = open(outputfilename, "w")

    println(file, "my cube file")
    println(file, "###")
    println(file, "1        0.0         0.0         0.0")
    println(file, x, "     ", step, "   0.0    0.0")
    println(file, y, "      0.0       ", step, "    0.0")
    println(file, z, "      0.0         0.0      ", step)
    println(file, "1 0.0 0.0 0.0")

    count = 0
    for i in 1:size(data)[1]

        count += 1

        @printf(file, "%.6e ", data[i,column])

        if count % 6 == 0 || count % z == 0
            println(file)
        end

        if count % z == 0
            count = 0
        end
    end

    close(file)
        
end