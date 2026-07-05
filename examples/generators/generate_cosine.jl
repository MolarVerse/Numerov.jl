function generate_cosine(N, D)

    file = open("potential.dat", "w")

    if D == 1
        for i in 0:N-1

            x = -π + i / N *(2π)

            println(file, x, "   ", cos(x) + 1)

        end
    elseif D == 2
        for i in 0:N-1
            for j in 0:N-1

                x = -π + i / N *(2π)
                y = -π + j / N *(2π)

                println(file, x, "   ", y, "   ", cos(x) + cos(y) + 2)
            end
        end
    elseif D == 3
        for i in 0:N-1
            for j in 0:N-1
                for k in 0:N-1

                    x = -π + i / N *(2π)
                    y = -π + j / N *(2π)
                    z = -π + k / N *(2π)

                    println(file, x, "   ", y, "   ", z, "   ", cos(x) + cos(y) + cos(z) + 3)
                end
            end
        end
    end

    close(file)
end

generate_cosine(parse(Int, ARGS[1]), parse(Int, ARGS[2]))