function generate_KronigPenney(N, D)

    file = open("potential.dat", "w")

    if D == 1
        for i in -4*N:4*N-1e-10

            x = i / (4*N)

            print(file, x, "   ")

            if -N <= i <= N
                println(file, "0.0")
            else
                println(file, "1.0")
            end
        end
    elseif D == 2
        for j in -4*N:4*N-1e-10
            for i in -4*N:4*N-1e-10

                x = j / (4*N)
                y = i / (4*N)

                print(file, x, "   ", y, "   ")

                if -N <= i <= N && -N <= j <= N
                    println(file, "0.0")
                else
                    println(file, "1.0")
                end
            end
        end
    elseif D == 3
        for k in -4*N:4*N-1e-10
            for j in -4*N:4*N-1e-10
                for i in -4*N:4*N-1e-10

                    x = k / (4*N)
                    y = j / (4*N)
                    z = i / (4*N)

                    print(file, x, "   ", y, "   ", z, "   ")

                    if -N <= i <= N && -N <= j <= N && -N <= k <= N
                        println(file, "0.0")
                    else
                        println(file, "1.0")
                    end
                end
            end
        end
    end

    close(file)
end

generate_KronigPenney(parse(Int, ARGS[1]), parse(Int, ARGS[2]))