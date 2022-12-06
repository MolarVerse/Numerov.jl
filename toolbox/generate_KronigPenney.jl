function generate_KronigPenney(N, D)

    file = open("potential.dat", "w")

    if D == 1
        for i in -4*N:4*N

            x = i / (4*N)

            print(file, x, "   ")

            if -N <= i <= N
                println(file, "0.0")
            else
                println(file, "1.0")
            end
        end
    end
    
end

generate_KronigPenney(parse(Int, ARGS[1]), parse(Int, ARGS[2]))