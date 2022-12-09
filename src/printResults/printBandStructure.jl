function printBandStructure(potential::Potential, k_points)

    file = open("bandstructure.dat", "w")

    data = readdlm("eigenvalues.dat"; skipstart=1)

    a = copy(k_points)
    b = copy(k_points)

    pushfirst!(a, Tuple(zeros(length(k_points[1]))))
    push!(b, Tuple(zeros(length(k_points[1]))))

    diff = (norm.([b[i] .- a[i] for i in eachindex(a)]))[1:end-1]

    brioullin_path = 0.0

    for (i, spacing) in enumerate(diff)
        brioullin_path += spacing
        @printf(file, "%lf ", brioullin_path)
        for j in length(k_points[1])+1:length(data[1,:])
            @printf(file, "%20.14lf ", data[i,j])
        end
        @printf(file, "\n")
        if i == length(diff) - potential.n_kpoints && potential.dimension == 3
            break
        end
    end

    for (i, spacing) in enumerate(diff)
        if i <= length(diff) - potential.n_kpoints || potential.dimension != 3
            continue
        end
        @printf(file, "%lf ", brioullin_path)
        for j in length(k_points[1])+1:length(data[1,:])
            @printf(file, "%20.14lf ", data[i,j])
        end
        @printf(file, "\n")
        brioullin_path += spacing
    end

    close(file)    
end