function normalize_eigenvectors(output::Output, intervall, dimension, potential)
    for i in 1:output.n_eigenvalues
        density = output.eigenvectors[i].^2
    
        norm = sum(density) * ustrip(uconvert(potential.coordsUnit, intervall*potential.internalElemCoords))^dimension

        output.eigenvectors[i] = output.eigenvectors[i] ./ sqrt(norm)
    end
end