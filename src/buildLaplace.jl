function buildLaplace(system::System1D)

    matrix = zeros(system.n_datapoints, system.n_datapoints)

    if system.stencil == 3

        coefficients = [1,-2,1]

        for i in 2:system.n_datapoints-1
            matrix[i,i-1:i+1] = coefficients
        end
        matrix[1,1:2] = coefficients[end-1:end]
        matrix[end,end-1:end] = coefficients[1:2]
    
    elseif system.stencil == 5

        coefficients = [-1/12, 4/3, -5/2, 4/3, -1/12]

    system.laplace = sparse(matrix)
    
end