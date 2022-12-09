function build_1d_stencil(system, n_datapoints, stencil, stencil_size)
    matrix = zeros(n_datapoints, n_datapoints)
    
    for j in 1:stencil_size

        sub_matrix_index = j - stencil_size÷2 - 1 

        matrix += spdiagm(sub_matrix_index => ones(n_datapoints - abs(sub_matrix_index)) * stencil[j])

        if system.periodic[end] && sub_matrix_index != 0
            matrix += spdiagm(sign(sub_matrix_index)*n_datapoints - sub_matrix_index => ones(abs(sub_matrix_index)) * stencil[j])
        end
    end

    return matrix
end

function build_2d_stencil(system, n_datapoints, stencil, stencil_size)

    total_points = prod(n_datapoints)

    matrix = spzeros(total_points, total_points)

    for i in 1:stencil_size

        sub_matrix_index = i - stencil_size÷2 - 1

        matrix_1d = build_1d_stencil(system, n_datapoints[2], stencil[i,:], stencil_size)

        for j in 1:n_datapoints[1]

            cell_index = j + sub_matrix_index

            range_i = ((j-1)*n_datapoints[2] + 1):j*n_datapoints[2]
            range_j = ((cell_index-1)*n_datapoints[2] + 1):cell_index*n_datapoints[2]

            if cell_index < 1 || cell_index > n_datapoints[1]
                if system.periodic[end-1]
                    range_j = (range_j[1] - sign(cell_index-1)*total_points):(range_j[end] - sign(cell_index-1)*total_points)
                else
                    continue
                end
            end

            matrix[range_i, range_j] = matrix_1d
        end
    end

    return matrix
end

function build_3d_stencil(system, n_datapoints, stencil, stencil_size)

    total_points = prod(n_datapoints)

    matrix = spzeros(total_points, total_points)

    for i in 1:stencil_size

        sub_matrix_index = i - stencil_size÷2 - 1

        matrix_2d = build_2d_stencil(system, n_datapoints[2:3], stencil[:,:,i], stencil_size)

        for j in 1:n_datapoints[1]

            cell_index = j + sub_matrix_index

            range_i = ((j-1)*prod(n_datapoints[2:3]) + 1):j*prod(n_datapoints[2:3])
            range_j = ((cell_index-1)*prod(n_datapoints[2:3]) + 1):cell_index*prod(n_datapoints[2:3])

            if cell_index < 1 || cell_index > n_datapoints[1]
                if system.periodic[end-2]
                    range_j = (range_j[1] - sign(cell_index-1)*total_points):(range_j[end] - sign(cell_index-1)*total_points)
                else
                    continue
                end
            end

            matrix[range_i, range_j] = matrix_2d
        end
    end

    return matrix
end