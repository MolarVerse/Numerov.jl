using SparseArrays
using LinearAlgebra

for x_1 in 1:n_datapoints[1]
	for stencil_index in 1:system.stencil

		super_matrix_index = x_1 - system.stencil ÷ 2 - 1 + stencil_index

		if super_matrix_index < 1 && system.periodic == false
			continue
		elseif super_matrix_index < 1
			super_matrix_index += n_datapoints[1]
		end

		matrix = zeros(n_datapoints[2], n_datapoints[2])

		for j in 1:system.stencil

			sub_matrix_index = j - system.stencil ÷ 2 - 1

			if sub_matrix_index == 0
				matrix += spdiagm(sub_matrix_index => ones(n_datapoints[2] - sub_matrix_index) * stencil[stencil_index, j])
			else
				matrix += spdiagm(sub_matrix_index => ones(n_datapoints[2] - abs(sub_matrix_index)) * stencil[stencil_index, j])
			end

			if system.periodic && sub_matrix_index != 0
                println(sub_matrix_index)
				matrix += spdiagm(sign(sub_matrix_index)*n_datapoints[2] - sub_matrix_index => ones(abs(sub_matrix_index)) * stencil[stencil_index, j])
				# matrix += spdiagm(-n_datapoints[2] + abs(sub_matrix_index) => ones(abs(sub_matrix_index)) * stencil[stencil_index, j])
			end
		end

		i_1 = 1 + (x_1 - 1) * n_datapoints[2]
		j_1 = x_1 * n_datapoints[2]

		i_2 = 1 + (super_matrix_index - 1) * n_datapoints[2]
		j_2 = (super_matrix_index) * n_datapoints[2]

        if(i_2 > n_datapoints[2]*n_datapoints[1])
            continue
        end

		laplace[i_1:j_1, i_2:j_2] += matrix #add system.laplace

	end
end

for i in 1:n_datapoints[1]*n_datapoints[2]
    for j in i:n_datapoints[1]*n_datapoints[2]
        laplace[j,i] = laplace[i,j]
    end
end
