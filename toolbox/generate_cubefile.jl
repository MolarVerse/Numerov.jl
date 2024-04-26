using DelimitedFiles

function gen_cubfile(input_file, output_file, npoints)

	infile = readdlm(input_file)

	outfile = open(output_file, "w")

	dimension = infile[2, 3] - infile[1, 3]
	boxsize = dimension * npoints / 2

	println(outfile, "CUBE FILE")
	println(outfile, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z")
	println(outfile, "   1   0.000000   0.000000   0.000000")
	println(outfile, "   $npoints   $dimension   0.000000   0.000000")
	println(outfile, "   $npoints   0.000000   $dimension   0.000000")
	println(outfile, "   $npoints   0.000000   0.000000   $dimension")
	println(outfile, "1 0.0 $boxsize $boxsize $boxsize")

	counter = 1
	for i in 1:npoints
		for j in 1:npoints
			for k in 1:npoints
				print(outfile, infile[counter, 4])
				if k % 6 == 0
					println(outfile)
				else
					print(outfile, " ")
				end
				counter += 1
			end
			println(outfile)
		end
	end

	close(outfile)
end

function main()
	input_file = ARGS[1]
	output_file = ARGS[2]
	npoints = parse(Int64, ARGS[3])
	gen_cubfile(input_file, output_file, npoints)
end

main()

