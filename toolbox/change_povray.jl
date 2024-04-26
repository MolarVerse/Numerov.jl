function change_povray(inputfile, outputfile)
	# Read the input file
	input = open(inputfile, "r")
	output = open(outputfile, "w")

	lines = readlines(input)

	centerBoxDone = false
	i = 1

	progress = 0
	progress_old = -1

	while i <= length(lines)

		if (progress > progress_old)
			println("Progress: ", progress * 10, "%")
			progress_old = progress
		else
			progress = Int(floor(10 * i / length(lines)))
		end

		if !startswith(lines[i], "// MoleculeID:")
			println(output, lines[i])
			i += 1
			continue
		end

		println(output, lines[i])
		i += 1

		if !centerBoxDone
			lines[i] = "#declare VMD_line_width=0.006;"
		else
			lines[i] = "#declare VMD_line_width=0.002;"
		end
		println(output, lines[i])

		i += 1

		for j in 1:3
			index = findfirst("rgbt", lines[i])
			lines[i] = lines[i][1:index[1]-1] * "rgbt<0.000,0.000,0.000,0.000>)"
			println(output, lines[i])
			i += 1
		end

		if !centerBoxDone
			lines[i] = "#declare VMD_line_width=0.006;"
			centerBoxDone = true
		else
			lines[i] = "#declare VMD_line_width=0.002;"
		end

		println(output, lines[i])
		i += 1

		for j in 1:9
			index = findfirst("rgbt", lines[i])
			lines[i] = lines[i][1:index[1]-1] * "rgbt<0.000,0.000,0.000,0.000>)"
			println(output, lines[i])
			i += 1
		end

	end

	println("Progress: 100%")

	close(output)
end

function main()
	inputfile = ARGS[1]
	outputfile = ARGS[2]
	change_povray(inputfile, outputfile)
end

main()
