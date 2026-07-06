"""
numerov - solve the Schrödinger equation using the Numerov method.

# Intro

Runs a full Numerov calculation: reads the input file, solves the 1D, 2D or 3D
time-independent Schrödinger equation on the grid potential defined therein and
writes the resulting eigenvalues, eigenvectors and frequencies - and, if
requested, the band structure - to `.dat` files together with a log file
(default: `Numerov.out`). All output files are written to the current working
directory; an existing `eigenvalues.dat` file is deleted before the new
eigenvalues are written.

# Args

- `inputFileName`: The name of the input file.
"""
Comonicon.@main function numerov(inputFileName::String)


	#########################################################################
	#                                                                       #
	# initialize structs and timeroutput + reset inputdictionary to default #
	#                                                                       #
	#########################################################################

	potential = Potential()
	system    = System() #default setup but gets overridden later!
	output    = Output()
	files     = Files()

	files.to = TimerOutput()

	@timeit files.to "main" begin

		[
			inputDictionary[key] = "" for
			key in keys(inputDictionary)
		] #to reset dict to default values if calculation started in same repl session

		files.inputFileName = inputFileName

		###################
		#                 #
		# read input file #
		#                 #
		###################

		readInputFile(inputFileName)

		#########################
		#                       #
		# parse and check input #
		#                       #
		#########################

		checkInput(potential)
		checkInput(system)
		checkInput(files)
		checkInput(output)

		#################################
		#                               #
		# print program info to logfile #
		#                               #
		#################################

		init_logfile(files)

		#######################
		#                     #
		# read potential file #
		#                     #
		#######################

		readPotential(potential, files)

		################
		#              #
		# setup system #
		#              #
		################

		setupSystem(potential, system)

		################################################################################
		#                                                                              #
		# build ∇ and Δ matrix - caution ∇ matrix later modified according to k-points #
		#                                                                              #
		################################################################################

		buildΔ(system, potential)
		build∇(system, potential)

		##########################################
		#                                        #
		# print sparsity information to log file # TODO: modify this comment box
		#                                        #
		##########################################

		files.eigenvalueFileName = "eigenvalues.dat"
		files.bandStructureFileName = "bandstructure.dat"

		isfile(files.eigenvalueFileName) &&
			rm(files.eigenvalueFileName) #rm eigenvalue file if it exists TODO: think of a way to restart calculation for different k

		fileInfo(files, potential, system)

		systemInfo(files, potential, system)

		sparseInfo(files, system)

		########################################################################################
		#                                                                                      #
		# loop over all k-points (for non reciprocal system only one single point calculation) #
		#                                                                                      #
		########################################################################################

		@timeit files.to "loop" begin
			for (i, k) in enumerate(potential.kpoints)

				#########################################
				#                                       #
				# shift potential to pot_min equals 0.0 #
				#                                       #
				#########################################

				potential.potential =
					potential.potential .- potential.shift

				##############################
				#                            #
				# solve Schrödinger equation #
				#                            #
				##############################

				@timeit files.to "solve" solve(
					potential,
					system,
					output,
					k,
					files,
				)

				####################################
				#                                  #
				# calculate k from mass weighted k #
				#                                  #
				####################################

				k = k .* sqrt.(potential.mass)

				#############################################################
				#                                                           #
				# setup all file names depending on the momentanous k-point #
				#                                                           #
				#############################################################

				k_string = join(
					ustrip.(
						uconvert.(
							potential.coordsUnit^(-1),
							k ./ potential.internalElemCoords,
						)
					),
					"_",
				)

				if system.reciprocal
					files.eigenvectorFileName             = "eigenvectors_k_$(k_string).dat"
					files.eigenvectorShiftedFileName      = "eigenvectors_shifted_k_$(k_string).dat"
					files.imag_eigenvectorFileName        = "imag_eigenvectors_k_$(k_string).dat"
					files.imag_eigenvectorShiftedFileName = "imag_eigenvectors_shifted_k_$(k_string).dat"
					files.frequencyFileName               = "frequencies_k_$(k_string).dat"
				else
					files.eigenvectorFileName        = "eigenvectors.dat"
					files.eigenvectorShiftedFileName = "eigenvectors_shifted.dat"
					files.frequencyFileName          = "frequencies.dat"
				end

				########################################
				#                                      #
				# shift potential back to input values #
				#                                      #
				########################################

				potential.potential =
					potential.potential .+ potential.shift

				#######################################
				#                                     #
				# convert k-values back to input unit #
				#                                     #
				#######################################

				k =
					ustrip.(
						uconvert.(
							potential.coordsUnit^(-1),
							k ./ potential.internalElemCoords,
						)
					)

				############################################################
				#                                                          #
				# print eigenvalues, eigenvectors and frequencies to files #
				#                                                          #
				############################################################

				printEigenvalues(potential, output, files, k)
				printEigenvectors(
					potential,
					system,
					output,
					files,
					k,
				)
				printFrequencies(
					potential,
					system,
					output,
					files,
					k,
				)

				println(
					i,
					"/",
					length(potential.kpoints),
					" Done",
				)
			end
		end

		################################################################
		#                                                              #
		# if bandstructure is requested than write band structure file #
		#                                                              #
		################################################################

		potential.bandStructure && printBandStructure(
			potential,
			files,
			potential.kpoints,
		)

	end

	################################
	#                              #
	# print timings of calculation #
	#                              #
	################################

	files.timingsFile = open(files.timingsFileName, "w")
	show(files.to)
	println()
	show(files.timingsFile, files.to)

	############
	#          #
	# clean up #
	#          #
	############

	close(files.timingsFile)
	close(files.logFile)

end
