function test_CLI()

    # --help and --version print and return success without exiting
    @test (@suppress Numerov.CLI.command_main(["--help"]))    == 0
    @test (@suppress Numerov.CLI.command_main(["--version"])) == 0

    # a full run through the CLI entry point produces the output files
    run_testcase("1DHarmonicOscillator") do
        @test (@suppress Numerov.CLI.command_main(["input.in"])) == 0
        @test isfile("eigenvalues.dat")
        @test isfile("Numerov.out")
    end

    # invalid input prints a single-line error and exits with status 1; the
    # error path calls exit(), so it must run in a subprocess
    mktempdir() do tmp
        cd(tmp) do
            write("bad.in", "bogus-keyword = 42\n")
            script = "using Numerov; exit(Numerov.CLI.command_main([\"bad.in\"]))"
            cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`
            err = IOBuffer()
            p = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = err))
            @test p.exitcode == 1
            errlines = filter(!isempty, split(String(take!(err)), '\n'))
            @test length(errlines) == 1
            @test startswith(errlines[1], "error: ArgumentError")
        end
    end
end
