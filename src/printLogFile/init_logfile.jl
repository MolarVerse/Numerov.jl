function init_logfile(files::Files)

    user     = ENV["USER"]
    hostname = gethostname()
    nthreads = Sys.CPU_THREADS
    ram      = Sys.total_memory() / 10.0^9
    kernel   = Sys.KERNEL
    machine  = Sys.MACHINE
    
    stringbuffer = "\n" *
                   "\n" *
    "   #####################################################################################" * "\n" *
    "   #                                                                                   #" * "\n" *
    "   #   #      #                                                                        #" * "\n" *
    "   #   ##     #                                                                        #" * "\n" *
    "   #   # #    #                                                                        #" * "\n" *
    "   #   #  #   #  #    #  #########  #####  # #### ######  #       #      #####  #      #" * "\n" *
    "   #   #   #  #  #    #  #   #   #  #   #  ##     #    #   #     #           #  #      #" * "\n" *
    "   #   #    # #  #    #  #   #   #  #####  #      #    #    #   #            #  #      #" * "\n" *
    "   #   #     ##  #    #  #   #   #  #      #      #    #     # #     ##      #  #   #  #" * "\n" *
    "   #   #      #  ######  #   #   #  #####  #      ######      #      ##  #####  #####  #" * "\n" *
    "   #                                                                                   #" * "\n" *
    "   #####################################################################################" * "\n" *
    "                                                                                        " * "\n" *
    "                                                                                        " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "       | USER AND SYSTEM INFORMATION                                                |   " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "                                                                                        " * "\n" *
    "         user             : $user                                                       " * "\n" *
    "         hostname         : $hostname                                                   " * "\n" *
    "         number of threads: $nthreads                                                   " * "\n" *
    "         total RAM        : $(ram)GB                                                    " * "\n" *
    "         kernel           : $kernel                                                     " * "\n" *
    "         machine          : $machine                                                    " * "\n" *
    "         julia version    : $VERSION                                                    " * "\n" *
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end