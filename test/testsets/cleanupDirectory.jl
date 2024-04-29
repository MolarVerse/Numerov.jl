function cleanup_directory(input::String, potential::String)
    input_files = [input, potential]
    rm.(filter(x -> x ∉ input_files, readdir()))
end