using DelimitedFiles

data = readdlm("eigenvalues.dat"; skipstart=1)

file = open("solution.dat", "w")

for i in 1:length(data[:,1])
    print(file, data[i,1], " ", data[i,2])
    for j in 1:length(data[1,:])-2
        print(file, " ", data[i,2+j])
    end
    println(file)
end
for i in 1:length(data[:,1])
    print(file, -data[i,1], " ", data[i,2])
    for j in 1:length(data[1,:])-2
        print(file, " ", data[i,2+j])
    end
    println(file)
end
for i in 1:length(data[:,1])
    print(file, data[i,1], " ", -data[i,2])
    for j in 1:length(data[1,:])-2
        print(file, " ", data[i,2+j])
    end
    println(file)
end
for i in 1:length(data[:,1])
    print(file, -data[i,1], " ", -data[i,2])
    for j in 1:length(data[1,:])-2
        print(file, " ", data[i,2+j])
    end
    println(file)
end

close(file)
