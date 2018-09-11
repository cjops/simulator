using DataFrames, CSV, Plots, Distributions
plotlyjs()
if (isdefined(:IJulia))
    IJulia.clear_output()
end

function get_genotypes(num_genotypes)
    #num_genotypes = 2^num_loci
    bits = floor(Int, log2(num_genotypes-1))+1
    genotypes = [bin(n, bits) for n = 0:num_genotypes-1]
    return genotypes
end

function get_dataset1()
    # import the three data tables included in the supplement
    dataset = []
    for i = 1:3
        df = CSV.read(
            "table$i.csv",
            types=Dict(1=>String),
            allowmissing=:none,
            rows=17)
        delete!(df, :Genotype)
        rename!(df, Symbol("No drug") => :NoDrug)
        push!(dataset, df)
    end
    return dataset
end

function get_dataset2()
    barlow = CSV.read(
        "barlow.csv",
        allowmissing=:none,
        transpose=true,
        types=Dict(1=>String))
    # the order of the genotypes needs to match our current code
    sort!(barlow, :Label)
    delete!(barlow, :Label)
    return barlow
end

function get_mutational_neighbors(genotypes::Array{String})
    all_neighbors = []
    for (i, genotype) in enumerate(genotypes)
        neighbors = []
        for j in eachindex(genotype)
            neighbor = collect(genotype)
            if neighbor[j] == '0'
                neighbor[j] = '1'
            else
                neighbor[j] = '0'
            end
            neighbor = findfirst(genotypes, join(neighbor))
            push!(neighbors, neighbor)
        end
        push!(all_neighbors, neighbors)
    end
    return all_neighbors
end

function growth!(r, N)
    for i in eachindex(N)
        N[i] = floor(N[i] * exp(r[i]))
        if rand() < (N[i] * exp(r[i]) - floor(N[i] * exp(r[i])))
            N[i] += 1
        end
    end
end

function mutation!(r, N, P_m, mutational_neighbors)
    oldN = copy(N)
    for i in eachindex(N)
        avg_num_mutants = oldN[i] * P_m
        # sample from a Poisson distribution
        num_mutants = rand(Poisson(avg_num_mutants))
        for m = 1:num_mutants
            N[mutational_neighbors[i, rand(1:end)]] += 1
        end
        N[i] -= num_mutants
    end
end

function death!(r, N, K)
    ΣN = sum(N)
    for i in eachindex(N)
        freq = N[i] / ΣN
        N[i] = floor(freq * K)

        if rand() < (freq * K - floor(freq * K))
            N[i] += 1
        end
    end
end

# given a fitness landscape and a seed genotype, simulate population growth
# for a certain number of timesteps (by default 1200). returns a tuple of
# critical times and a trace of all population sizes at all timesteps.
function run_simulation(
        landscape::Array{T,1},
        seed,
        population = [];
        t::Number = 1200, # timesteps
        P_m::Number = 1.0e-8, # probability of mutation
        K::Number = 1.0e+9 # carrying capacity
    ) where {T<:Number}
    ordering = get_genotypes(length(landscape))
    mutational_neighbors = get_mutational_neighbors(ordering)
    if typeof(seed) == String
        seed = findfirst(ordering, seed)
    end
    if isempty(population)
        population=zeros(Int64, length(landscape))
        population[seed] = K
    end
    optimal = indmax(landscape) # index of globally optimum genotype
    dominant = indmax(population) # index of currently dominant genotype
    T_1 = 0 # time of first appearance of optimal genotype
    T_d = 0 # time to dominance
    T_f = 0 # time to fixation
    trace = Array{Int64}(t, length(population))
    trace[1,:] = population # first timestep is the starting population
    trajectory = [dominant] # a potential trajectory
    if seed == optimal
        T_1 = T_d = T_f = 1
    end
    for t = 2:t
        # go through each stage of the simulation
        growth!(landscape, population)
        mutation!(landscape, population, P_m, mutational_neighbors)
        death!(landscape, population, K)

        trace[t,:] = population #copy final population to the trace

        # check for a change in the dominant genotype
        if indmax(population) != dominant
            dominant = indmax(population)
            push!(trajectory, dominant)
        end

        # check for critical times
        if T_1 == 0 && population[optimal] > 0
            T_1 = t
        end
        if T_d == 0 && population[optimal] > 0.5 * K
            T_d = t
        end
        if T_f == 0 && population[optimal] > 0.99 * K
            T_f = t
        end
    end
    return (T_1, T_d, T_f), trace, trajectory
end

# run simulation over multiple landscapes -- needs some work
function run_simulation(
        landscapes::Array{Tuple{Array{T,1},X},1},
        seed,
        population = [];
        t::Number = 1200, # total timesteps
        P_m::Number = 1.0e-8, # probability of mutation
        K::Number = 1.0e+9 # carrying capacity
    ) where {T<:Number, X<:Integer}
    total_timesteps = sum([landscapes[i][2] for i in eachindex(landscapes)])
    if t > total_timesteps
        push!(landscapes, (landscapes[end][1], t - total_timesteps))
    end
    criticaltimes, trace, trajectory = run_simulation(
        landscapes[1][1],
        seed,
        population,
        t=landscapes[1][2],
        P_m=P_m,
        K=K)
    tsofar = landscapes[1][2]
    for (landscape, timesteps) in landscapes[2:end]
        newcrit, newtrace, newtraj = run_simulation(
            landscape,
            seed,
            trace[end,:],
            t=timesteps,
            P_m=P_m,
            K=K)
        trace = vcat(trace[1:end-1,:], newtrace)
        append!(trajectory, newtraj)
        tempcrit = []
        for i in eachindex(criticaltimes)
            if criticaltimes[i] == 0 && newcrit[i] != 0 && indmax(landscape) == indmax(landscapes[1][1])
                push!(tempcrit, tsofar + newcrit[i])
            else
                push!(tempcrit, criticaltimes[i])
            end
        end
        criticaltimes = Tuple(tempcrit)
        tsofar += timesteps
    end
    criticaltimes, trace, trajectory
end

# given a trace, the critical times, and an array of genotypes (by string or
# index), output a nice plot
function plot_simulation(
        trace,
        criticaltimes,
        genotypestoplot=1:size(trace, 2);
        vlines=3
    )
    ordering = get_genotypes(size(trace, 2))

    # if given strings, convert to indices
    if typeof(genotypestoplot) == Array{String,1}
        genotypestoplot = [findfirst(ordering, x) for x in genotypestoplot]
    end

    # sort by first appearance and only display genotypes that appeared
    firstappearances = [findfirst(trace[:,x]) for x in genotypestoplot]
    p = sortperm(firstappearances)
    genotypestoplot = genotypestoplot[p]
    genotypestoplot = genotypestoplot[length(genotypestoplot)-length(find(firstappearances))+1:length(genotypestoplot)]

    data = [trace[:,n] for n in genotypestoplot]

    #plot() only accepts a row vector of strings, hence the reshape
    labels = reshape([ordering[n] for n in genotypestoplot], 1, :)

    # plot abundances
    plot(   data,
            label=labels,
            ylims=(10^0,10^10),
            yscale=:log10,
            xlims=(0,size(trace,1)),
            xlabel="Timestep",
            ylabel="Abundance of each genotype",
            legendtitle="Genotype")

    # add vertical lines for each of the critical times to our plot
    sub = ("1", "d", "f")
    labels = ("First appearance", "Time to dominance", "Time to fixation")
    vrange = 1:0
    if vlines == 3
        vrange = 1:3
    elseif vlines == 1
        vrange = 3:3
    end
    for n=vrange
        if criticaltimes[n] != 0
            ann =  [(   criticaltimes[n],10,
                        text("<i>T</i><sub>" * sub[n] * "</sub>",10,:left))]
            # seriestype=:vline suddenly stopped working with plotly, so we
            # are using Shapes
            plot!(  [Shape([criticaltimes[n], criticaltimes[n]], [1, 10^10])],
                    label=labels[n],
                    primary=false,
                    annotations=ann)
        end
    end
    return plot!()
end
