using BenchmarkTools, GLMakie

# quick and dirty julia implementation
function ising!(lat::Matrix{Int8}, t::Float64, nsteps::Int)
    β = 1 / t
    array = lat
    m, n = size(array)
    prob = [(exp(-2 * β * k) for k in -4:4)...,]
    @inbounds for _ in 1:nsteps, i in 1:n
        @simd for j in 1:m
            site = array[j, i]
            NN = array[j == m ? 1 : j + 1, i]
            SS = array[j == 1 ? m : j - 1, i]
            EE = array[j, i == n ? 1 : i + 1]
            WW = array[j, i == 1 ? n : i - 1]
            energy = site * (NN + SS + EE + WW)
            array[j, i] = rand() < prob[energy+5] ? -array[j, i] : array[j, i]
        end
    end
    return array
end

function plot_ising(lat::Matrix{Int8})
    array = lat
    heatmap(array, ticks=false)
end

T = 2 / log(1 + sqrt(2))
side = 1000
steps = 1000
array = rand(Int8[-1,1], side, side)

@elapsed ising!(array, T, steps) 