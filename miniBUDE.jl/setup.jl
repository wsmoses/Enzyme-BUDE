#!/bin/env julia
import Pkg
Pkg.activate("./Threaded")

Enzyme_path = Base.find_package("Enzyme")
if Enzyme_path === nothing
    error("First instantiate Project")
end
dir = joinpath(dirname(dirname(Enzyme_path)), "deps")

@info "Building Enzyme from inc"

run(`$(Base.julia_cmd()) --project=$(dir) -e 'import Pkg; Pkg.instantiate()'`)
run(`$(Base.julia_cmd()) --project=$(dir) $(dir)/build_local.jl --branch inc`)

cp(joinpath(dirname(dir), "LocalPreferences.toml"), "Threaded/LocalPreferences.toml", force=true)
