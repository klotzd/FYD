#!/usr/bin/env bash

echo "starting solution sequence"

julia fluidization.jl 
julia masstransfer.jl 
julia -i SOLVE.jl
