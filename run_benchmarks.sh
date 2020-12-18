#!/usr/bin/env bash

set -v

# Julia
(
    cd julia
    cd OverlapIntegrals
    julia --project=. test/bench.jl
)
# Nim
(
    cd nim
    nim r -d:release --opt:speed bench.nim
)

set +v
