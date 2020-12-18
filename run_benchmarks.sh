#!/usr/bin/env bash

set -v

# C
(
    cd c
)
# Clojure
(
    cd clojure
)
# Common Lisp
(
    cd cl
)
# C++
(
    cd cpp
)
# Go
(
    cd go
)
# Hy
(
    cd hy
)
# Julia
(
    cd julia
    cd OverlapIntegrals
    julia --project=. test/bench.jl
)
# Kotlin
(
    cd kotlin
)
# Nim
(
    cd nim
    nim r -d:release --opt:speed bench.nim
)
# Python
(
    cd python
)
# Racket
(
    cd racket
)
# Ruby
(
    cd ruby
)
# Rust
(
    cd rust/overlapintegrals
)

set +v
