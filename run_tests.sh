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
    clisp overlapintegrals.lisp
    sbcl --script overlapintegrals.lisp
)
# C++
(
    cd cpp
    mkdir -p build
    cd build
    cmake ..
    make
    ./main.x
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
    julia --project -e "using Pkg; Pkg.test()"
)
# Kotlin
(
    cd kotlin
)
# Nim
(
    cd nim
    nim r -d:release overlapintegrals.nim
)
# Python
(
    cd python
    python -m pytest -v
)
# Racket
(
    cd racket
    raco test .
)
# Ruby
(
    cd ruby
)
# Rust
(
    cd rust/overlapintegrals
    cargo test
)

set +v
