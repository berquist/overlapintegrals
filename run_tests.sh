#!/usr/bin/env bash

set -v

# C
(
    cd gcc
)
# Clojure
(
    cd clojure
)
# C++
(
    cd g++
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
    julia test_overlapintegrals.jl
)
# Kotlin
(
    cd kotlin
)
# Nim
(
    cd nim
)
# Python
(
    cd python
    python -m pytest -v
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
    cd rust
    cargo test
)

set +v
