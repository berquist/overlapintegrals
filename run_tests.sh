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
# C++
(
    cd c++
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
    nim c -r overlapintegrals.nim
)
# Python
(
    cd python
    python -m pytest -v
)
# Racket
(
    cd racket
    racket overlapintegrals-test.rkt
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
