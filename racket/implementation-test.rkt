#lang racket/base

(require rackunit math)

(test-case
    "binomial"
  (check-eq? (binomial 0 0) 1)
  (check-eq? (binomial 1 0) 1)
  (check-eq? (binomial 2 0) 1)
  (check-eq? (binomial 3 0) 1)
  (check-eq? (binomial 0 1) 0)
  (check-eq? (binomial 0 2) 0)
  (check-eq? (binomial 0 3) 0)
  (check-eq? (binomial 4 2) 6)
  (check-eq? (binomial 10 3) 120))
