#lang racket/base

(require rackunit "overlapintegrals.rkt")

(test-case
    "fact2"
  (check-equal? (fact2 0) 1)
  (check-equal? (fact2 1) 1)
  (check-equal? (fact2 2) 2)
  (check-equal? (fact2 3) (* 3 1))
  (check-equal? (fact2 4) (* 4 2))
  (check-equal? (fact2 5) (* 5 3 1))
  (check-equal? (fact2 6) (* 6 4 2))
  (check-equal? (fact2 7) (* 7 5 3 1)))

(test-case
    "binomial-prefactor"
  (check-= (binomial-prefactor 1 1 1 0.1 0.2) 0.3 1e-7)
  (check-= (binomial-prefactor 1 1 1 0.3 0.4) 0.7 1e-7)
  (check-= (binomial-prefactor 1 3 1 0.1 0.2) 0.007 1e-7)
  (check-= (binomial-prefactor 2 3 1 0.1 0.2) 0.09 1e-7))

(test-case
    "overlap1d"
  (check-= (overlap1d 1 1 0.1 0.2 1.0) 0.52 1e-5)
  (check-= (overlap1d 3 1 0.1 0.2 1.0) 0.7952 1e-5))

(test-case
    "tho66"
  (define za 1.8)
  (define zb 2.8)
  (define ra '(0.0 0.0 0.0))
  (define rb '(0.5 0.8 -0.2))
  (check-= (tho66 za zb ra rb '(0 0 0) '(0 0 0)) 0.20373275913014607 1e-16)
  (check-= (tho66 za zb ra rb '(1 0 0) '(0 0 0)) 0.062005622343957505 1e-16)
  (check-= (tho66 za zb ra rb '(1 1 0) '(1 1 0)) -0.00043801221837779696 1e-16)
  (check-= (tho66 za zb ra rb '(2 1 0) '(1 1 0)) -0.0002385994651113168 1e-16))
