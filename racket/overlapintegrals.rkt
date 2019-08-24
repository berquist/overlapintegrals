#lang racket/base

(require math)

(define (fact2 n)
  (if (<= n 0)
      1
      (* n (fact2 (- n 2)))))

(define (binomial-prefactor s ia ib xpa xpb)
  (for/fold ([acc 0])
            ([t (in-range 0 (+ s 1))])
    (if (and (<= (- s ia) t) (<= t ib))
        (+ acc
           (*
            (binomial ia (- s t))
            (binomial ib t)
            (expt xpa (+ (- ia s) t))
            (expt xpb (- ib t))))
        acc)))

(define (overlap1d l1 l2 pax pbx gamma)
  (for/fold ([acc 0])
            ([i (in-range 0 (add1 (floor (* 0.5 (+ l1 l2)))))])
    (+ acc
       (/ (* (binomial-prefactor (* 2 i) l1 l2 pax pbx)
             (fact2 (sub1 (* 2 i))))
          (expt (* 2 gamma) i)))))

(define (dist2 xa ya za xb yb zb)
  (+
   (expt (- xa xb) 2)
   (expt (- ya yb) 2)
   (expt (- za zb) 2)))

(define (product-center-1d za xa zb xb)
  (/ (+ (* za xa) (* zb xb))
     (+ za zb)))

(define (tho66 alpha1 alpha2 ra rb la lb)
  (define gamma (+ alpha1 alpha2))
  (define xa (list-ref ra 0))
  (define ya (list-ref ra 1))
  (define za (list-ref ra 2))
  (define xb (list-ref rb 0))
  (define yb (list-ref rb 1))
  (define zb (list-ref rb 2))
  (define rab2 (dist2 xa ya za xb yb zb))
  (define l1 (list-ref la 0))
  (define m1 (list-ref la 1))
  (define n1 (list-ref la 2))
  (define l2 (list-ref lb 0))
  (define m2 (list-ref lb 1))
  (define n2 (list-ref lb 2))
  (define xp (product-center-1d alpha1 xa alpha2 xb))
  (define yp (product-center-1d alpha1 ya alpha2 yb))
  (define zp (product-center-1d alpha1 za alpha2 zb))
  (define pre (* (exp (/ (* (- alpha1) alpha2 rab2) gamma))
                 (expt (/ pi gamma) 1.5)))
  (define wx (overlap1d l1 l2 (- xp xa) (- xp xb) gamma))
  (define wy (overlap1d m1 m2 (- yp ya) (- yp yb) gamma))
  (define wz (overlap1d n1 n2 (- zp za) (- zp zb) gamma))
  (* pre wx wy wz))

(provide fact2
         binomial-prefactor
         overlap1d
         tho66)
