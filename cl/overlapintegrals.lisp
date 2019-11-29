(defun fact2 (n)
  (if (<= n 0)
      1
      (* n (fact2 (- n 2)))))

(defun dist2 (xa ya za xb yb zb)
  (+
   (expt (- xa xb) 2)
   (expt (- ya yb) 2)
   (expt (- za zb) 2)))

(defun fact (n)
  (if (<= n 0)
      1
      (* n (fact (1- n)))))

(defun binomial (n k)
  "Returns the binomial coefficient."
  (/ (fact n) (* (fact k) (fact (- n k)))))

(defun binomial-prefactor (s ia ib xpa xpb)
  (loop
    for v from 0 to s
    if (and (<= (- s ia) v) (<= v ib))
      summing
      (*
       (binomial ia (- s v))
       (binomial ib v)
       (expt xpa (+ (- ia s) v))
       (expt xpb (- ib v)))))

(defun overlap1d (l1 l2 pax pbx gamma)
  (loop
    for i from 0 to (floor (* 0.5 (+ l1 l2)))
    summing
    (/ (* (binomial-prefactor (* 2 i) l1 l2 pax pbx)
          (fact2 (1- (* 2 i))))
       (expt (* 2 gamma) i))))

(defun product-center-1d (za xa zb xb)
  (/ (+ (* za xa) (* zb xb))
     (+ za zb)))

(defun tho66 (alpha1 alpha2 ra rb la lb)
  (let* ((gamma (+ alpha1 alpha2))
         (xa (nth 0 ra))
         (ya (nth 1 ra))
         (za (nth 2 ra))
         (xb (nth 0 rb))
         (yb (nth 1 rb))
         (zb (nth 2 rb))
         (rab2 (dist2 xa ya za xb yb zb))
         (l1 (nth 0 la))
         (m1 (nth 1 la))
         (n1 (nth 2 la))
         (l2 (nth 0 lb))
         (m2 (nth 1 lb))
         (n2 (nth 2 lb))
         (xp (product-center-1d alpha1 xa alpha2 xb))
         (yp (product-center-1d alpha1 ya alpha2 yb))
         (zp (product-center-1d alpha1 za alpha2 zb))
         (pre (* (exp (/ (* (- alpha1) alpha2 rab2) gamma))
                 (expt (/ pi gamma) 1.5)))
         (wx (overlap1d l1 l2 (- xp xa) (- xp xb) gamma))
         (wy (overlap1d m1 m2 (- yp ya) (- yp yb) gamma))
         (wz (overlap1d n1 n2 (- zp za) (- zp zb) gamma)))
    (* pre wx wy wz)))

(assert (= (fact2 0) 1))
(assert (= (fact2 1) 1))
(assert (= (fact2 2) 2))
(assert (= (fact2 3) (* 3 1)))
(assert (= (fact2 4) (* 4 2)))
(assert (= (fact2 5) (* 5 3 1)))
(assert (= (fact2 6) (* 6 4 2)))
(assert (= (fact2 7) (* 7 5 3 1)))

(defun approx (actual expected thresh)
  (<= (abs (- actual expected)) thresh))

(assert (approx (binomial-prefactor 1 1 1 0.1 0.2) 0.3 1e-7))
(assert (approx (binomial-prefactor 1 1 1 0.3 0.4) 0.7 1e-7))
(assert (approx (binomial-prefactor 1 3 1 0.1 0.2) 0.007 1e-7))
(assert (approx (binomial-prefactor 2 3 1 0.1 0.2) 0.09 1e-7))

(assert (approx (overlap1d 1 1 0.1 0.2 1.0) 0.52 1e-5))
(assert (approx (overlap1d 3 1 0.1 0.2 1.0) 0.7952 1e-5))

(defparameter za 1.8)
(defparameter zb 2.8)
(defparameter ra '(0.0 0.0 0.0))
(defparameter rb '(0.5 0.8 -0.2))
(assert (approx (tho66 za zb ra rb '(0 0 0) '(0 0 0)) 0.20373275913014607 1e-7))
;; 16 for CLISP, 8 for SBCL ¯\_(ツ)_/¯
(assert (approx (tho66 za zb ra rb '(1 0 0) '(0 0 0)) 0.062005622343957505 1e-8))
(assert (approx (tho66 za zb ra rb '(1 1 0) '(1 1 0)) -0.00043801221837779696 1e-10))
(assert (approx (tho66 za zb ra rb '(2 1 0) '(1 1 0)) -0.0002385994651113168 1e-10))
