(ns adaptive-quadrature.test.core
  (:use [adaptive-quadrature.core])
  (:use [clojure.math.numeric-tower])
  (:use [clojure.test]))

(defn approximately [a b eps]
  (< (abs (- a b)) eps))

(defn header [desc a b n eps sigma actual]
  (print "problem " 
          desc 
	     (format "n = %s, eps = %s, sigma = %s\n actual = %s\n" 
                  n eps sigma actual)
	     "---------------------\n"))

(deftest problem-a
  (let [a 0 
        b 1 
        n 22 
        eps 1.0E-6 
        actual (/ 2 3.)
	    sigma (time (adaptive-quadrature sqrt a b eps n))]
    (header "a) Integrate x^1/2 over [0 1].\n" 
            a b n eps sigma actual)
    (is (approximately sigma actual eps) 
        "Problem a is outside tolerance.")))

(deftest problem-b
  (let [a 0 
        b 1 
        n 29 
        eps 1.0E-4 
        actual (/ 2 3.)
        f (fn [x] (sqrt (- 1 x)))
        sigma (time (adaptive-quadrature f a b eps n))]
    (header "b) Integrate (1 - x)^1/2 over [0 1].\n" 
             a b n eps sigma actual)
    (is (approximately sigma actual eps) 
        "Problem b is outside tolerance.")))

(deftest problem-c
  (let [a 0 
        b 1 
        n 47 
        eps 1.0E-6 
        actual (/ 4 5.)
	    f (fn [x] (expt (- 1 x) 0.25))
	    sigma (time (adaptive-quadrature f a b eps n))]
    (header "c) Integrate (1 - x)^1/4 over [0 1].\n" 
            a b n eps sigma actual)
    (is (approximately sigma actual eps) 
        "Problem c is outside tolerance.")))

(deftest problem-c2
  (let [a 0 
        b 1 
        n 57 
        eps 1.0E-16 
        actual (/ 4 5.)
	    sigma (time (adaptive-quadrature (fn [x] (expt (- 1 x) 0.25)) a b eps n))]
    (header "c2) Integrate (1 - x)^1/4 over [0 1].\n" 
             a b n eps sigma actual)
    (is (approximately sigma actual eps) 
        "Problem c2 is outside tolerance.")))
