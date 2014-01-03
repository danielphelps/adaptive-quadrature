(ns adaptive-quadrature.core
  "Implements the adaptive quadrature described on page 511 of Numerical 
   Analysis Kinkade et al."
  (:require [clojure.math.numeric-tower :as math]))

(defn- simpsons-estimate
  "Equation '8.5' page 509 - approximates the integral of f over [a b]."
  [f a b h]
  (* (/ h 3.) (+ (f a) (* 4 (f (+ a h))) (f b))))

(defn- close-enough?
  "Finds if |a - b| < |error|."
  [a b error]
    (< (math/abs (- a b)) (math/abs error)))

(defn- insured-approximation
  "Equation 7 page 509 in Kinkade et al."
  [S* S** S]
  (+ S* S** (* (/ 1. 15) (+ S* S** (* -1 S)))))

(defn- adapt-quad-internal
  "Do not call this fn directly.  Start with adaptive-quadrature instead."
  [f delta eps n k sigma a h fa fc fb S]
  (let [b (+ a (* 2. h)) c (+ a h) h (/ h 2.)
	S-left  (simpsons-estimate f a c h)
	S-right (simpsons-estimate f c b h)]
    (cond 
      (close-enough? (+ S-left S-right) S (/ (* 60. eps h) delta))
                          (+ sigma (insured-approximation S-left S-right S))
      (>= k n) (throw (Exception. (str "Failure:  k >= n.  sigma = " sigma)))
      :else (+ (adapt-quad-internal   ;From a to the midpoint
                 f delta eps n (inc k) sigma a h fa (f (+ a h)) fc S-left)
               (adapt-quad-internal              
                 f delta eps n (inc k) ;From the midpoint to b
                 sigma (+ a (* 2. h)) h fc (f (+ a (* 3. h))) fb S-right)))))
        
(defn adaptive-quadrature
  "Approximates the definite integral of f over [a b] with an error less 
  or equal than eps.  f is a real valued function of one real argument.  
  The parameter n specifies how many recursive calls are allowed.  An 
  exception is thrown before the n+1st recursive call."
  [f a b eps n]
  (let [delta (- b a)  sigma 0  
        h (/ delta 2.) c (/ (+ a b) 2.) 
        k 1            fa (f a)  
        fb (f b)       fc (f c)  
        S (simpsons-estimate f a b h)]
    (adapt-quad-internal f delta eps n k sigma a h fa fc fb S)))
