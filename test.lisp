;;;; test.lisp — Tests for Conservation Spectral SDK
;;;; Tests numerical computation, symbolic algebra, and theorem proving

;;; Load all modules
(load "matrix.lisp")
(load "symbolic.lisp")
(load "conservation.lisp")
(load "prover.lisp")

;;; Use a dedicated test package to avoid name conflicts
(defpackage :conservation-test
  (:use :cl)
  (:shadowing-import-from :conservation-spectral
    :make-graph :graph-add-edge :graph-add-vertex :graph-n :graph-vertices :graph-edges
    :build-laplacian :build-normalized-laplacian :eigendecompose
    :conservation-ratio :spectral-gap :cheeger-constant :analyze :full-analysis
    :make-tracker :tracker-feed :tracker-check
    :fingerprint-compute :fingerprint-compare
    :conservation-report :conservation-report-spectral-gap :conservation-report-cheeger
    :conservation-report-eigenvalues :conservation-report-ratios
    :conservation-report-fingerprint :conservation-report-graph-info
    :symbolic-laplacian-from-edges :symbolic-conservation
    :vertex :make-vertex :vertex-id :vertex-attribute
    :edge :make-edge :edge-from :edge-to :edge-weight)
  (:shadowing-import-from :conservation-matrix
    :make-matrix :matrix-ref :matrix-set :matrix-rows :matrix-cols
    :make-vec :matrix-vector-multiply :vector-dot :vector-norm
    :matrix-copy :matrix-add :matrix-sub :matrix-scale :matrix-multiply
    :matrix-transpose :matrix-identity)
  (:shadowing-import-from :conservation-symbolic
    :sym+ :sym* :sym- :sym/ :sym-expt :sym-negate :sym-simplify
    :sym->string :sym-expand :sym-substitute :sym-diff :sym-eval
    :collect-variables)
  (:shadowing-import-from :conservation-prover
    :prove-conservation-bounded :prove-conservation-for-all-attributes
    :conservation-proof-result-p :conservation-proof-result-verified
    :conservation-proof-result-description :conservation-proof-result-evidence))

(in-package :conservation-test)

(defvar *tests-passed* 0)
(defvar *tests-failed* 0)

(defmacro test (name &body body)
  `(progn
     (format t "~%TEST: ~a ... " ,name)
     (handler-case
         (progn ,@body
                (incf *tests-passed*)
                (format t "OK"))
       (error (c)
         (incf *tests-failed*)
         (format t "FAILED: ~a" c)))))

(defmacro assert-approx (a b &optional (tol 0.01))
  `(unless (< (abs (- ,a ,b)) ,tol)
     (error "Expected ~f ≈ ~f (tol=~f)" ,a ,b ,tol)))

;;; ============================================================
;;; Matrix tests
;;; ============================================================

(test "Matrix creation and access"
  (let ((m (make-matrix 3 3)))
    (matrix-set m 0 0 5.0d0)
    (assert (= (matrix-ref m 0 0) 5.0d0))
    (assert (= (matrix-rows m) 3))
    (assert (= (matrix-cols m) 3))))

(test "Matrix addition"
  (let ((a (make-matrix 2 2))
        (b (make-matrix 2 2)))
    (matrix-set a 0 0 1.0d0) (matrix-set a 0 1 2.0d0)
    (matrix-set a 1 0 3.0d0) (matrix-set a 1 1 4.0d0)
    (matrix-set b 0 0 5.0d0) (matrix-set b 0 1 6.0d0)
    (matrix-set b 1 0 7.0d0) (matrix-set b 1 1 8.0d0)
    (let ((c (matrix-add a b)))
      (assert (= (matrix-ref c 0 0) 6.0d0))
      (assert (= (matrix-ref c 1 1) 12.0d0)))))

(test "Matrix multiply"
  (let ((a (make-matrix 2 3))
        (b (make-matrix 3 2)))
    (matrix-set a 0 0 1.0d0) (matrix-set a 0 1 2.0d0) (matrix-set a 0 2 3.0d0)
    (matrix-set a 1 0 4.0d0) (matrix-set a 1 1 5.0d0) (matrix-set a 1 2 6.0d0)
    (matrix-set b 0 0 7.0d0) (matrix-set b 0 1 8.0d0)
    (matrix-set b 1 0 9.0d0) (matrix-set b 1 1 10.0d0)
    (matrix-set b 2 0 11.0d0) (matrix-set b 2 1 12.0d0)
    (let ((c (matrix-multiply a b)))
      ;; [1 2 3] * [7  8 ]   = [58  64 ]
      ;; [4 5 6]   [9  10]     [139 154]
      ;;           [11 12]
      (assert (= (matrix-ref c 0 0) 58.0d0))
      (assert (= (matrix-ref c 1 1) 154.0d0)))))

(test "Matrix-vector multiply"
  (let ((m (make-matrix 2 2))
        (v (make-vec 2)))
    (matrix-set m 0 0 1.0d0) (matrix-set m 0 1 2.0d0)
    (matrix-set m 1 0 3.0d0) (matrix-set m 1 1 4.0d0)
    (setf (aref v 0) 1.0d0 (aref v 1) 2.0d0)
    (let ((r (matrix-vector-multiply m v)))
      (assert (= (aref r 0) 5.0d0))
      (assert (= (aref r 1) 11.0d0)))))

(test "Identity matrix"
  (let ((I (matrix-identity 3)))
    (dotimes (i 3)
      (dotimes (j 3)
        (assert (= (matrix-ref I i j) (if (= i j) 1.0d0 0.0d0)))))))

(test "Matrix transpose"
  (let ((m (make-matrix 2 3)))
    (matrix-set m 0 0 1.0d0) (matrix-set m 0 1 2.0d0) (matrix-set m 0 2 3.0d0)
    (matrix-set m 1 0 4.0d0) (matrix-set m 1 1 5.0d0) (matrix-set m 1 2 6.0d0)
    (let ((t_ (matrix-transpose m)))
      (assert (= (matrix-rows t_) 3))
      (assert (= (matrix-ref t_ 0 1) 4.0d0))
      (assert (= (matrix-ref t_ 2 1) 6.0d0)))))

;;; ============================================================
;;; Symbolic algebra tests
;;; ============================================================

(test "Symbolic simplification"
  (let ((expr (sym+ 'x 0)))
    (assert (equal (sym-simplify expr) 'x)))
  (let ((expr (sym* 'x 1)))
    (assert (equal (sym-simplify expr) 'x)))
  (let ((expr (sym* 'x 0)))
    (assert (equal (sym-simplify expr) 0)))
  (let ((expr (sym+ 3 4)))
    (assert (equal (sym-simplify expr) 7))))

(test "Symbolic differentiation"
  ;; d/dx(x²) = 2x
  (let* ((expr (sym-expt 'x 2))
         (deriv (sym-diff expr 'x)))
    (format t " d/dx(x²) = ~a" (sym->string (sym-simplify deriv))))
  ;; d/dx(x*y + 3) = y
  (let* ((expr (sym+ (sym* 'x 'y) 3))
         (deriv (sym-diff expr 'x)))
    (format t " d/dx(xy+3) = ~a" (sym->string (sym-simplify deriv)))))

(test "Symbolic substitution"
  (let ((expr '(* x y))
        (result (sym-substitute '(* x y) 'x 5)))
    (assert (equal result '(* 5 y)))))

(test "Symbolic evaluation"
  (let ((result (sym-eval '(+ (* 2 x) 3) '((x . 5)))))
    (assert (= result 13))))

(test "Symbolic pretty printing"
  (let ((expr (sym+ (sym* 'x 'y) 3)))
    (format t " ~a" (sym->string expr))))

;;; ============================================================
;;; Graph + Laplacian tests
;;; ============================================================

(test "Build graph Laplacian (triangle)"
  ;; Triangle graph: 3 vertices, edges (0-1), (1-2), (0-2) all weight 1
  (let ((g (make-graph 3)))
    (graph-add-edge g 0 1 :weight 1.0d0)
    (graph-add-edge g 1 2 :weight 1.0d0)
    (graph-add-edge g 0 2 :weight 1.0d0)
    (let ((L (build-laplacian g)))
      ;; Diagonal should be 2 (degree 2 for each vertex)
      (assert-approx (matrix-ref L 0 0) 2.0d0)
      (assert-approx (matrix-ref L 1 1) 2.0d0)
      (assert-approx (matrix-ref L 2 2) 2.0d0)
      ;; Off-diagonal should be -1
      (assert-approx (matrix-ref L 0 1) -1.0d0)
      (assert-approx (matrix-ref L 1 0) -1.0d0)
      (assert-approx (matrix-ref L 0 2) -1.0d0)
      (assert-approx (matrix-ref L 2 0) -1.0d0))))

(test "Build Laplacian (path graph)"
  ;; Path: 0-1-2-3
  (let ((g (make-graph 4)))
    (graph-add-edge g 0 1 :weight 1.0d0)
    (graph-add-edge g 1 2 :weight 1.0d0)
    (graph-add-edge g 2 3 :weight 1.0d0)
    (let ((L (build-laplacian g)))
      ;; Degrees: 1, 2, 2, 1
      (assert-approx (matrix-ref L 0 0) 1.0d0)
      (assert-approx (matrix-ref L 1 1) 2.0d0)
      (assert-approx (matrix-ref L 2 2) 2.0d0)
      (assert-approx (matrix-ref L 3 3) 1.0d0)
      ;; Row sums should be 0
      (dotimes (i 4)
        (let ((row-sum 0.0d0))
          (dotimes (j 4)
            (incf row-sum (matrix-ref L i j)))
          (assert-approx row-sum 0.0d0 1e-10))))))

(test "Laplacian is positive semi-definite (eigenvalues ≥ 0)"
  (let ((g (make-graph 5)))
    (graph-add-edge g 0 1) (graph-add-edge g 1 2)
    (graph-add-edge g 2 3) (graph-add-edge g 3 4)
    (graph-add-edge g 0 4)
    (let ((L (build-laplacian g)))
      (multiple-value-bind (eigenvalues eigenvectors)
          (eigendecompose L)
        (declare (ignore eigenvectors))
        ;; All eigenvalues should be ≥ 0
        (dotimes (i (length eigenvalues))
          (assert (>= (aref eigenvalues i) -0.1))) ;; small tolerance
        ;; Smallest eigenvalue should be ~0
        (assert-approx (aref eigenvalues 0) 0.0d0 0.1)))))

;;; ============================================================
;;; Conservation analysis tests
;;; ============================================================

(test "Conservation ratio computation"
  (let ((g (make-graph 4)))
    (graph-add-vertex g 0 :attribute 1.0d0)
    (graph-add-vertex g 1 :attribute 2.0d0)
    (graph-add-vertex g 2 :attribute 3.0d0)
    (graph-add-vertex g 3 :attribute 4.0d0)
    (graph-add-edge g 0 1) (graph-add-edge g 1 2) (graph-add-edge g 2 3)
    (let ((L (build-laplacian g)))
      (multiple-value-bind (eigenvalues eigenvectors)
          (eigendecompose L)
        (declare (ignore eigenvalues))
        ;; Conservation ratio for each eigenvector should be a number
        (dotimes (k 4)
          (let* ((attrs (make-vec 4))
                 (ratio (progn
                          (dotimes (i 4) (setf (aref attrs i)
                                                (vertex-attribute (aref (graph-vertices g) i))))
                          (conservation-ratio (nth k eigenvectors) attrs))))
            (assert (numberp ratio))
            (assert (>= ratio -0.001))))))))

(test "Spectral gap"
  (let ((g (make-graph 4)))
    (graph-add-edge g 0 1) (graph-add-edge g 1 2) (graph-add-edge g 2 3)
    (let ((L (build-laplacian g)))
      (multiple-value-bind (eigenvalues evs)
          (eigendecompose L)
        (declare (ignore evs))
        (let ((gap (spectral-gap eigenvalues)))
          (assert (> gap 0)))))))

(test "Cheeger constant"
  (let ((g (make-graph 6)))
    (graph-add-edge g 0 1) (graph-add-edge g 1 2) (graph-add-edge g 2 3)
    (graph-add-edge g 3 4) (graph-add-edge g 4 5)
    (let ((L (build-laplacian g)))
      (multiple-value-bind (eigenvalues eigenvectors)
          (eigendecompose L)
        (declare (ignore eigenvalues))
        (let ((h (cheeger-constant L (second eigenvectors))))
          (assert (> h 0))
          (assert (<= h 1.0)))))))

(test "Full analysis"
  (let ((g (make-graph 8)))
    ;; Musical chord progression: C-Am-F-G (I-vi-IV-V)
    (graph-add-vertex g 0 :attribute 1.0d0)   ; C
    (graph-add-vertex g 1 :attribute 0.8d0)    ; Am
    (graph-add-vertex g 2 :attribute 0.6d0)    ; F
    (graph-add-vertex g 3 :attribute 0.7d0)    ; G
    (graph-add-vertex g 4 :attribute 0.5d0)    ; Em
    (graph-add-vertex g 5 :attribute 0.4d0)    ; Dm
    (graph-add-vertex g 6 :attribute 0.3d0)    ; Bdim
    (graph-add-vertex g 7 :attribute 0.9d0)    ; C (resolution)
    ;; Harmonic connections
    (graph-add-edge g 0 1 :weight 0.8d0)  ; C→Am
    (graph-add-edge g 1 2 :weight 0.7d0)  ; Am→F
    (graph-add-edge g 2 3 :weight 0.9d0)  ; F→G
    (graph-add-edge g 3 0 :weight 0.6d0)  ; G→C
    (graph-add-edge g 0 4 :weight 0.5d0)  ; C→Em
    (graph-add-edge g 4 5 :weight 0.6d0)  ; Em→Dm
    (graph-add-edge g 5 2 :weight 0.7d0)  ; Dm→F
    (graph-add-edge g 3 6 :weight 0.3d0)  ; G→Bdim
    (graph-add-edge g 6 7 :weight 0.4d0)  ; Bdim→C
    (graph-add-edge g 7 0 :weight 1.0d0)  ; resolution
    (let ((report (analyze g)))
      (assert (conservation-report-p report))
      (assert (> (conservation-report-spectral-gap report) 0))
      (assert (>= (conservation-report-cheeger report) 0))
      (assert (> (length (conservation-report-fingerprint report)) 0)))))

;;; ============================================================
;;; Tracker tests
;;; ============================================================

(test "Tracker nominal"
  (let ((t_ (make-tracker 5)))
    ;; Fill baseline
    (dotimes (i 5)
      (tracker-feed t_ 1.0d0))
    ;; Nominal observation
    (assert (eq (tracker-feed t_ 1.1d0) :nominal))))

(test "Tracker anomaly detection"
  (let ((t_ (make-tracker 5)))
    ;; Fill baseline with consistent values
    (dotimes (i 5)
      (tracker-feed t_ 1.0d0))
    ;; Extreme observation should trigger
    (let ((result (tracker-feed t_ 100.0d0)))
      (assert (member result '(:warning :critical))))))

;;; ============================================================
;;; Fingerprint tests
;;; ============================================================

(test "Fingerprint same graph"
  (let* ((g (make-graph 4)))
    (graph-add-edge g 0 1) (graph-add-edge g 1 2) (graph-add-edge g 2 3)
    (let ((r1 (full-analysis g))
          (r2 (full-analysis g)))
      (let ((sim (fingerprint-compare
                  (conservation-report-fingerprint r1)
                  (conservation-report-fingerprint r2))))
        (assert-approx sim 1.0d0 0.01)))))

(test "Fingerprint different graphs"
  (let* ((g1 (make-graph 4))
         (g2 (make-graph 4)))
    (graph-add-edge g1 0 1) (graph-add-edge g1 1 2) (graph-add-edge g1 2 3)
    (graph-add-edge g2 0 1) (graph-add-edge g2 1 2)
    (graph-add-edge g2 2 3) (graph-add-edge g2 0 3) (graph-add-edge g2 1 3)
    (let ((r1 (full-analysis g1))
          (r2 (full-analysis g2)))
      (let ((sim (fingerprint-compare
                  (conservation-report-fingerprint r1)
                  (conservation-report-fingerprint r2))))
        (assert (< sim 1.0d0))))))

;;; ============================================================
;;; Symbolic conservation tests
;;; ============================================================

(test "Symbolic Laplacian"
  (let* ((sym-L (symbolic-laplacian-from-edges
                 '((0 1 w01) (1 2 w12) (0 2 w02)) 3)))
    ;; Check diagonal entries
    (format t " L[0][0] = ~a" (sym->string (sym-simplify (aref sym-L 0 0))))
    (format t " L[0][1] = ~a" (sym->string (sym-simplify (aref sym-L 0 1))))
    ;; Diagonal should have both weights
    (assert (not (zerop (aref sym-L 0 0))))))

(test "Symbolic conservation ratio"
  (let* ((sym-L (symbolic-laplacian-from-edges
                 '((0 1 w01) (1 2 w12)) 3))
         (sym-attrs (make-array 3 :initial-contents '(a0 a1 a2)))
         (sym-ratio (symbolic-conservation sym-L sym-attrs)))
    (format t "~%  Symbolic conservation: ~a" (sym->string (sym-simplify sym-ratio)))))

;;; ============================================================
;;; Prover tests
;;; ============================================================

(test "Prove conservation bounded"
  (let* ((g (make-graph 5))
         (report (progn
                   (graph-add-edge g 0 1) (graph-add-edge g 1 2)
                   (graph-add-edge g 2 3) (graph-add-edge g 3 4)
                   (graph-add-edge g 0 4)
                   (prove-conservation-bounded g))))
    (assert (conservation-proof-result-p report))
    (format t "~%  ~a" (conservation-proof-result-description report))))

(test "Prove conservation for all attributes"
  (let* ((g (make-graph 4))
         (report (progn
                   (graph-add-edge g 0 1) (graph-add-edge g 1 2)
                   (graph-add-edge g 2 3)
                   (prove-conservation-for-all-attributes g))))
    (assert (conservation-proof-result-verified report))
    (format t "~%  ~a" (conservation-proof-result-description report))))

;;; Summary
;;; ============================================================

(format t "~2%=====================================~%")
(format t "RESULTS: ~d passed, ~d failed, ~d total~%"
        *tests-passed* *tests-failed*
        (+ *tests-passed* *tests-failed*))
(format t "=====================================~%")

;; Exit with status
(sb-ext:exit :code (if (= *tests-failed* 0) 0 1))
