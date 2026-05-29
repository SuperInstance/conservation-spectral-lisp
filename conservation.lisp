;;;; conservation.lisp — Core Conservation Spectral Analysis
;;;; Numerical + Symbolic computation of conservation ratios
;;;; The ONLY implementation that can do ALGEBRA on conservation

(defpackage :conservation-spectral
  (:use :cl)
  (:export :build-laplacian :build-normalized-laplacian
           :eigendecompose :power-iteration
           :conservation-ratio :spectral-gap :cheeger-constant
           :analyze :full-analysis
           :symbolic-laplacian :symbolic-conservation
           :symbolic-laplacian-from-edges
           :graph :make-graph :graph-add-edge :graph-add-vertex
           :graph-n :graph-vertices :graph-edges
           :vertex :make-vertex :vertex-id :vertex-attribute
           :edge :make-edge :edge-from :edge-to :edge-weight
           :tracker :make-tracker :tracker-feed :tracker-check
           :fingerprint-compute :fingerprint-compare
           :conservation-report :conservation-report-spectral-gap :conservation-report-cheeger
           :report-eigenvalues :report-ratios :report-fingerprint
           :report-graph-info))

(in-package :conservation-spectral)

;; Import from matrix and symbolic packages
(use-package :conservation-matrix)
(use-package :conservation-symbolic)

;;; ============================================================
;;; Graph representation
;;; ============================================================

(defstruct vertex
  (id 0 :type fixnum)
  (attribute 0.0d0 :type double-float))

(defstruct edge
  (from 0 :type fixnum)
  (to 0 :type fixnum)
  (weight 1.0d0 :type double-float))

(defstruct (graph (:constructor %make-graph))
  (vertices nil :type (simple-array vertex (*)))
  (edges nil :type list)
  (n 0 :type fixnum))

(defun make-graph (n)
  "Create graph with n vertices"
  (let ((vs (make-array n :element-type t :initial-element nil)))
    (dotimes (i n)
      (setf (aref vs i) (make-vertex :id i :attribute 0.0d0)))
    (%make-graph :vertices vs :edges nil :n n)))

(defun graph-add-vertex (g id &key (attribute 0.0d0))
  "Set vertex attribute"
  (assert (< id (graph-n g)) (id) "Vertex id out of bounds")
  (setf (vertex-attribute (aref (graph-vertices g) id))
        (coerce attribute 'double-float))
  g)

(defun graph-add-edge (g from to &key (weight 1.0d0))
  "Add undirected edge"
  (push (make-edge :from from :to to :weight (coerce weight 'double-float))
        (graph-edges g))
  g)

;;; ============================================================
;;; Numerical Laplacian
;;; ============================================================

(defun build-laplacian (g)
  "Build graph Laplacian L = D - W from graph"
  (let* ((n (graph-n g))
         (lap (make-matrix n n)))
    (dolist (e (graph-edges g) lap)
      (let ((u (edge-from e))
            (v (edge-to e))
            (w (edge-weight e)))
        (matrix-set lap u v (- (matrix-ref lap u v) w))
        (matrix-set lap v u (- (matrix-ref lap v u) w))
        (matrix-set lap u u (+ (matrix-ref lap u u) w))
        (matrix-set lap v v (+ (matrix-ref lap v v) w))))))

(defun build-normalized-laplacian (g)
  "Build symmetric normalized Laplacian L = I - D^{-1/2} W D^{-1/2}"
  (let* ((n (graph-n g))
         (lap (build-laplacian g))
         (deg (conservation-matrix:make-vec n)))
    (dotimes (i n)
      (setf (aref deg i) (matrix-ref lap i i)))
    (let ((nlap (make-matrix n n)))
      (dotimes (i n nlap)
        (dotimes (j n)
          (let ((di (aref deg i))
                (dj (aref deg j)))
            (if (or (< di 1e-15) (< dj 1e-15))
                (matrix-set nlap i j (if (= i j) 1.0d0 0.0d0))
                (let ((val (/ (matrix-ref lap i j)
                              (* (sqrt di) (sqrt dj)))))
                  (matrix-set nlap i j val)))))))))

;;; ============================================================
;;; Eigendecomposition via power iteration + deflation
;;; ============================================================

(defun power-iteration (matrix &key (max-iter 2000) (tolerance 1e-12))
  "Find largest eigenvalue and eigenvector via power iteration"
  (let* ((n (matrix-rows matrix))
         (v (conservation-matrix:make-vec n))
         (eigenvalue 0.0d0))
    (dotimes (i n)
      (setf (aref v i) (/ 1.0d0 (+ 1.0d0 (coerce i 'double-float)))))
    (let ((norm (vector-norm v)))
      (when (> norm 1e-30)
        (dotimes (i n)
          (setf (aref v i) (/ (aref v i) norm)))))
    (dotimes (iter max-iter
              (values eigenvalue v))
      (let* ((w (matrix-vector-multiply matrix v))
             (norm (vector-norm w)))
        (when (< norm 1e-30) (return-from power-iteration (values 0.0d0 v)))
        (dotimes (i n)
          (setf (aref v i) (/ (aref w i) norm)))
        (let* ((mv (matrix-vector-multiply matrix v))
               (new-eigenvalue (vector-dot v mv)))
          (when (< (abs (- new-eigenvalue eigenvalue)) tolerance)
            (setf eigenvalue new-eigenvalue)
            (return-from power-iteration (values eigenvalue v)))
          (setf eigenvalue new-eigenvalue))))))

(defun eigendecompose (lap &key (k nil))
  "Compute eigenvalues and eigenvectors using shifted power iteration + deflation.
Returns (values eigenvalues eigenvectors)."
  (let* ((n (matrix-rows lap))
         (nk (or k n))
         (eigenvalues (conservation-matrix:make-vec nk))
         (eigenvectors nil))
    (let ((shift 0.0d0))
      (dotimes (i n)
        (let ((d (matrix-ref lap i i)))
          (when (> d shift) (setf shift d))))
      (let ((M (make-matrix n n)))
        (dotimes (i n)
          (dotimes (j n)
            (let ((val (- (matrix-ref lap i j))))
              (matrix-set M i j (if (= i j) (+ val shift) val)))))
        (let ((R (matrix-copy M)))
          (dotimes (ev nk)
            (multiple-value-bind (lambda-m v)
                (power-iteration R :max-iter 2000 :tolerance 1e-12)
              (let ((lambda-l (- shift lambda-m)))
                (setf (aref eigenvalues ev) lambda-l)
                (push v eigenvectors)
                (dotimes (i n)
                  (dotimes (j n)
                    (matrix-set R i j
                                (- (matrix-ref R i j)
                                   (* lambda-m (aref v i) (aref v j)))))))))
          (let ((pairs (sort (loop for i below nk
                                   collect (cons (aref eigenvalues i)
                                                 (nth i eigenvectors)))
                             #'< :key #'car)))
            (values (let ((sorted-ev (conservation-matrix:make-vec nk)))
                      (dotimes (i nk sorted-ev)
                        (setf (aref sorted-ev i) (car (nth i pairs)))))
                    (mapcar #'cdr pairs))))))))

;;; ============================================================
;;; Conservation analysis (numerical)
;;; ============================================================

(defun conservation-ratio (eigenvector attributes)
  "Compute conservation ratio for one eigenvector against attributes.
Ratio = variance of gradient of projected attribute."
  (let* ((n (length attributes))
         (projection (conservation-matrix:make-vec n)))
    (dotimes (i n)
      (setf (aref projection i) (* (aref attributes i) (aref eigenvector i))))
    (when (< n 2) (return-from conservation-ratio 0.0d0))
    (let* ((grad-len (1- n))
           (gradient (conservation-matrix:make-vec grad-len)))
      (dotimes (i grad-len)
        (setf (aref gradient i) (- (aref projection (1+ i)) (aref projection i))))
      (let ((mean 0.0d0))
        (dotimes (i grad-len)
          (incf mean (aref gradient i)))
        (setf mean (/ mean (coerce grad-len 'double-float)))
        (let ((var 0.0d0))
          (dotimes (i grad-len)
            (let ((d (- (aref gradient i) mean)))
              (incf var (* d d))))
          (/ var (coerce grad-len 'double-float)))))))

(defun spectral-gap (eigenvalues)
  "Largest gap between consecutive eigenvalues"
  (let ((max-gap 0.0d0)
        (n (length eigenvalues)))
    (dotimes (i (1- n) max-gap)
      (let ((gap (- (aref eigenvalues (1+ i)) (aref eigenvalues i))))
        (when (> gap max-gap)
          (setf max-gap gap))))))

(defun cheeger-constant (laplacian fiedler-vector)
  "Approximate Cheeger constant from Fiedler vector"
  (let* ((n (matrix-rows laplacian))
         (in-s (make-array n :element-type 'bit :initial-element 0))
         (cut 0.0d0)
         (vol-s 0.0d0)
         (total-vol 0.0d0))
    (dotimes (i n)
      (when (< (aref fiedler-vector i) 0.0d0)
        (setf (aref in-s i) 1)))
    (dotimes (i n)
      (let ((deg-i (matrix-ref laplacian i i)))
        (incf total-vol deg-i)
        (when (= (aref in-s i) 1)
          (incf vol-s deg-i)
          (dotimes (j n)
            (when (and (/= i j) (= (aref in-s j) 0))
              (incf cut (- (matrix-ref laplacian i j))))))))
    (let ((min-vol (min vol-s (- total-vol vol-s))))
      (if (< min-vol 1e-15)
          0.0d0
          (/ cut min-vol)))))

;;; ============================================================
;;; Symbolic Laplacian
;;; ============================================================

(defun symbolic-laplacian-from-edges (edges n)
  "Build Laplacian with SYMBOLIC edge weights.
EDGES is a list of (from to weight-symbol).
Returns a 2D array of symbolic expressions."
  (let ((L (make-array (list n n) :initial-element nil)))
    (dolist (e edges L)
      (destructuring-bind (u v w) e
        (setf (aref L u u)
              (if (null (aref L u u))
                  w (sym+ (aref L u u) w)))
        (setf (aref L v v)
              (if (null (aref L v v))
                  w (sym+ (aref L v v) w)))
        (setf (aref L u v)
              (if (null (aref L u v))
                  (sym-negate w) (sym- (aref L u v) w)))
        (setf (aref L v u)
              (if (null (aref L v u))
                  (sym-negate w) (sym- (aref L v u) w)))))))

(defun symbolic-matrix-vector-multiply (matrix vector)
  "Symbolic matrix-vector multiplication."
  (let* ((n (length vector))
         (result (make-array n :initial-element 0)))
    (dotimes (i n result)
      (let ((sum 0))
        (dotimes (j n)
          (let ((product (if (or (zerop (aref matrix i j))
                                 (zerop (aref vector j)))
                             0
                             (sym* (aref matrix i j) (aref vector j)))))
            (setf sum (if (zerop sum) product (sym+ sum product)))))
        (setf (aref result i) sum)))))

(defun symbolic-dot-product (u v)
  "Symbolic dot product"
  (let ((sum 0))
    (dotimes (i (length u) sum)
      (let ((product (if (or (zerop (aref u i)) (zerop (aref v i)))
                         0
                         (sym* (aref u i) (aref v i)))))
        (setf sum (if (zerop sum) product (sym+ sum product)))))))

(defun symbolic-conservation (laplacian attribute)
  "Compute conservation ratio SYMBOLICALLY.
Returns an algebraic expression: ||La||²/||a||²"
  (let* ((La (symbolic-matrix-vector-multiply laplacian attribute))
         (numerator (symbolic-dot-product La La))
         (denominator (symbolic-dot-product attribute attribute)))
    (sym-simplify (sym/ numerator denominator))))

;;; ============================================================
;;; Tracker (sliding window anomaly detection)
;;; ============================================================

(defstruct (tracker (:constructor %make-tracker))
  (window-size 10 :type fixnum)
  (history nil :type (or null (simple-array double-float (*))))
  (count 0 :type fixnum)
  (baseline-mean 0.0d0 :type double-float)
  (baseline-std 0.0d0 :type double-float)
  (baseline-set nil :type boolean))

(defun make-tracker (window-size)
  "Create sliding window tracker"
  (%make-tracker :window-size window-size
                :history (make-array window-size :element-type 'double-float
                                     :initial-element 0.0d0)))

(defun tracker-check (tracker)
  "Check current state"
  (if (not (tracker-baseline-set tracker))
      :nominal
      (let* ((latest (aref (tracker-history tracker)
                           (1- (tracker-count tracker))))
             (std (tracker-baseline-std tracker)))
        (if (< std 1e-15)
            :nominal
            (let ((zscore (/ (abs (- latest (tracker-baseline-mean tracker))) std)))
              (cond ((> zscore 3.0) :critical)
                    ((> zscore 2.0) :warning)
                    (t :nominal)))))))

(defun tracker-feed (tracker observation)
  "Feed observation. Returns :nominal, :warning, or :critical"
  (let* ((obs (coerce observation 'double-float))
         (ws (tracker-window-size tracker))
         (hist (tracker-history tracker)))
    (if (< (tracker-count tracker) ws)
        (progn
          (setf (aref hist (tracker-count tracker)) obs)
          (incf (tracker-count tracker)))
        (progn
          (dotimes (i (1- ws))
            (setf (aref hist i) (aref hist (1+ i))))
          (setf (aref hist (1- ws)) obs)))
    (when (and (= (tracker-count tracker) ws)
               (not (tracker-baseline-set tracker)))
      (let ((sum 0.0d0))
        (dotimes (i ws) (incf sum (aref hist i)))
        (setf (tracker-baseline-mean tracker) (/ sum (coerce ws 'double-float)))
        (let ((var 0.0d0))
          (dotimes (i ws)
            (let ((d (- (aref hist i) (tracker-baseline-mean tracker))))
              (incf var (* d d))))
          (setf (tracker-baseline-std tracker) (sqrt (/ var (coerce ws 'double-float))))
          (setf (tracker-baseline-set tracker) t)))
      (return-from tracker-feed :nominal))
    (tracker-check tracker)))

;;; ============================================================
;;; Spectral fingerprint
;;; ============================================================

(defun fingerprint-compute (eigenvalues)
  "Compute hex fingerprint from eigenvalues"
  (with-output-to-string (s)
    (dotimes (i (length eigenvalues))
      (let* ((val (aref eigenvalues i))
             (bits (logand (floor (* 1e15 (abs val))) #xFFFFFFFFFFFFFFFF))
             (h1 (logxor bits (ash bits -33)))
             (h2 (logand (* h1 #xff51afd7ed558ccd) #xFFFFFFFFFFFFFFFF))
             (h3 (logxor h2 (ash h2 -33)))
             (h4 (logand (* h3 #xc4ceb9fe1a85ec53) #xFFFFFFFFFFFFFFFF))
             (h5 (logxor h4 (ash h4 -33))))
        (format s "~16,'0x" h5)))))

(defun fingerprint-compare (fp1 fp2)
  "Compare two fingerprints, return similarity [0,1]"
  (let* ((l1 (length fp1))
         (l2 (length fp2))
         (max-len (max l1 l2)))
    (if (= max-len 0)
        1.0d0
        (let ((matches 0)
              (min-len (min l1 l2)))
          (dotimes (i min-len)
            (when (char= (char fp1 i) (char fp2 i))
              (incf matches)))
          (/ (coerce matches 'double-float) (coerce max-len 'double-float))))))

;;; ============================================================
;;; Full analysis report
;;; ============================================================

(defstruct conservation-report
  (spectral-gap 0.0d0 :type double-float)
  (cheeger 0.0d0 :type double-float)
  (eigenvalues nil :type (or null (simple-array double-float (*))))
  (ratios nil :type list)
  (fingerprint "" :type string)
  (graph-info nil :type list))

(defun full-analysis (g)
  "Run full conservation spectral analysis on graph"
  (let* ((lap (build-laplacian g))
         (n (graph-n g)))
    (multiple-value-bind (eigenvalues eigenvectors)
        (eigendecompose lap)
      (let* ((attrs (conservation-matrix:make-vec n))
             (ratios nil))
        (dotimes (i n)
          (setf (aref attrs i) (vertex-attribute (aref (graph-vertices g) i))))
        (dotimes (k n)
          (push (cons k (conservation-ratio (nth k eigenvectors) attrs))
                ratios))
        (setf ratios (nreverse ratios))
        (let ((cheeger (if (>= n 2)
                           (cheeger-constant lap (second eigenvectors))
                           0.0d0)))
          (make-conservation-report
           :spectral-gap (spectral-gap eigenvalues)
           :cheeger cheeger
           :eigenvalues eigenvalues
           :ratios ratios
           :fingerprint (fingerprint-compute eigenvalues)
           :graph-info (list :vertices n
                             :edges (length (graph-edges g))
                             :spectral-gap (spectral-gap eigenvalues)
                             :cheeger cheeger)))))))

(defun analyze (g)
  "Convenience: run analysis and print report"
  (let ((report (full-analysis g)))
    (format t "~%=== Conservation Spectral Analysis ===~%")
    (format t "Vertices: ~d, Edges: ~d~%"
            (getf (conservation-report-graph-info report) :vertices)
            (getf (conservation-report-graph-info report) :edges))
    (format t "Spectral Gap: ~f~%" (conservation-report-spectral-gap report))
    (format t "Cheeger Constant: ~f~%" (conservation-report-cheeger report))
    (format t "Eigenvalues: ~{~f~^, ~}~%"
            (coerce (conservation-report-eigenvalues report) 'list))
    (format t "Conservation Ratios:~%")
    (dolist (r (conservation-report-ratios report))
      (format t "  eigenvector ~d: ~f~%" (car r) (cdr r)))
    (format t "Fingerprint: ~a~%" (conservation-report-fingerprint report))
    (format t "===========================~%~%")
    report))
