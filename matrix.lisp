;;;; matrix.lisp — Numerical matrix operations for Conservation Spectral SDK
;;;; Dense row-major matrix math, power iteration, eigendecomposition

(defpackage :conservation-matrix
  (:use :cl)
  (:export :make-matrix :matrix-dim :matrix-ref :matrix-set
           :matrix-copy :matrix-identity
           :matrix-add :matrix-sub :matrix-scale
           :matrix-multiply :matrix-transpose
           :matrix-vector-multiply :vector-norm :vector-scale
           :vector-dot :vector-add :vector-sub
           :matrix-trace
           :make-vec :matrix-rows :matrix-cols :matrix-data))

(in-package :conservation-matrix)

;;; === Matrix creation and access ===

(defstruct (matrix (:constructor %make-matrix))
  (rows 0 :type fixnum)
  (cols 0 :type fixnum)
  (data nil :type (simple-array double-float (*))))

(defun make-matrix (rows cols &key (initial-element 0.0d0))
  "Create a rows×cols matrix initialized to initial-element"
  (%make-matrix :rows rows :cols cols
                :data (make-array (* rows cols)
                                  :element-type 'double-float
                                  :initial-element (coerce initial-element 'double-float))))

(defun matrix-dim (m)
  "Return (rows . cols)"
  (cons (matrix-rows m) (matrix-cols m)))

(defun matrix-ref (m i j)
  "Access element at row i, col j"
  (aref (matrix-data m) (+ (* i (matrix-cols m)) j)))

(defun matrix-set (m i j val)
  "Set element at row i, col j"
  (setf (aref (matrix-data m) (+ (* i (matrix-cols m)) j))
        (coerce val 'double-float)))

(defun matrix-copy (m)
  "Deep copy of matrix"
  (let ((copy (make-matrix (matrix-rows m) (matrix-cols m))))
    (replace (matrix-data copy) (matrix-data m))
    copy))

(defun matrix-identity (n)
  "Create n×n identity matrix"
  (let ((m (make-matrix n n)))
    (dotimes (i n m)
      (matrix-set m i i 1.0d0))))

;;; === Matrix arithmetic ===

(defun matrix-add (a b)
  "Element-wise addition A + B"
  (assert (and (= (matrix-rows a) (matrix-rows b))
               (= (matrix-cols a) (matrix-cols b)))
          (a b) "Dimension mismatch")
  (let ((result (make-matrix (matrix-rows a) (matrix-cols a))))
    (dotimes (i (* (matrix-rows a) (matrix-cols a)) result)
      (setf (aref (matrix-data result) i)
            (+ (aref (matrix-data a) i) (aref (matrix-data b) i))))))

(defun matrix-sub (a b)
  "Element-wise subtraction A - B"
  (assert (and (= (matrix-rows a) (matrix-rows b))
               (= (matrix-cols a) (matrix-cols b)))
          (a b) "Dimension mismatch")
  (let ((result (make-matrix (matrix-rows a) (matrix-cols a))))
    (dotimes (i (* (matrix-rows a) (matrix-cols a)) result)
      (setf (aref (matrix-data result) i)
            (- (aref (matrix-data a) i) (aref (matrix-data b) i))))))

(defun matrix-scale (m scalar)
  "Scale all elements by scalar"
  (let ((s (coerce scalar 'double-float))
        (result (make-matrix (matrix-rows m) (matrix-cols m))))
    (dotimes (i (* (matrix-rows m) (matrix-cols m)) result)
      (setf (aref (matrix-data result) i)
            (* s (aref (matrix-data m) i))))))

(defun matrix-multiply (a b)
  "Matrix multiplication A × B"
  (assert (= (matrix-cols a) (matrix-rows b)) (a b) "Dimension mismatch")
  (let* ((m (matrix-rows a))
         (n (matrix-cols b))
         (k (matrix-cols a))
         (result (make-matrix m n)))
    (dotimes (i m result)
      (dotimes (j n)
        (let ((sum 0.0d0))
          (dotimes (p k)
            (incf sum (* (matrix-ref a i p) (matrix-ref b p j))))
          (matrix-set result i j sum))))))

(defun matrix-transpose (m)
  "Transpose matrix"
  (let ((result (make-matrix (matrix-cols m) (matrix-rows m))))
    (dotimes (i (matrix-rows m) result)
      (dotimes (j (matrix-cols m))
        (matrix-set result j i (matrix-ref m i j))))))

(defun matrix-trace (m)
  "Trace of a square matrix"
  (assert (= (matrix-rows m) (matrix-cols m)) (m) "Not square")
  (let ((sum 0.0d0))
    (dotimes (i (matrix-rows m) sum)
      (incf sum (matrix-ref m i i)))))

;;; === Vector operations (1D arrays of double-float) ===

(defun make-vec (n &key (initial-element 0.0d0))
  (make-array n :element-type 'double-float
              :initial-element (coerce initial-element 'double-float)))

(defun vector-dot (u v)
  "Dot product of two vectors"
  (let ((sum 0.0d0))
    (dotimes (i (length u) sum)
      (incf sum (* (aref u i) (aref v i))))))

(defun vector-norm (v)
  "Euclidean norm"
  (sqrt (vector-dot v v)))

(defun vector-scale (v scalar)
  "Scale vector by scalar"
  (let ((s (coerce scalar 'double-float))
        (result (make-vec (length v))))
    (dotimes (i (length v) result)
      (setf (aref result i) (* s (aref v i))))))

(defun vector-add (a b)
  (let ((result (make-vec (length a))))
    (dotimes (i (length a) result)
      (setf (aref result i) (+ (aref a i) (aref b i))))))

(defun vector-sub (a b)
  (let ((result (make-vec (length a))))
    (dotimes (i (length a) result)
      (setf (aref result i) (- (aref a i) (aref b i))))))

;;; === Matrix-vector multiply ===

(defun matrix-vector-multiply (m v)
  "Compute M × v"
  (assert (= (matrix-cols m) (length v)) (m v) "Dimension mismatch")
  (let ((result (make-vec (matrix-rows m))))
    (dotimes (i (matrix-rows m) result)
      (let ((sum 0.0d0))
        (dotimes (j (matrix-cols m))
          (incf sum (* (matrix-ref m i j) (aref v j))))
        (setf (aref result i) sum)))))

;;; === Print matrix (for debugging) ===

(defmethod print-object ((m matrix) stream)
  (print-unreadable-object (m stream :type t)
    (format stream "~d×~d" (matrix-rows m) (matrix-cols m))))
