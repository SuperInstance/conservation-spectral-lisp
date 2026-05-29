(load "matrix.lisp")
(load "symbolic.lisp")
(load "conservation.lisp")

(defpackage :ct (:use :cl)
  (:shadowing-import-from :conservation-spectral 
    :make-graph :graph-add-edge :build-laplacian 
    :graph-edges :graph-vertices :graph-n
    :vertex-attribute :edge-from :edge-to :edge-weight
    :vertex :edge)
  (:shadowing-import-from :conservation-matrix 
    :make-matrix :matrix-ref :matrix-set :make-vec
    :matrix-vector-multiply :vector-dot))

(in-package :ct)

(format t "Creating graph~%")
(let ((g (conservation-spectral:make-graph 3)))
  (format t "Graph created: ~a~%" g)
  (conservation-spectral:graph-add-edge g 0 1 :weight 1.0d0)
  (format t "Added edge~%")
  (conservation-spectral:graph-add-edge g 1 2 :weight 1.0d0)
  (conservation-spectral:graph-add-edge g 0 2 :weight 1.0d0)
  (format t "Edges: ~a~%" (conservation-spectral:graph-edges g))
  (format t "Building Laplacian~%")
  (let ((L (conservation-spectral:build-laplacian g)))
    (format t "L[0][0] = ~f~%" (conservation-matrix:matrix-ref L 0 0))
    (format t "L[0][1] = ~f~%" (conservation-matrix:matrix-ref L 0 1))
    (format t "L[1][1] = ~f~%" (conservation-matrix:matrix-ref L 1 1))))

(sb-ext:exit)
