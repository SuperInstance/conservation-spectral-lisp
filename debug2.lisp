(load "matrix.lisp")
(load "symbolic.lisp")
(load "conservation.lisp")

(format t "make-graph lambda: ~a~%" 
  (function-lambda-expression (symbol-function 'conservation-spectral:make-graph)))

;; Try calling directly
(handler-case 
  (let ((g (conservation-spectral:make-graph 3)))
    (format t "graph: ~a~%" g))
  (error (c) (format t "ERROR: ~a~%" c)))

;; Try %make-graph
(format t "%make-graph fbound: ~a~%" (fboundp 'conservation-spectral:%make-graph))

(sb-ext:exit)
