;;;; symbolic.lisp — Symbolic algebra engine for Conservation Spectral SDK
;;;; The unique advantage: Lisp treats CODE as DATA (S-expressions)
;;;; We manipulate algebraic expressions the same way we manipulate lists.

(defpackage :conservation-symbolic
  (:use :cl)
  (:export :sym :sym-var :sym-num :sym-expr :sym+
           :sym* :sym- :sym/ :sym-expt :sym-negate
           :sym-constant-p :sym-variable-p :sym-expr-p
           :sym-args :sym-op
           :sym-simplify :sym-expand :sym-substitute
           :sym->string :sym->latex
           :sym-diff :sym-eval
           :collect-variables))

(in-package :conservation-symbolic)

;;; === Symbolic expression representation ===
;;; S-expressions ARE the symbolic expressions!
;;; (+ (* x y) 3) is both valid Lisp and valid algebra.

(defun sym-constant-p (expr)
  "Is this a numeric constant?"
  (numberp expr))

(defun sym-variable-p (expr)
  "Is this a symbolic variable?"
  (symbolp expr))

(defun sym-expr-p (expr)
  "Is this a compound expression?"
  (and (consp expr) (member (car expr) '(+ - * / expt neg sqrt sum))))

(defun sym-op (expr)
  "Get operator of symbolic expression"
  (when (consp expr) (car expr)))

(defun sym-args (expr)
  "Get arguments of symbolic expression"
  (when (consp expr) (cdr expr)))

;;; === Constructors ===

(defun sym (&rest args)
  "Create a symbolic expression"
  args)

(defun sym-var (name)
  "Create a symbolic variable"
  (declare (ignore name))
  (gentemp "V"))

(defun sym-num (value)
  "Wrap a number (identity for numbers)"
  value)

(defun sym+ (&rest args)
  "Symbolic addition"
  (let ((nonzero (remove 0 args)))
    (cond ((null nonzero) 0)
          ((= (length nonzero) 1) (first nonzero))
          (t (cons '+ nonzero)))))

(defun sym* (&rest args)
  "Symbolic multiplication"
  (let ((nonzero (remove-if (lambda (x) (and (numberp x) (= x 0))) args)))
    (cond ((null nonzero) 0)
          ((member 0 args) 0)
          ((and (= (length nonzero) 1)) (first nonzero))
          ((member 1 args)
           (apply #'sym* (remove 1 nonzero)))
          (t (cons '* nonzero)))))

(defun sym- (a &optional b)
  "Symbolic subtraction"
  (cond ((null b) (if (numberp a) (- a) (list 'neg a)))
        ((and (numberp a) (numberp b)) (- a b))
        ((and (numberp b) (= b 0)) a)
        (t (list '- a b))))

(defun sym/ (a b)
  "Symbolic division"
  (cond ((and (numberp a) (numberp b) (not (zerop b))) (/ a b))
        ((and (numberp b) (= b 1)) a)
        (t (list '/ a b))))

(defun sym-expt (base power)
  "Symbolic exponentiation"
  (cond ((and (numberp base) (numberp power)) (expt base power))
        ((and (numberp power) (= power 0)) 1)
        ((and (numberp power) (= power 1)) base)
        (t (list 'expt base power))))

(defun sym-negate (expr)
  "Negate symbolically"
  (if (numberp expr)
      (- expr)
      (list 'neg expr)))

;;; === Simplification ===

(defun sym-simplify (expr)
  "Simplify symbolic expression using algebraic rules"
  (cond
    ;; Atom: return as-is
    ((atom expr) expr)
    ;; Compound: simplify then apply rules
    (t
     (let ((op (car expr))
           (args (mapcar #'sym-simplify (cdr expr))))
       (case op
         ;; Negation
         (neg
          (let ((a (first args)))
            (cond ((numberp a) (- a))
                  ((and (consp a) (eq (car a) 'neg)) (second a))
                  (t (list 'neg a)))))
         ;; Addition
         (+
          (let ((flat (flatten-nary '+ args)))
            (simplify-sum flat)))
         ;; Subtraction
         (-
          (let ((a (first args))
                (b (second args)))
            (cond ((null b) (sym-simplify (list 'neg a)))
                  ((and (numberp a) (numberp b)) (- a b))
                  ((and (numberp b) (= b 0)) a)
                  (t (list '- a b)))))
         ;; Multiplication
         (*
          (let ((flat (flatten-nary '* args)))
            (simplify-product flat)))
         ;; Division
         (/
          (let ((num (first args))
                (den (second args)))
            (cond ((and (numberp num) (numberp den) (not (zerop den)))
                   (/ num den))
                  ((and (numberp den) (= den 1)) num)
                  (t (list '/ num den)))))
         ;; Exponentiation
         (expt
          (let ((base (first args))
                (power (second args)))
            (cond ((and (numberp base) (numberp power)) (expt base power))
                  ((and (numberp power) (= power 0)) 1)
                  ((and (numberp power) (= power 1)) base)
                  (t (list 'expt base power)))))
         ;; Default: rebuild
         (otherwise (cons op args)))))))

(defun flatten-nary (op args)
  "Flatten nested applications of associative operator"
  (let ((result nil))
    (dolist (a args)
      (if (and (consp a) (eq (car a) op))
          (setf result (nconc result (cdr a)))
          (push a result)))
    (nreverse result)))

(defun simplify-sum (terms)
  "Simplify a sum by collecting constants"
  (let ((constants 0)
        (vars nil))
    (dolist (t_ terms)
      (if (numberp t_)
          (incf constants t_)
          (push t_ vars)))
    (let ((nonzero-vars (nreverse vars)))
      (cond ((and (null nonzero-vars) (= constants 0)) 0)
            ((null nonzero-vars) constants)
            ((= constants 0)
             (if (= (length nonzero-vars) 1)
                 (first nonzero-vars)
                 (cons '+ nonzero-vars)))
            (t (cons '+ (append nonzero-vars (list constants))))))))

(defun simplify-product (terms)
  "Simplify a product by collecting constants"
  (let ((constants 1)
        (vars nil)
        (has-zero nil))
    (dolist (t_ terms)
      (cond ((and (numberp t_) (= t_ 0)) (setf has-zero t))
            ((numberp t_) (setf constants (* constants t_)))
            (t (push t_ vars))))
    (if has-zero
        0
        (let ((nonzero-vars (nreverse vars)))
          (cond ((and (= constants 1) (null nonzero-vars)) 1)
                ((null nonzero-vars) constants)
                ((= constants 1)
                 (if (= (length nonzero-vars) 1)
                     (first nonzero-vars)
                     (cons '* nonzero-vars)))
                (t (cons '* (append nonzero-vars (list constants)))))))))

;;; === Expansion ===

(defun sym-expand (expr)
  "Expand products: a*(b+c) → a*b + a*c"
  (cond ((atom expr) expr)
        (t
         (let ((op (car expr))
               (args (cdr expr)))
           (case op
             (*
              (let ((expanded-args (mapcar #'sym-expand args)))
                (expand-product expanded-args)))
             (+
              (apply #'sym+ (mapcar #'sym-expand args)))
             (otherwise
              (cons op (mapcar #'sym-expand args))))))))

(defun expand-product (factors)
  "Distribute multiplication over addition"
  (cond ((null factors) 1)
        ((= (length factors) 1) (first factors))
        (t
         (let ((a (first factors))
               (rest-expanded (expand-product (rest factors))))
           (distribute a rest-expanded)))))

(defun distribute (a b)
  "Distribute a over b if b is a sum"
  (cond ((and (consp b) (eq (car b) '+))
         (apply #'sym+ (mapcar (lambda (term) (sym* a term)) (cdr b))))
        ((and (consp a) (eq (car a) '+))
         (apply #'sym+ (mapcar (lambda (term) (sym* term b)) (cdr a))))
        (t (sym* a b))))

;;; === Substitution ===

(defun sym-substitute (expr var replacement)
  "Substitute all occurrences of VAR with REPLACEMENT"
  (cond ((equal expr var) replacement)
        ((atom expr) expr)
        (t (cons (car expr)
                 (mapcar (lambda (a) (sym-substitute a var replacement))
                         (cdr expr))))))

;;; === Differentiation ===

(defun sym-diff (expr var)
  "Symbolic differentiation with respect to VAR"
  (sym-simplify
   (cond
     ;; Constant
     ((numberp expr) 0)
     ;; Variable
     ((symbolp expr)
      (if (eq expr var) 1 0))
     ;; Compound
     ((consp expr)
      (let ((op (car expr))
            (args (cdr expr)))
        (case op
          (+ (apply #'sym+ (mapcar (lambda (a) (sym-diff a var)) args)))
          (- (apply #'sym- (mapcar (lambda (a) (sym-diff a var)) args)))
          (* (let ((u (first args))
                   (v (second args)))
               ;; Product rule: d(uv) = u'dv + v'du
               ;; For n-ary: iterative
               (if (= (length args) 2)
                   (sym+ (sym* (sym-diff u var) v)
                         (sym* u (sym-diff v var)))
                   (sym-diff (sym* u (apply #'sym* (rest args))) var))))
          (/ (let ((u (first args))
                   (v (second args)))
               ;; Quotient rule
               (sym/ (sym- (sym* (sym-diff u var) v)
                           (sym* u (sym-diff v var)))
                     (sym* v v))))
          (expt (let ((base (first args))
                      (power (second args)))
                  ;; d/dx(a^n) = n * a^(n-1) * da/dx (power rule)
                  (sym* power
                        (sym-expt base (sym- power 1))
                        (sym-diff base var))))
          (neg (sym-negate (sym-diff (first args) var)))
          (otherwise 0))))
     (t 0))))

;;; === Evaluation ===

(defun sym-eval (expr bindings)
  "Evaluate symbolic expression with variable bindings (alist)"
  (cond ((numberp expr) expr)
        ((symbolp expr)
         (let ((binding (assoc expr bindings)))
           (if binding (cdr binding) expr)))
        ((consp expr)
         (let ((op (car expr))
               (args (mapcar (lambda (a) (sym-eval a bindings)) (cdr expr))))
           ;; Check if all args are numbers
           (when (every #'numberp args)
             (case op
               (+ (apply #'+ args))
               (- (if (= (length args) 1) (- (first args)) (apply #'- args)))
               (* (apply #'* args))
               (/ (apply #'/ args))
               (expt (apply #'expt args))
               (neg (- (first args)))
               (sqrt (sqrt (first args)))
               (otherwise (cons op args)))))))
  )

;;; === Collect variables ===

(defun collect-variables (expr)
  "Collect all free variables in a symbolic expression"
  (cond ((numberp expr) nil)
        ((symbolp expr) (list expr))
        ((consp expr)
         (remove-duplicates
          (reduce #'append (mapcar #'collect-variables (cdr args))
                  :initial-value nil)))))

;;; === Pretty printing ===

(defun sym->string (expr)
  "Convert symbolic expression to human-readable string"
  (cond ((numberp expr) (format nil "~f" expr))
        ((symbolp expr) (string expr))
        ((consp expr)
         (let ((op (car expr))
               (args (cdr expr)))
           (case op
             (neg (format nil "(-~a)" (sym->string (first args))))
             (+ (format nil "(~{~a~^ + ~})" (mapcar #'sym->string args)))
             (- (if (= (length args) 1)
                    (format nil "(0 - ~a)" (sym->string (first args)))
                    (format nil "(~a - ~a)" (sym->string (first args))
                            (sym->string (second args)))))
             (* (format nil "(~{~a~^ · ~})" (mapcar #'sym->string args)))
             (/ (format nil "(~a / ~a)" (sym->string (first args))
                        (sym->string (second args))))
             (expt (format nil "(~a)^~a" (sym->string (first args))
                           (sym->string (second args))))
             (sum (format nil "Σ[~a](~a)" (sym->string (first args))
                          (sym->string (second args))))
             (otherwise (format nil "(~a ~{~a~^ ~})"
                                op (mapcar #'sym->string args))))))))

(defun sym->latex (expr)
  "Convert symbolic expression to LaTeX string"
  (cond ((numberp expr) (format nil "~f" expr))
        ((symbolp expr) (string expr))
        ((consp expr)
         (let ((op (car expr))
               (args (cdr expr)))
           (case op
             (neg (format nil "-~a" (sym->latex (first args))))
             (+ (format nil "~{~a~^ + ~}" (mapcar #'sym->latex args)))
             (- (if (= (length args) 1)
                    (format nil "-~a" (sym->latex (first args)))
                    (format nil "\\frac{~a - ~a}{}"
                            (sym->latex (first args))
                            (sym->latex (second args)))))
             (* (format nil "~{~a~^ \\cdot ~}" (mapcar #'sym->latex args)))
             (/ (format nil "\\frac{~a}{~a}"
                        (sym->latex (first args))
                        (sym->latex (second args))))
             (expt (format nil "{~a}^{~a}"
                           (sym->latex (first args))
                           (sym->latex (second args))))
             (sqrt (format nil "\\sqrt{~a}" (sym->latex (first args))))
             (sum (format nil "\\sum ~a" (sym->latex (second args))))
             (otherwise (format nil "\\text{~a}(~{~a~^, ~})"
                                op (mapcar #'sym->latex args))))))))
