;(restrict-compiler-policy 'speed  0 0)
;(restrict-compiler-policy 'debug  3 3)
;(restrict-compiler-policy 'safety 3 3)
(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

;; (ql:quickload :cl-mpm)
;; (ql:quickload :cl-mpm/damage)
;; (ql:quickload :cl-mpm/mpi)
(ql:quickload :cl-mpm/examples/slump)

(asdf:compile-system :cl-mpm :force t)
(asdf:compile-system :cl-mpm/setup :force t)
(asdf:compile-system :cl-mpm/particle :force t)
(asdf:compile-system :cl-mpm/mesh :force t)
(asdf:compile-system :cl-mpm/constitutive :force t)
(asdf:compile-system :cl-mpm/damage :force t)
(asdf:compile-system :cl-mpm/utils :force t)
(asdf:compile-system :cl-mpm/bc :force t)
(asdf:compile-system :cl-mpm/forces :force t)
(asdf:compile-system :cl-mpm/ext :force t)
(asdf:compile-system :cl-mpm/output :force t)
(asdf:compile-system :cl-mpm/fastmath :force t)
(asdf:compile-system :cl-mpm/mpi :force t)
 
(ql:quickload :cl-mpm-worker)
(ql:quickload :cl-mpm/examples/slump)
(ql:quickload :cl-mpm/mpi)
(in-package :cl-mpm-worker)
;; (in-package :cl-mpm/examples/chalk)
;; ;(asdf:compile-system :cl-mpm :force t)
;; ;(asdf:compile-system :cl-mpm/damage :force t)
;; ;(asdf:compile-system :cl-mpm/mpi :force t)
;; (asdf:compile-system :cl-mpm/examples/chalk :force t)
;; (in-package :cl-mpm-worker)

(sb-ext:save-lisp-and-die
   "mpi-worker"
    :executable t
    :toplevel #'main
    :save-runtime-options t)
(uiop:quit)
