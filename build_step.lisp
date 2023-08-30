(ql:quickload :cl-mpm-worker)
(in-package :cl-mpm-worker)
;(ql:quickload "cl-mpm/examples/slump")
;(require 'cl-mpm-worker)

;(defun cl-mpm-worker::primary-main ()
;  (format t "Running MPI with ~D jobs~%"  (cl-mpi::mpi-comm-size))
;  (cl-mpm/examples/slump::mpi-run (cl-mpi::mpi-comm-size)))


;(build)
;(asdf:make :cl-mpm-worker)
;lfarm-server 
;(ql:quickload :lfarm-server)
;(defun main (&rest args)
;  (print "hello")
;   (lfarm-server:start-server "127.0.0.1" 11111)
;  )
(sb-ext:save-lisp-and-die
   "mpi-worker"
    :executable t
    :toplevel #'main
    :save-runtime-options t)
(uiop:quit)
