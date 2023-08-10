(ql:quickload :cl-mpm-worker)
(ql:quickload "cl-mpm/examples/slump")

(defun cl-mpm-worker::primary-main ()
  (format t "Running MPI with ~D jobs~%"  (cl-mpi::mpi-comm-size))
  (cl-mpm/examples/slump::mpi-run (cl-mpi::mpi-comm-size)))

(cl-mpm-worker::build)
