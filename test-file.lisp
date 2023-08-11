(format t "Hello test file~%")
(ql:quickload :cl-mpm/examples/slump)
(cl-mpm/examples/slump::mpi-run 1)
