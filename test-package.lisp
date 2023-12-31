(format t "Hello test file~%")
(sb-ext:restrict-compiler-policy 'speed  0 0)
(sb-ext:restrict-compiler-policy 'debug  3 3)
(sb-ext:restrict-compiler-policy 'safety 3 3)
(ql:quickload :cl-mpi :silent t)
(ql:quickload :cl-mpm :silent t)
(ql:quickload :cl-mpm/setup :silent t)
(ql:quickload :cl-mpm/particle :silent t)
(ql:quickload :cl-mpm/output :silent t)
(ql:quickload :cl-mpm/buoyancy :silent t)
(ql:quickload :cl-mpm/penalty :silent t)
(ql:quickload :cl-mpm/eigenerosion :silent t)
(ql:quickload :cl-mpm/damage :silent t)
(ql:quickload :cl-mpm/mpi :silent t)
(ql:quickload :vgplot :silent t)
(ql:quickload :magicl :silent t)
;; (ql:quickload :cl-mpm/examples/slump)
;; (in-package :cl-mpm/examples/slump)

(declaim (optimize (debug 3) (safety 3) (speed 0)))
(defparameter *sim* nil)
;;shim to make it work in a repl without mpi
;; (defun cl-mpi::mpi-comm-rank ()
;;   0)
;; (defun cl-mpi::mpi-comm-size ()
;;   4)



(defun length-from-def (sim mp dim)
  (let* ((mp-scale 2)
         (h-initial (magicl:tref (cl-mpm/particle::mp-domain-size mp) dim 0)))
    h-initial
    ))
(defun max-stress (mp)
  (declare (optimize (speed 0) (debug 3)))
  (multiple-value-bind (l v) (magicl:eig (cl-mpm::voight-to-matrix (cl-mpm/particle:mp-stress mp)))
    (cl-mpm/particle::mp-time-averaged-visc mp)
    )
  )
(defun local-dist (sim mp)
  (with-accessors ((ll cl-mpm/particle::mp-true-local-length)) mp
    (with-accessors ((llm cl-mpm/particle::mp-true-local-length)) *dist-mp*
      (cl-mpm/damage::diff-squared mp *dist-mp*)
      )))
(declaim (notinline plot))
(defun plot (sim &optional (plot :damage))
  (declare (optimize (speed 0) (debug 3)))
  (vgplot:format-plot t "set palette defined (0 'blue', 1 'red')")
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:format-plot t "set object 1 rect from 0,0 to ~f,~f fc rgb 'blue' fs transparent solid 0.5 noborder behind" ms-x *water-height*))

  (multiple-value-bind (x y c stress-y lx ly e density temp vx ybar dist)
    (loop for mp across (cl-mpm:sim-mps sim)
          collect (magicl:tref (cl-mpm::mp-position mp) 0 0) into x
          collect (magicl:tref (cl-mpm::mp-position mp) 1 0) into y
          collect (magicl:tref (cl-mpm::mp-velocity mp) 0 0) into vx
          collect (length-from-def sim mp 0) into lx
          collect (length-from-def sim mp 1) into ly
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage) (cl-mpm/particle:mp-damage mp) 0) into c
          collect (if (slot-exists-p mp 'cl-mpm/particle::damage-ybar) (cl-mpm/particle::mp-damage-ybar mp) 0) into ybar
          collect (if (slot-exists-p mp 'cl-mpm/particle::temperature) (cl-mpm/particle:mp-temperature mp) 0) into temp
          collect (if (slot-exists-p mp 'cl-mpm/particle::strain-energy-density) (cl-mpm/particle::mp-strain-energy-density mp) 0) into e
          collect (/ (cl-mpm/particle:mp-mass mp) (cl-mpm/particle:mp-volume mp)) into density
          ;; collect (cl-mpm/particle:mp-volume mp) into density
          collect (max-stress mp) into stress-y
          collect 0d0 into dist
          finally (return (values x y c stress-y lx ly e density temp vx ybar dist)))

    (let* ((node-x '())
           (node-y '())
           (node-c '())
           (mesh (cl-mpm:sim-mesh *sim*))
           (nodes (cl-mpm/mesh:mesh-nodes mesh))
           )
      (dotimes (i (array-total-size nodes))
        (let ((n (row-major-aref nodes i)))
          (with-accessors ((index cl-mpm/mesh:node-index)
                           (boundary cl-mpm/mesh::node-boundary-node)
                           (boundary-s cl-mpm/mesh::node-boundary-scalar)
                           (active cl-mpm/mesh::node-active)
                           )

              n
            (let ((n-ratio (/ (cl-mpm/mesh::node-volume n) (cl-mpm/mesh::node-volume-true n))))
              (when active
                (if boundary
                    (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
                      (push x node-x)
                      (push y node-y)
                      (push 2 node-c)
                      )
                    (destructuring-bind (x y) (cl-mpm/mesh:index-to-position mesh index)
                      (push x node-x)
                      (push y node-y)
                      (push 0 node-c)
                      ))
                )))))
      (cond
        ((eq plot :point)
         (if node-x
             (progn
               (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min node-c) (+ 0.01 (apply #'max node-c)))
               (vgplot:plot
                x y lx ly ";;with ellipses"
                node-x node-y node-c ";;with points pt 7 lc palette"))
             (vgplot:plot
              x y lx ly ";;with ellipses")))
        ((eq plot :contact)
         (let* ((contact-points (cl-mpm/penalty::collect-contact-points-bc (cl-mpm:sim-mesh *sim*) (cl-mpm:sim-mps *sim*) *floor-bc*))
                (c-x (mapcar (lambda (p) (magicl:tref p 0 0)) contact-points))
                (c-y (mapcar (lambda (p) (magicl:tref p 1 0)) contact-points))
                (c-c (mapcar (lambda (p) 1d0) contact-points))
           )
           (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-40 (apply #'max c)))
         (if c-x
             (vgplot:plot
              x y c ";;with points pt 7 lc palette"
              c-x c-y ";;with points pt 7")
             (vgplot:plot
              x y c ";;with points pt 7 lc palette"))))
        ((eq plot :damage)
         (vgplot:format-plot t "set style fill solid")
         (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max c)))
         (vgplot:plot x y lx ly c ";;with ellipses lc palette"))
        ((eq plot :dist)
         (vgplot:format-plot t "set cbrange [~f:~f]" 0d0 (+ 1e-6 (apply #'max dist)))
         (vgplot:plot x y dist ";;with points pt 7 lc palette"))
        ((eq plot :velocity)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min vx) (+ 0.01 (apply #'max vx)))
         (vgplot:plot x y vx ";;with points pt 7 lc palette"))
        ((eq plot :temperature)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min temp) (+ 0.01 (apply #'max temp)))
         (vgplot:plot x y temp ";;with points pt 7 lc palette"))
        ((eq plot :energy)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min e) (+ 0.01 (apply #'max e)))
         (vgplot:plot x y e ";;with points pt 7 lc palette"))
        ((eq plot :stress)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 1e-40 (apply #'max stress-y)))
         (vgplot:plot x y stress-y ";;with points pt 7 lc palette"))
        ((eq plot :damage-ybar)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min ybar) (+ 1e-20 (apply #'max ybar)))
         (vgplot:plot x y ybar ";;with points pt 7 lc palette"))
        ((eq plot :deformed)
         ;; (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min stress-y) (+ 0.01 (apply #'max stress-y)))
         (vgplot:plot x y lx ly ";;with ellipses"))
        ((eq plot :density)
         (vgplot:format-plot t "set cbrange [~f:~f]" (apply #'min density) (+ 0.01 (apply #'max density)))
         (vgplot:plot x y e ";;with points pt 7 lc palette")
         )))
    )
  ;; (vgplot:format-plot t "replot (~f*x + ~f)~%" *sliding-slope* *sliding-offset*)
  (let* ((ms (cl-mpm/mesh:mesh-mesh-size (cl-mpm:sim-mesh sim)))
         (ms-x (first ms))
         (ms-y (second ms))
         )
    (vgplot:axis (list 0 ms-x
                       0 ms-y))
    (vgplot:format-plot t "set size ratio ~f" (/ ms-y ms-x)))
    (let ((h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*))))
      (vgplot:format-plot t "set ytics ~f" h)
      (vgplot:format-plot t "set xtics ~f" h))

  (vgplot:replot))

(declaim (notinline setup))
(declaim (notinline setup-test-column))
(defun setup-test-column (sim-type
                          size block-size block-offset slope &optional (e-scale 1d0) (mp-scale 1d0)
                          &rest mp-args)
  (declare (optimize (speed 0)))
  (let* ((sim (cl-mpm/setup::make-block (/ 1 e-scale)
                                        (mapcar (lambda (s) (* s e-scale)) size)
                                        #'cl-mpm/shape-function:make-shape-function-bspline
                                        sim-type
                                        ;'cl-mpm/damage::mpm-sim-damage
                                        ;'cl-mpm/mpi::mpm-sim-mpi-stress
                                        ))

         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         ;(e-scale 1)
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (angle 0d0)
         (density *ice-density*)
         ;; (mass (/ (* 900 h-x h-y) (expt mp-scale 2)))
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (progn
      (let ((block-position
              (mapcar #'+ (list (* h-x (- (+ (/ 1 (* 2 mp-scale))) 0))
                                (* h-y (+ (/ 1d0 (* 2d0 mp-scale)))))
                      block-offset)))
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-sloped-list
                block-position
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                'cl-mpm/particle::particle-viscoplastic-damage
                :E 1d9
                :nu 0.3250d0

                :visc-factor 111d6
                :visc-power 3d0

                :initiation-stress 0.7d6
                :damage-rate 1d-5
                :critical-damage 0.50d0
                :local-length 50d0
                :local-length-damaged 0.1d0
                :damage 0.0d0

                :gravity -9.8d0
                ;:index 0

                :slope slope
                )))
        )
      (let ((mass-scale 1d6))
        (setf (cl-mpm::sim-mass-scale sim) mass-scale)
        (setf (cl-mpm:sim-damping-factor sim)
              ;; 0.0d0
              100d0
              )
        )
      (setf (cl-mpm:sim-mass-filter sim) 1d-15)
      (setf (cl-mpm::sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) t)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm:sim-dt sim) 1d-4)
      (setf (cl-mpm:sim-bcs sim) (make-array 0))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0)))
             ))

      ;; (format t "Bottom level ~F~%" h-y)
       (let* ((terminus-size (+ (second block-size) (* slope (first block-size))))
              (ocean-x 1000)
             (ocean-y (+ h-y (* 0.90d0 0.6d0 terminus-size)))
             (angle 0d0)
             )

      ;;   (loop for mp across (cl-mpm:sim-mps sim)
      ;;         do
      ;;            (with-accessors ((pos cl-mpm/particle:mp-position)
      ;;                             (stress cl-mpm/particle::mp-stress-kirchoff)
      ;;                             (undamaged-stress cl-mpm/particle::mp-undamaged-stress)
      ;;                             (stress-cauchy cl-mpm/particle::mp-stress)
      ;;                             )
      ;;                mp
      ;;              (setf stress
      ;;                    (cl-mpm/utils:matrix-to-voight
      ;;                     (magicl:eye 2 :value (* 1d0 (cl-mpm/buoyancy::pressure-at-depth (magicl:tref pos 1 0) ocean-y *water-density*))))
      ;;                    stress-cauchy (magicl:scale stress 1d0)
      ;;                    undamaged-stress (magicl:scale stress 1d0)
      ;;                    )))

      ;;   (format t "Ocean level ~a~%" ocean-y)
      ;;   (defparameter *water-height* ocean-y)
         (defparameter *floor-bc*
           (cl-mpm/penalty::make-bc-penalty-point-normal
            sim
            (cl-mpm/utils:vector-from-list '(0d0 1d0))
            ;; (magicl:from-list (list (sin (- (* pi (/ angle 180d0))))
            ;;                         (cos (+ (* pi (/ angle 180d0))))) '(2 1))
            (magicl:from-list (list 00d0 (+ 1d0 h-y)) '(2 1))
            (* *ice-density* 1d3)
            0.9d0
            ))
         (setf (cl-mpm::sim-bcs-force-list sim)
               (list
                ;(cl-mpm/bc:make-bcs-from-list
                ; (list
                ;  (cl-mpm/buoyancy::make-bc-buoyancy
                ;   sim
                ;   ocean-y
                ;   *water-density*
                ;   )
                ;  ))
                (cl-mpm/bc:make-bcs-from-list
                 (list *floor-bc*)
                 )))
         )
      ;; (let ((normal (magicl:from-list (list (sin (- (* pi (/ angle 180d0))))
      ;;                                       (cos (+ (* pi (/ angle 180d0))))) '(2 1))))
      ;;   (defparameter *sliding-slope* (/ (- (magicl:tref normal 0 0))
      ;;                                    (magicl:tref normal 1 0)
      ;;                                    ))
      ;;   (defparameter *sliding-offset* (- h-y (* 0d0 *sliding-slope*))))
      sim)))

(defparameter *ice-density* 900)
(defparameter *water-density* 1000)
;; (defparameter *ice-density* 900)
;; (defparameter *water-density* 1000)
;Setup
(defun setup (sim-type)
  (declare (optimize (speed 0)))
  (defparameter *run-sim* nil)
  (let* ((mesh-size 50)
         (mps-per-cell 2)
         (slope -0.02)
         (shelf-height 400)
         (shelf-aspect 12)
         (shelf-length (* shelf-height shelf-aspect))
         (shelf-end-height (+ shelf-height (* (- slope) shelf-length)))
         (shelf-height-terminus shelf-height)
         (shelf-height shelf-end-height)
         (offset (list 0 0))
         )
    (defparameter *sim*
      (setup-test-column
        sim-type
        (list
            (+ shelf-length (* 2 shelf-height))
            (+ shelf-height 100))
        (list shelf-length shelf-height)
        (mapcar #'+ offset (list 000 (* 1 mesh-size)))
        slope
        (/ 1 mesh-size) mps-per-cell)
      )

    (format t "Type of sim ~a~%" (type-of *sim*))
    ;;Delete all the plotted frames
    (loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))

  (format t "Estimated dt ~f~%" (cl-mpm:sim-dt *sim*))
  (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
  (defparameter *velocity* '())
  (defparameter *time* '())
  (defparameter *t* 0)
  (defparameter *x* 0d0)
  (defparameter *x-pos* '())
  (defparameter *cfl-max* '())
  (defparameter *sim-step* 0)
  (defparameter *water-height* 0d0)
  ))


(defparameter *run-sim* nil)
(defun calculate-dt (courent c-value target-step)
  (let* (
         (cfl-dt (if (> courent 0d0) (/ c-value (/ courent (cl-mpm:sim-dt *sim*))) nil))
         (new-dt (/ c-value (if (> courent 0d0)
                                (/ courent (cl-mpm:sim-dt *sim*))
                                (/ c-value 1e-6))))
         (max-steps 1000)
         (sub-steps (max (min (floor (/ target-step new-dt)) max-steps) 1)))
    (when (> (floor (/ target-step new-dt)) max-steps)
        (format t "CFL requires more steps than max-steps~%"))
    (format t "C: ~f - steps: ~D - %dt: ~f~%" courent sub-steps new-dt)
    (format t "Cfl derived dt:~f~%" cfl-dt)
                                        ;(setf (cl-mpm:sim-dt *sim*) new-dt)
    (values new-dt sub-steps)))

(defun run ()
  ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk")
  ;;                         *sim*)
  (defparameter *run-sim* t)
    ;(vgplot:close-all-plots)
    ;(vgplot:figure)
  (sleep 1)
 (let* ((target-time 1d2)
         (dt (cl-mpm:sim-dt *sim*))
         (dt-scale 0.1d0)
         (substeps (floor target-time dt))
        (rank (cl-mpi::mpi-comm-rank)))

    (cl-mpm::update-sim *sim*)
    (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
           (substeps-e (floor target-time dt-e)))
      (format t "CFL dt estimate: ~f~%" dt-e)
      ;; (format t "Estimated dt ~f~%" (cl-mpm:sim-dt *sim*))
      (format t "CFL step count estimate: ~D~%" substeps-e)
      (setf (cl-mpm:sim-dt *sim*) dt-e)
      (setf substeps substeps-e))
    (format t "Substeps ~D~%" substeps)
    (time (loop for steps from 0 to 100
                while *run-sim*
                do
                   (progn
                     ;; (if (> steps 3)
                     ;;     (setf (cl-mpm::sim-enable-damage *sim*) t)
                     ;;     (setf (cl-mpm::sim-enable-damage *sim*) nil)
                     ;;     )
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d_~2,'0d.vtk" *sim-step* rank)) *sim*)
                     ;(cl-mpm/output::save-vtk-nodes (merge-pathnames (format nil "output/sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     ;(cl-mpm/output:save-csv (merge-pathnames (format nil "output/sim_~5,'0d.csv" *sim-step*)) *sim*)
                     (let ((cfl 0))
                       (time (dotimes (i substeps)
                               (cl-mpm::update-sim *sim*)
                               (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))

                       (format t "CFL: ~f~%" cfl)
                       (push cfl *cfl-max*)
                       (let* ((dt-e (* dt-scale (cl-mpm::calculate-min-dt *sim*)))
                              (substeps-e (floor target-time dt-e)))
                         (format t "CFL dt estimate: ~f~%" dt-e)
                         (format t "CFL step count estimate: ~D~%" substeps-e)
                         (setf (cl-mpm:sim-dt *sim*) dt-e)
                         (setf substeps substeps-e))
                         )
                     (incf *sim-step*)
                     ;; (plot *sim*)
                     ;; (vgplot:print-plot (merge-pathnames (format nil "outframes/frame_~5,'0d.png" *sim-step*))
                     ;;                    :terminal "png size 1920,1080"
                     ;;                    )
                     (sleep .01)
                     ))))
  (cl-mpm/output:save-vtk (merge-pathnames (format nil "output/sim_~5,'0d.vtk" *sim-step*)) *sim*)
  )

(defmacro time-form (it form)
  `(progn
     (declaim (optimize speed))
     (let* ((iterations ,it)
            (start (get-internal-real-time)))
       (dotimes (i ,it)
         ,form)
       (let* ((end (get-internal-real-time))
              (units internal-time-units-per-second)
              (dt (/ (- end start) (* iterations units)))
              )
         (format t "Total time: ~f ~%" (/ (- end start) units)) (format t "Time per iteration: ~f~%" (/ (- end start) (* iterations units)))
         (format t "Throughput: ~f~%" (/ 1 dt))
         dt))))

(defun mpi-main-loop ()
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (setf lparallel:*kernel* (lparallel:make-kernel 16 :name "custom-kernel"))
    (format t "Test kernel~%")
    (time 
      (lparallel:pdotimes (i 100000)
                         (loop for j from 0 to 100000
                               do 
                               (sqrt i))))
    (format t "Setup~%")
    (setup 'cl-mpm/mpi::mpm-sim-mpi-stress)
    ;; (setf (cl-mpm::sim-dt *sim*) 1d0)
    (when (= rank 0)
      (cl-mpm/output:save-vtk-mesh (merge-pathnames "output/mesh.vtk") *sim*))
    (format t "Decompose~%")
    (cl-mpm/mpi::domain-decompose *sim*)

    (let* ((target-time 1d3)
           (dt-scale 1d0)
           (substeps (floor target-time (cl-mpm::sim-dt *sim*))))
      (cl-mpm::update-sim *sim*)
      (time (let* ((dt-e (* dt-scale (cl-mpm/mpi::calculate-min-dt *sim*)))
                  (substeps-e (floor target-time dt-e)))
              (format t "CFL dt estimate: ~f~%" dt-e)
              (format t "CFL step count estimate: ~D~%" substeps-e)
             (setf (cl-mpm:sim-dt *sim*) dt-e)
             (setf substeps substeps-e))
           )
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Run~%")

      (dotimes (step 100)
        (when (= rank 0)
          (format t "Update step ~D~%" step))
         ;(if (> step 3)
         ;    (progn (setf (cl-mpm::sim-enable-damage *sim*) t))
         ;    (progn (setf (cl-mpm::sim-enable-damage *sim*) nil)))

        (time (dotimes (i substeps)
                (cl-mpm::update-sim *sim*)))
        (let* ((dt-e (* dt-scale (cl-mpm/mpi::calculate-min-dt *sim*)))
               (substeps-e (floor target-time dt-e)))
          (format t "CFL dt estimate: ~f~%" dt-e)
          (format t "CFL step count estimate: ~D~%" substeps-e)
          (setf (cl-mpm:sim-dt *sim*) dt-e)
          (setf substeps substeps-e))
        (incf *sim-step*)
        (cl-mpm/output:save-vtk (merge-pathnames (format nil
                                                         "output/sim_rank_~2,'0d_~5,'0d.vtk" rank *sim-step* 
                                                         )) *sim*)
        )
      )
    (format t "rank: ~D Finished~%" rank))
  )

(declaim (notinline mpi-run))
(defun mpi-run (total-rank-count)
  (setf lparallel:*kernel* (lparallel:make-kernel 16 :name "custom-kernel"))
  (format t "Collecting servers~%")
  ;(when (> total-rank-count 0)
  (cl-mpm/mpi::collect-servers total-rank-count)
  (push (lambda () (print "Killing servers")
          (cl-mpm/mpi::kill-servers)
          )
        sb-ext:*exit-hooks*)
    ;)
  ;; (lfarm:broadcast-task (lambda ()
  ;;                         (progn
  ;;                           (load "test-file.lisp")
  ;;                           t)))
  (lfarm:broadcast-task (lambda ()
    (progn
      ;; (ql:quickload :cl-mpm/examples/slump)
      ;; (in-package :cl-mpm/examples/slump)
      (print *package*)
      (setf lparallel:*kernel* (lparallel:make-kernel 16))
      ;; (setup)
      (mpi-main-loop)
    t)))
  ;; (let ((n 10))
  ;;   (format t "Setup ~%")
  ;;   (setup 'cl-mpm/damage::mpm-sim-damage)
  ;;   (format t "Sim type ~A~%" (type-of *sim*))
  ;;   (cl-mpm::update-sim *sim*)
  ;;   (format t "Single core ~%")
  ;;   (time-form n
  ;;              (cl-mpm::update-sim *sim*))
  ;;   (setup 'cl-mpm/mpi::mpm-sim-mpi-stress)
  ;;   (format t "Sim type ~A~%" (type-of *sim*))
  ;;   (format t "MPI ~%")
  ;;   (cl-mpm::update-sim *sim*)
  ;;   (time-form n
  ;;              (cl-mpm::update-sim *sim*)))
  ;; (uiop:quit)
  (format t "Done ~%")
  (format t "Kill servers ~%")
  (cl-mpm/mpi::kill-servers)
  0
  )
 (mpi-main-loop)
