(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/shear-box)

(defun setup-test-column (size offset block-size &optional (e-scale 1) (mp-scale 1) &key (angle 0d0) (friction 0.1d0) (surcharge-load 72.5d3))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm::mpm-sim-usf
               ;:sim-type 'cl-mpm/damage::mpm-sim-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         ;; (floor-offset (* h-y 2))
         (floor-offset 2d0)
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((angle-rad (* angle (/ pi 180)))
             ;; (init-stress 60d3)
             (init-stress 50d3)
             (gf 48d0)
             (length-scale 1.5d-2)
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             )
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-chalk-brittle
                ;; 'cl-mpm/particle::particle-vm
                'cl-mpm/particle::particle-mc
                :E 1d9
                :nu 0.24d0
                :enable-plasticity t
                :psi 0d0
                :phi (* 42d0 (/ pi 180))
                :c 131d3
                :index 0
                :gravity 0.0d0
                ))))
      (let* ((sur-height (* 0.5 (second block-size)))
             (sur-size (list 0.06d0 sur-height))
             (load surcharge-load)
             (gravity (/ load (* density sur-height))))
        (format t "Gravity ~F~%" gravity)
        (cl-mpm::add-mps
         sim
         (cl-mpm/setup::make-mps-from-list
          (cl-mpm/setup::make-block-mps-list
           (mapcar #'+ offset (list 0d0 (second block-size)))
           sur-size
           (mapcar (lambda (e) (* e e-scale 2)) sur-size)
           density
           'cl-mpm/particle::particle-elastic-damage
           :E 1d9
           :nu 0.24d0
           :initiation-stress 1d20
           :index 1
           :gravity (- gravity))))
        )
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-enable-fbar sim) t)
      ;; (setf (cl-mpm::sim-mass-filter sim) 1d0)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (let ((ms 1d0))
        (setf (cl-mpm::sim-mass-scale sim) ms)
        (setf (cl-mpm:sim-damping-factor sim)
              (* 0.1d0 (cl-mpm/setup::estimate-critical-damping sim))))
      (let ((dt-scale 1d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))
      (format t "Estimated dt ~F~%" (cl-mpm:sim-dt sim))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))
      sim)))
;;Ad-hoc uncoupled strain-softening plasticity
(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
                   (c cl-mpm/particle::mp-c)
                   (phi cl-mpm/particle::mp-phi)
                   )
      mp
    (let ((phi_0 (* 42d0 (/ pi 180)))
          (phi_1 (* 30d0 (/ pi 180)))
          (c_0 131d3)
          (soft 10d0))
      (setf
       c (* c_0 (exp (- (* soft ps))))
       phi (+ phi_1 (* (- phi_0 phi_1) (exp (- (* soft ps)))))))))

(defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
  (with-accessors ((damping cl-mpm::sim-damping-factor)
                   (mass-scale cl-mpm::sim-mass-scale)
                   (mesh cl-mpm::sim-mesh)
                   (dt cl-mpm::sim-dt))
      sim
    (cl-mpm::iterate-over-nodes
     mesh
     (lambda (node)
       (cl-mpm::calculate-forces node damping dt mass-scale)
       ;(cl-mpm::calculate-forces-cundall-conservative node damping dt mass-scale)
       ))))

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
      (time
        (progn
          ,@body))
      (progn
        ,@body)))

(defun run (&optional (output-directory "./output/"))
  (cl-mpm/output:save-vtk-mesh (merge-pathnames output-directory "mesh.vtk")
                          *sim*)

  (defparameter *data-t* (list))
  (defparameter *data-disp* (list))
  (defparameter *data-v* (list))
  (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
    (format stream "disp,load~%"))
  (vgplot:close-all-plots)
  (let* ((dt (cl-mpm:sim-dt *sim*))
         (displacment 6d-3)
         (total-time (* displacment 1d1)) 
         (dt-scale 0.1d0)
         (load-steps 200)
         (target-time (/ displacment load-steps)) 
         (substeps (floor target-time dt))
         (disp-inc (/ displacment load-steps)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 10d0 (cl-mpm/setup::estimate-critical-damping *sim*))
          (cl-mpm::sim-enable-damage *sim*) nil)

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d-2
     :substeps 20
     :conv-steps 200)

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))

    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d0 
             (sqrt (cl-mpm:sim-mass-scale *sim*))
             (cl-mpm/setup::estimate-critical-damping *sim*))
          (cl-mpm::sim-enable-damage *sim*) t
          )
    (format t "Substeps ~D~%" substeps)
    (setf cl-mpm/penalty::*debug-force* 0)
    (time (loop for steps from 0 below load-steps
                while *run-sim*
                do
                   (progn
                     (format t "Step ~d ~%" steps)
                     (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~5,'0d.vtk" *sim-step*)) *sim*)
                     (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~5,'0d.vtk" *sim-step*)) *sim*)
                     (let ((load-av 0d0)
                           (disp-av 0d0))
                       (time
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          (incf load-av (/ cl-mpm/penalty::*debug-force* substeps))
                          (incf disp-av (/ *displacement-increment* substeps))
                          (incf *displacement-increment* (/ disp-inc substeps))
                          (incf *t* (cl-mpm::sim-dt *sim*))))

                       ;; (setf load-av cl-mpm/penalty::*debug-force*)
                       ;; (setf disp-av *displacement-increment*)

                       (push *t* *data-t*)
                       (push disp-av *data-disp*)
                       (push load-av *data-v*)
                       (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                         (format stream "~f,~f~%" 
                                 disp-av 
                                 ;(* load-av (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh *sim*)))
                                 load-av)))
                     (incf *sim-step*)
                     (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                       (format t "CFL dt estimate: ~f~%" dt-e)
                       (format t "CFL step count estimate: ~D~%" substeps-e)
                       (setf substeps substeps-e))
                     (swank.live:update-swank)))))
  )

(defun mpi-loop ()
  (let* ((refine (if (uiop:getenv "REFINE") (parse-integer (uiop:getenv "REFINE")) 2))
         (load (if (uiop:getenv "LOAD") (parse-float:parse-float (uiop:getenv "LOAD")) 72.5d3))
         (output-dir (format nil "./output-~D-~E/" refine load)))
    (format t "Refine: ~A~%" refine)
    (format t "Load: ~A~%" load)
    (ensure-directories-exist (merge-pathnames output-dir))
    (setup :refine refine :mps 4 :surcharge-load load :friction 0.0d0)
    ;(setup :refine refine :mps 2 :surcharge-load load :friction 0.0d0)
    (run output-dir)))

(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "1"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
(lparallel:end-kernel)
