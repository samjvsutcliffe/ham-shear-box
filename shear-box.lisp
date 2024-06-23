(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/shear-box)
(ql:quickload :parse-float)


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
  (let* ((target-time 0.001d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (dt-scale 0.9d0)
         (load-steps 100)
         (displacment 8d-3)
         (disp-inc (/ displacment load-steps)))

    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (cl-mpm::update-sim *sim*)
    (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                    (format t "CFL dt estimate: ~f~%" dt-e)
                    (format t "CFL step count estimate: ~D~%" substeps-e)
                    (setf substeps substeps-e))
    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d0 (cl-mpm/setup::estimate-critical-damping *sim*)))

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

    (cl-mpm/dynamic-relaxation:converge-quasi-static
     *sim*
     :energy-crit 1d-2
     :oobf-crit 1d0
     :substeps 10
     :conv-steps 200)

    (loop for mp across (cl-mpm:sim-mps *sim*)
          do (when (= (cl-mpm/particle::mp-index mp) 0)
               (setf (cl-mpm/particle::mp-enable-plasticity mp) t)))

    (setf (cl-mpm:sim-damping-factor *sim*)
          (* 1d-1 (cl-mpm/setup::estimate-critical-damping *sim*))
          ;; (cl-mpm::sim-enable-damage *sim*) t
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
                         (format stream "~f,~f~%" disp-av load-av)))
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
         (output-dir (format nil "./output-~D-~E/" 0 0)))
    (format t "Refine: ~A~%" refine)
    (format t "Load: ~A~%" load)
    (ensure-directories-exist (merge-pathnames output-dir))
    (setup :refine refine :mps 4 :surcharge-load load)
    (run output-dir)))

(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "1"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
(lparallel:end-kernel)
