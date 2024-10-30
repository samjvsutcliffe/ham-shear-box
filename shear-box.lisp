(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/shear-box)

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-delayed) dt)
  (let ((damage-increment 0d0))
    (with-accessors ((stress cl-mpm/particle::mp-undamaged-stress)
                     (strain cl-mpm/particle::mp-strain)
                     (damage cl-mpm/particle:mp-damage)
                     (init-stress cl-mpm/particle::mp-initiation-stress)
                     (critical-damage cl-mpm/particle::mp-critical-damage)
                     (damage-rate cl-mpm/particle::mp-damage-rate)
                     (pressure cl-mpm/particle::mp-pressure)
                     (ybar cl-mpm/particle::mp-damage-ybar)
                     (def cl-mpm/particle::mp-deformation-gradient)
                     (angle cl-mpm/particle::mp-friction-angle)
                     (c cl-mpm/particle::mp-coheasion)
                     (nu cl-mpm/particle::mp-nu)
                     (ft cl-mpm/particle::mp-ft)
                     (fc cl-mpm/particle::mp-fc)
                     (E cl-mpm/particle::mp-e)
                     (de cl-mpm/particle::mp-elastic-matrix)
                     (kc-r cl-mpm/particle::mp-k-compressive-residual-ratio)
                     (kt-r cl-mpm/particle::mp-k-tensile-residual-ratio)
                     (g-r cl-mpm/particle::mp-shear-residual-ratio)
                     ) mp
      (declare (double-float pressure damage))
      (progn
        (setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
        ;(setf damage-increment
        ;      (max 0d0
        ;           (cl-mpm/damage::drucker-prager-criterion
        ;            (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
		;(setf damage-increment
        ;      (max 0d0
        ;           (cl-mpm/damage::criterion-dp
        ;            stress
        ;            ;(magicl:scale stress (/ 1d0 (magicl:det def))) 
        ;            (* angle (/ pi 180d0)))))
        ;(setf damage-increment (cl-mpm/damage::criterion-max-principal-stress stress))

        ;(incf damage-increment
        ;      (* E (cl-mpm/particle::mp-strain-plastic-vm mp)))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

(defun domain-decompose (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (let ((dsize (floor (cl-mpi:mpi-comm-size)))
          (dsize-square (floor (sqrt (cl-mpi:mpi-comm-size)))))
      (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count sim)
            (list dsize-square dsize-square 1)
            ;(list dsize 1 1)
            ))
    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps sim)))
      (format t "Decompose~%"))
    (let ((mp-0 (aref (cl-mpm:sim-mps *sim*) 0)))
      (when (slot-exists-p mp-0 'cl-mpm/particle::local-length)
        (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length mp-0))))
	        (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size)) ) )
    (let ((size 0.06d0))
      (cl-mpm/mpi::domain-decompose
       sim
       :domain-scaler
       (lambda (domain)
         (destructuring-bind (x y z) domain
           (let ((dnew (list (mapcar (lambda (p)
                                       (when (and (> p 0d0)
                                                  (< p 1d0))
                                         (setf p (min 1d0 (+ (/ p 3) 1/3))))
                                       p) x)
                             y z)))
             (format t "Domain ~A ~A~%" (list x y z) dnew)
             dnew)))))
    (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps sim)))))


(defun setup-test-column (size offset block-size &optional (e-scale 1) (mp-scale 1) &key (angle 0d0) (friction 0.1d0) (surcharge-load 72.5d3))
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               ;:sim-type 'cl-mpm::mpm-sim-usf
               ;:sim-type 'cl-mpm/damage::mpm-sim-damage
               ;:sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes
               :sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes-damage
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
             (init-stress 131d3)
             (gf 5d0)
             ;(length-scale (* 1 h))
             (length-scale (* 1 h))
             ;(length-scale (* 0.015d0 2))
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress 1d9)))
        (format t "Ductility ~E~%" ductility)
        (make-mps-damage))
      ;; (let* ((sur-height h-x)
      ;;        (sur-height (* 0.5 (second block-size)))
      ;;        (sur-size (list 0.06d0 sur-height))
      ;;        (ld surcharge-load)
      ;;        (gravity (if (> ld 0d0) (/ ld (* density sur-height)) 0d0))
      ;;        ;(gravity (/ ld (* density sur-height)))
      ;;        ;(gravity 0d0)
      ;;        )
      ;;   (format t "Gravity ~F~%" gravity)
      ;;   (cl-mpm::add-mps
      ;;    sim
      ;;    (cl-mpm/setup::make-mps-from-list
      ;;     (cl-mpm/setup::make-block-mps-list
      ;;      (mapcar #'+ offset (list 0d0 (second block-size)))
      ;;      sur-size
      ;;      (mapcar (lambda (e) (* e e-scale 2)) sur-size)
      ;;      density
      ;;      'cl-mpm/particle::particle-elastic-damage
      ;;      :E (* 1d9 1d0)
      ;;      :nu 0.24d0
      ;;      :initiation-stress 1d20
      ;;      :local-length 0d0
      ;;      :index 1
      ;;      :gravity (- gravity)
      ;;      ))))

      (defparameter *mesh-resolution* h-x)
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
					(cl-mpm/bc::make-outside-bc-varfix
					 (cl-mpm:sim-mesh sim)
					 '(0 nil nil)
					 '(0 nil nil)
					 '(nil 0 nil)
					 '(nil 0 nil)
					 '(nil nil 0)
					 '(nil nil 0)))
      sim)))

(defun setup (&key (refine 1d0) (mps 4) (friction 0.0d0) (surcharge-load 72.5d3)
                (epsilon-scale 1d2)
                (piston-scale 1d0))
  (defparameter *displacement-increment* 0d0)
  (let* ((mps-per-dim mps)
         (mesh-size (/ 0.03d0 refine))
         (sunk-size 0.03d0)
         (box-size (* 2d0 sunk-size))
         (domain-size (* 3d0 box-size))
         (box-offset (* mesh-size 2d0))
         (offset (list box-size box-offset))
         (rank (cl-mpi:mpi-comm-rank))
         )
    (setf *box-size* box-size)
    (defparameter *sim* (setup-test-column
                         (list domain-size (+ (* 2 box-size) box-offset))
                         offset
                         (list box-size box-size)
                         (/ 1d0 mesh-size)
                         mps-per-dim
                         :friction friction
                         :surcharge-load surcharge-load))
    (make-penalty-box *sim* box-size (* 2d0 box-size) sunk-size friction box-offset
                      :epsilon-scale epsilon-scale
                      :corner-size (* mesh-size 0.25d0))
    (make-piston box-size box-offset surcharge-load epsilon-scale piston-scale)
    (domain-decompose *sim*)
    (defparameter *true-load-bc* *shear-box-left-dynamic*)
    (when (= rank 0)
      (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (defparameter *run-sim* t)
    (defparameter *t* 0)
    (defparameter *sim-step* 0)))

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
      (time
        (progn
          ,@body))
      (progn
        ,@body)))

(defun get-load ()
  (cl-mpm/mpi:mpi-sum
   (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*)))


;(defun get-load ()
;    ;(cl-mpm/penalty::bc-penalty-load *true-load-bc*)
;      (- 
;        (cl-mpm/penalty::bc-penalty-load *shear-box-left-dynamic*)
;        (cl-mpm/penalty::bc-penalty-load *shear-box-right-dynamic*)))

(defun run (&key (output-directory "./output/") 
              (refine 1)
              (displacement 0.1d-3)
              (time-scale 1d0)
              (dt-scale 0.25d0)
              (damage-time-scale 1d0)
              (sample-scale 1d0)
              )
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (when (= rank 0)
      (format t "Output dir ~A~%" output-directory)
      (ensure-directories-exist (merge-pathnames output-directory))
      ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames output-directory "mesh.vtk") *sim*)
      (cl-mpm/output::save-simulation-parameters
       (merge-pathnames output-directory "settings.json")
       *sim*)
      (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
        (format stream "disp,load,plastic,damage,energy~%")))
    (vgplot:close-all-plots)
    (let* ((displacment displacement)
           (time-per-mm (* 100d0 time-scale))
           (total-time (* time-per-mm displacment))
           (load-steps (round (* sample-scale 500 (/ displacment 1d-3))))
           (target-time (/ total-time load-steps))
           (dt (cl-mpm:sim-dt *sim*))
           (substeps (floor target-time dt))
           (dt-scale 0.5d0)
           ;(enable-plasticity (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
           (enable-plasticity nil)
           (enable-damage t)
           (disp-inc (/ displacment load-steps)))

      (when (= rank 0)
        (format t "Plasticity: ~A~%Damage: ~A~%" enable-plasticity enable-damage))

      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.05d0
               (sqrt (cl-mpm:sim-mass-scale *sim*))
               (cl-mpm/setup::estimate-critical-damping *sim*))
            (cl-mpm::sim-enable-damage *sim*) nil)

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (= (cl-mpm/particle::mp-index mp) 0)
                 (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

      (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
      (time (cl-mpm:update-sim *sim*))
      (time (cl-mpm:update-sim *sim*))
      (time (cl-mpm:update-sim *sim*))


      (setf *enable-box-friction* nil)
      (cl-mpm/dynamic-relaxation:converge-quasi-static
       *sim*
       :energy-crit 1d-2
       :oobf-crit 1d-2
       :dt-scale dt-scale
       :substeps 10
       :conv-steps 400
       :post-iter-step
       (lambda (i e o)
         (when (= rank 0)
           (format t "Surcharge load ~E~%" (/ *piston-confinement* 10))
           (setf *piston-confinement* 0d0)))
       )

      (setf *enable-box-friction* t)

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (= (cl-mpm/particle::mp-index mp) 0)
                 (setf (cl-mpm/particle::mp-enable-plasticity mp) 
                       enable-plasticity))) 

      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 1d-2
               (cl-mpm/setup::estimate-critical-damping *sim*)))

      (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
      (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (typep mp 'cl-mpm/particle::particle-damage) 
                 (when (= (cl-mpm/particle::mp-index mp) 0)
                   (setf (cl-mpm/particle::mp-delay-time mp) (* target-time 1d-1)))))

      (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
        (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) substeps))

      (when (= rank 0)
        (format t "Substeps ~D~%" substeps))

      (cl-mpm:update-sim *sim*)
      (let ((disp-av 0d0)
            (load-av 0d0)
            (p-av 0d0)
            (d-av 0d0)
            (e-av 0d0)
            )
        ;; (push *t* *data-t*)
        ;; (push disp-av *data-disp*)
        ;; (push load-av *data-v*)
        (setf load-av (get-load))
        (setf disp-av *displacement-increment*)
        (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
          (format stream "~f,~f,~f,~f~%" disp-av load-av p-av d-av e-av)))
      (setf cl-mpm/penalty::*debug-force* 0
            (cl-mpm::sim-enable-damage *sim*) enable-damage)
      (time (loop for steps from 0 below load-steps
                  while *run-sim*
                  do
                     (progn
                       (when (= rank 0)
                         (format t "Step ~d ~%" steps))
                       (when (= (mod steps 10) 0)
                         (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                        (when (= rank 0)
                          (save-json-penalty-box (merge-pathnames output-directory (format nil "sim_pb_~5,'0d.json" *sim-step*)) *sim*) )
                         ;(cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                         )
                       (let ((load-av 0d0)
                             (disp-av 0d0)
                             (p-av 0d0)
                             (d-av 0d0)
                             (e-av 0d0))
                         (time
                          (dotimes (i substeps)
                            (cl-mpm::update-sim *sim*)
                            (incf load-av (/ (get-load) substeps))
                            (incf disp-av (/ *displacement-increment* substeps))
                            (incf *displacement-increment* (/ disp-inc substeps))
                            (incf e-av (/ (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*) substeps))
                            (incf *t* (cl-mpm::sim-dt *sim*))))
                         (setf d-av (get-damage))
                         (setf p-av (get-plastic))
                         (when (= rank 0)
                           (format t "Surcharge load ~E~%" (/ *piston-confinement* 10))
                           (setf *piston-confinement* 0d0))
                         (when (= rank 0)
                           (format t "Disp ~E - Load ~E~%" disp-av load-av)
                           (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                             (format stream "~f,~f,~f,~f,~f~%" disp-av load-av p-av d-av e-av))))

                       (incf *sim-step*)
                       (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                         (when (= rank 0)
                           (format t "CFL dt estimate: ~f~%" dt-e)
                           (format t "CFL step count estimate: ~D~%" substeps-e))
                         (setf substeps substeps-e))
                       (swank.live:update-swank)))))))
(defun run-static (&optional (output-directory "./output/"))
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (when (= rank 0)
      (format t "Output dir ~A~%" output-directory)
      (ensure-directories-exist (merge-pathnames output-directory))
      (cl-mpm/output::save-simulation-parameters
       (merge-pathnames output-directory "settings.json")
       *sim*)
      (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
        (format stream "disp,load~%")))

    (defparameter *data-t* (list))
    (defparameter *data-disp* (list))
    (defparameter *data-v* (list))
    
    (vgplot:close-all-plots)
    (let* ( (dt (cl-mpm:sim-dt *sim*))
           (dt-scale 0.5d0)
           (displacment 1d-3)
           (load-steps (* 500 (/ displacment  1d-3)))
           (enable-plasticity t)
           (disp-inc (/ displacment load-steps)))
      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (change-class mp 'cl-mpm/particle::particle-chalk-brittle))

      (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.05d0
               ;; (sqrt (cl-mpm:sim-mass-scale *sim*))
               (cl-mpm/setup::estimate-critical-damping *sim*)))

       (cl-mpm/dynamic-relaxation:converge-quasi-static
                          *sim*
                          :energy-crit 1d-2
                          :oobf-crit 1d-2
                          :substeps 10
                          :conv-steps 200
                          :post-iter-step
                          (lambda (i energy oobf)))

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (= (cl-mpm/particle::mp-index mp) 0)
                 (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))
      (setf *enable-box-friction* t)
      (defparameter *displacement-increment* 0d0)
      (vgplot:figure)
      (setf cl-mpm/penalty::*debug-force* 0)
      (time (loop for steps from 0 below load-steps
                  while *run-sim*
                  do
                     (progn
                       (when (= rank 0)
                         (format t "Step ~d ~%" steps))
                       (when (= (mod steps 10) 0)
                         (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                         (when (= rank 0)
                           (save-vtk-penalty-box (merge-pathnames output-directory (format nil "sim_box_~5,'0d.vtk" *sim-step*)) *sim*))
                                        ;(cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                         )
                       (let ((load-av 0d0)
                             (disp-av 0d0))
                         (incf *displacement-increment* disp-inc)
                         (cl-mpm/dynamic-relaxation:converge-quasi-static
                          *sim*
                          :energy-crit 1d-2
                          :oobf-crit 1d-2
                          :substeps 10
                          :conv-steps 200
                          :post-iter-step
                          (lambda (i energy oobf)))
                         (cl-mpm/damage::calculate-damage *sim*)
                         (setf load-av (cl-mpm/mpi:mpi-sum (get-load)))
                         (setf disp-av *displacement-increment*)
                         (when (= rank 0)
                           (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                             (format stream "~f,~f~%" disp-av load-av))))
                       (incf *sim-step*)
                       )))
      ))
  )

(defparameter *damage* 0d0)
(defun mpi-loop ()
  (let* ((refine (if (uiop:getenv "REFINE") (parse-integer (uiop:getenv "REFINE")) 2))
         (load (if (uiop:getenv "LOAD") (parse-float:parse-float (uiop:getenv "LOAD")) 72.5d3))
         (damage (if (uiop:getenv "DAMAGE") (parse-float:parse-float (uiop:getenv "DAMAGE")) 0d0))
         (output-dir (format nil "./output-~F-~f/" refine load)))
    (setf *damage* damage)
    (format t "Refine: ~A~%" refine)
    (format t "Load: ~A~%" load)
    (format t "Damage: ~A~%" damage)
    (ensure-directories-exist (merge-pathnames output-dir))
    (setup :refine refine :mps 3 :surcharge-load load)
    ;; (let ((*output-directory* output-dir)
    ;;       (*sim-step* 0)
    ;;       (rank (cl-mpi:mpi-comm-rank)))
    ;;   (cl-mpm/output:save-vtk (merge-pathnames *output-directory* (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
    ;;   ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames *output-directory* (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
    ;;   )
    (run :output-directory output-dir 
         :refine refine
         :scale 1d0
         :sample-scale 1d0)
    ;(run-static output-dir)
    ))

(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "1"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
(lparallel:end-kernel)
