(restrict-compiler-policy 'speed 3 3)
(restrict-compiler-policy 'debug 0 0)
(restrict-compiler-policy 'safety 0 0)
(setf *block-compile-default* t)

(in-package :cl-mpm/examples/shear-box)

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-chalk-brittle) dt)
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
        (when (< damage 1d0)
          ;(setf damage-increment (cl-mpm/damage::tensile-energy-norm strain E de))
          ;(setf damage-increment
          ;      (max 0d0
          ;           (cl-mpm/damage::drucker-prager-criterion
          ;            (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
		  (setf damage-increment
                (max 0d0
                     (cl-mpm/damage::criterion-dp
                      (magicl:scale stress (/ 1d0 (magicl:det def))) (* angle (/ pi 180d0)))))
          ;(incf damage-increment
          ;      (* E (cl-mpm/particle::mp-strain-plastic-vm mp)))
          )
        (when (>= damage 1d0)
          (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

(defun domain-decompose (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (let ((dsize (floor (cl-mpi:mpi-comm-size))))
      (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count sim) (list dsize 1 1)))
    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps sim)))
      (format t "Decompose~%"))
    (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
	    (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))
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
             ;; (init-stress 60d3)
             (init-stress 300d3)
             ;(init-stress 300d3)
             ;(gf 48d0)
             (gf 100d0)
             (length-scale (* 1 h))
             ;; (length-scale 0.015d0)
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             )
        (format t "Ductility ~E~%" ductility)
        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                ;; 'cl-mpm/particle::particle-chalk-brittle
                ;; 'cl-mpm/particle::particle-vm

                'cl-mpm/particle::particle-chalk-delayed
                :E 1d9
                :nu 0.24d0
                ;:nu 0.00d0
                :kt-res-ratio 1d-9
                :kc-res-ratio 1d0
                :g-res-ratio 1d-9
                ;:friction-angle 43d0
                :friction-angle 43d0
                :initiation-stress init-stress;18d3
                :delay-time 1d-3
                :delay-exponent 1d0
                :ductility ductility
                :local-length length-scale
                :local-length-damaged 10d-10

                :damage (coerce *damage* 'double-float)

                :enable-plasticity t
                :psi 0d0
                ;:phi (* 42d0 (/ pi 180))
                ;:c (* 131d3 10d0)
                :phi (* 50d0 (/ pi 180))
                :c (* 131d3 10d0)

                ;'cl-mpm/particle::particle-mc
                ;:E 1d9
                ;:nu 0.24d0
                ;:enable-plasticity t
                ;:psi 0d0
                ;:phi (* 42d0 (/ pi 180))
                ;:c 131d3
                ;:phi-r (* 30d0 (/ pi 180))
                ;:c-r 0d0
                ;:softening 10d0

                :index 0
                :gravity 0.0d0
                ))))
      ;; (let* ((sur-height (* 0.5 (second block-size)))
      ;;        (sur-size (list 0.06d0 sur-height))
      ;;        (ld surcharge-load)
      ;;        (gravity (/ ld (* density sur-height))))
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
      ;;      :E 1d9
      ;;      :nu 0.24d0
      ;;      :initiation-stress 1d20
      ;;      :local-length 0d0
      ;;      :index 1
      ;;      :gravity (- gravity))))
      ;;   )
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
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc::make-bc-fixed i '(nil nil 0)))
            ))
      (defparameter *pressure-bc*
        (cl-mpm/buoyancy::make-bc-pressure
         sim
         0d0
         (- surcharge-load)
         :clip-func
         (lambda (pos)
           (and
            (> (cl-mpm/utils:varef pos 1)
               (+ (second offset) (* 0.5d0 (second block-size))))
            (> (cl-mpm/utils:varef pos 0)
               (first offset))
            (< (cl-mpm/utils:varef pos 0)
               (+ (first offset) (first block-size)))
            )

           )))
      (cl-mpm:add-bcs-force-list
       sim
       *pressure-bc*)
      sim)))

(defun setup (&key (refine 1d0) (mps 2) (friction 0.0d0) (surcharge-load 72.5d3))
  (defparameter *displacement-increment* 0d0)
  (let* ((mps-per-dim mps)
         (mesh-size (/ 0.03d0 refine))
         (sunk-size 0.03d0)
         (box-size (* 2d0 sunk-size))
         (domain-size (* 3d0 box-size))
         (box-offset (* mesh-size 1d0))
         (offset (list box-size box-offset))
         (rank (cl-mpi:mpi-comm-rank))
         )
    (setf *box-size* box-size)
    (defparameter *sim* (setup-test-column
                         (list domain-size (+ (* 1.5 box-size) box-offset))
                         offset
                         (list box-size box-size)
                         (/ 1d0 mesh-size)
                         mps-per-dim
                         :friction friction
                         :surcharge-load surcharge-load))
    (make-penalty-box *sim* box-size (* 2d0 box-size) sunk-size friction box-offset)
    (domain-decompose *sim*)
    (defparameter *true-load-bc* *shear-box-left-dynamic*)
    ;; (setf (cl-mpm::sim-bcs-force-list *sim*)
    ;;       (list
    ;;        (cl-mpm/bc:make-bcs-from-list
    ;;         (list
    ;;          (cl-mpm/bc::make-bc-closure
    ;;           nil
    ;;           (lambda ()
    ;;             (apply-penalty-box box-size (* 2d0 box-size) sunk-size friction)))))))
    (when (= rank 0)
      (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*))))
    (defparameter *run-sim* t)
    (defparameter *t* 0)
    (defparameter *sim-step* 0)))

;;Ad-hoc uncoupled strain-softening plasticity
;(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt) 
;  (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
;                   (c cl-mpm/particle::mp-c)
;                   (phi cl-mpm/particle::mp-phi)
;                   )
;      mp
;      ;(setf phi (* 30d0 (/ pi 180))
;      ;      c 0d0)
;    (let ((phi_0 (* 42d0 (/ pi 180)))
;          (phi_1 (* 30d0 (/ pi 180)))
;          (c_0 131d3)
;          (soft 000d0))
;      (setf
;       c (* c_0 (exp (- (* soft ps))))
;       phi (+ phi_1 (* (- phi_0 phi_1) (exp (- (* soft ps)))))))))

(defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt)
  ;;Do sweet nothing 
  )

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
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (when (= rank 0)
      (format t "Output dir ~A~%" output-directory)
      (ensure-directories-exist (merge-pathnames output-directory))
      ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames output-directory "mesh.vtk") *sim*)
      (cl-mpm/output::save-simulation-parameters
       (merge-pathnames output-directory "settings.json")
       *sim*)
      (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :supersede)
        (format stream "disp,load~%")))
    (vgplot:close-all-plots)
    (let* ((displacment 6d-3)
           (total-time (* 10d0 displacment))
           (load-steps (floor (* 100 (/ displacment 1d-3))))
           (target-time (/ total-time load-steps))
           (dt (cl-mpm:sim-dt *sim*))
           (substeps (floor target-time dt))
           (dt-scale 0.25d0)
           (enable-plasticity (cl-mpm/particle::mp-enable-plasticity (aref (cl-mpm:sim-mps *sim*) 0)))
           (disp-inc (/ displacment load-steps)))

      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 0.1d0 
               (sqrt (cl-mpm:sim-mass-scale *sim*))
               (cl-mpm/setup::estimate-critical-damping *sim*))
            (cl-mpm::sim-enable-damage *sim*) nil)

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (= (cl-mpm/particle::mp-index mp) 0)
                 (setf (cl-mpm/particle::mp-enable-plasticity mp) nil)))

      (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))

      (cl-mpm/dynamic-relaxation:converge-quasi-static
       *sim*
       :energy-crit 1d-2
       :oobf-crit 1d-2
       :dt-scale 0.5d0
       :substeps 20
       :conv-steps 200)

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (= (cl-mpm/particle::mp-index mp) 0)
                 (setf (cl-mpm/particle::mp-enable-plasticity mp) enable-plasticity)))

      (setf (cl-mpm:sim-damping-factor *sim*)
            (* 1d-2
               (sqrt (cl-mpm:sim-mass-scale *sim*))
               (cl-mpm/setup::estimate-critical-damping *sim*))
            )

      (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
      (setf substeps (round target-time (cl-mpm:sim-dt *sim*)))

      (loop for mp across (cl-mpm:sim-mps *sim*)
            do (when (typep mp 'cl-mpm/particle::particle-damage) 
                 (when (= (cl-mpm/particle::mp-index mp) 0)
                   (setf (cl-mpm/particle::mp-delay-time mp) (* target-time 1d-2)))))

      (when (slot-exists-p *sim* 'cl-mpm/damage::delocal-counter-max)
        (setf (cl-mpm/damage::sim-damage-delocal-counter-max *sim*) substeps))

      (when (= rank 0)
        (format t "Substeps ~D~%" substeps))

      (cl-mpm:update-sim *sim*)
      (let ((disp-av 0d0)
            (load-av 0d0))
        ;; (push *t* *data-t*)
        ;; (push disp-av *data-disp*)
        ;; (push load-av *data-v*)
        (setf load-av cl-mpm/penalty::*debug-force*)
        (setf disp-av *displacement-increment*)
        (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
          (format stream "~f,~f~%" disp-av load-av)))
      (setf cl-mpm/penalty::*debug-force* 0
            (cl-mpm::sim-enable-damage *sim*) t)
      (time (loop for steps from 0 below load-steps
                  while *run-sim*
                  do
                     (progn
                       (when (= rank 0)
                         (format t "Step ~d ~%" steps))
                       (when (= (mod steps 0) 0)

                         (cl-mpm/output:save-vtk (merge-pathnames output-directory (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                         (cl-mpm/output::save-vtk-nodes (merge-pathnames output-directory (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                         )
                       (let ((load-av 0d0)
                             (disp-av 0d0))
                         (time
                          (dotimes (i substeps)
                            (cl-mpm::update-sim *sim*)
                            ;; (incf load-av (/ (get-load) substeps))
                            ;; (incf disp-av (/ *displacement-increment* substeps))
                            (incf *displacement-increment* (/ disp-inc substeps))
                            (incf *t* (cl-mpm::sim-dt *sim*))))

                         (setf load-av (cl-mpm/mpi:mpi-sum (get-load)))
                         (setf disp-av *displacement-increment*)

                         ;; (push *t* *data-t*)
                         ;; (push disp-av *data-disp*)
                         ;; (push load-av *data-v*)
                         (when (= rank 0)
                           (format t "Disp ~E - Load ~E~%" disp-av load-av)
                           (with-open-file (stream (merge-pathnames output-directory "disp.csv") :direction :output :if-exists :append)
                             (format stream "~f,~f~%" disp-av load-av))))
                       (incf *sim-step*)
                       (multiple-value-bind (dt-e substeps-e) (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                         (when (= rank 0)
                           (format t "CFL dt estimate: ~f~%" dt-e)
                           (format t "CFL step count estimate: ~D~%" substeps-e))
                         (setf substeps substeps-e))
                       (swank.live:update-swank)))))))

(defparameter *damage* 0d0)
(defun mpi-loop ()
  (let* ((refine (if (uiop:getenv "REFINE") (parse-integer (uiop:getenv "REFINE")) 2))
         (load (if (uiop:getenv "LOAD") (parse-float:parse-float (uiop:getenv "LOAD")) 72.5d3))
         (damage (if (uiop:getenv "DAMAGE") (parse-float:parse-float (uiop:getenv "DAMAGE")) 0d0))
         (output-dir (format nil "./output-~D-~f/" refine load)))
    (setf *damage* damage)
    (format t "Refine: ~A~%" refine)
    (format t "Load: ~A~%" load)
    (format t "Damage: ~A~%" damage)
    (ensure-directories-exist (merge-pathnames output-dir))
    (setup :refine refine :mps 4 :surcharge-load load :friction 0.0d0)
    ;; (let ((*output-directory* output-dir)
    ;;       (*sim-step* 0)
    ;;       (rank (cl-mpi:mpi-comm-rank)))
    ;;   (cl-mpm/output:save-vtk (merge-pathnames *output-directory* (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
    ;;   ;; (cl-mpm/output::save-vtk-nodes (merge-pathnames *output-directory* (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
    ;;   )
    (run output-dir)
    ))

(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "1"))))
  (setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
(lparallel:end-kernel)
