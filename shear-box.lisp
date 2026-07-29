;; (sb-ext::restrict-compiler-policy 'speed  0 0)
;; (sb-ext::restrict-compiler-policy 'debug  3 3)
;; (sb-ext::restrict-compiler-policy 'safety 3 3)
;; (setf *block-compile-default* t)
(sb-ext:restrict-compiler-policy 'speed 3 3)
(sb-ext:restrict-compiler-policy 'debug 0 0)
(sb-ext:restrict-compiler-policy 'safety 0 0)
;(setf *block-compile-default* t)

(in-package :cl-mpm/examples/shear-box)

(declaim (optimize (debug 0) (safety 0) (speed 3)))


(defparameter *damage-model* :MC)

(defun domain-decompose (sim)
  (when (typep sim 'cl-mpm/mpi::mpm-sim-mpi)
    (let ((rank (cl-mpi:mpi-comm-rank)))
      (let ((height 1))
        (when (> (cl-mpi:mpi-comm-size) 8)
          (setf height 2))
        (let* ((dsize (ceiling (cl-mpi:mpi-comm-size) height))
               (dsize-square (ceiling (sqrt (cl-mpi:mpi-comm-size)) 2)))
          (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count sim)
                                        ;(list dsize-square dsize-square 1)
                (list dsize height 1)
                )))
      (when (= rank 0)
        (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps sim)))
        (format t "Decompose~%"))
      (let ((mp-0 (aref (cl-mpm:sim-mps *sim*) 0)))
        (when (slot-exists-p mp-0 'cl-mpm/particle::local-length)
          (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length mp-0))))
	          (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))))
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
                               (mapcar (lambda (p)
                                         (when (and (> p 0d0)
                                                    (< p 1d0))
                                           (setf p (min 1d0 (+ (/ p 2) 0))))
                                         p) y)
                               z)))
               (format t "Domain ~A ~A~%" (list x y z) dnew)
               dnew)))))
      (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps sim))))))

(defparameter *overscale* (if (uiop:getenv "OVER") (parse-float:parse-float (uiop:getenv "OVER")) 0d0))


;; (defmpgen make-mps-plastic-damage
;;   'cl-mpm/particle::particle-chalk-brittle
;;   :E *elastic-constant*
;;   :nu 0.24d0
;;   :kt-res-ratio 1d0
;;   :kc-res-ratio 0d0
;;   :friction-angle (cl-mpm/utils:deg-to-rad angle)
;;   :residual-friction (cl-mpm/utils:deg-to-rad angle-r)
;;   :initiation-stress init-stress
;;   :ductility ductility
;;   :local-length length-scale
;;   :enable-damage t
;;   :enable-plasticity t
;;   :psi (cl-mpm/utils:deg-to-rad 5d0)
;;   :oversize 0d0;pd-inflection
;;   )

(defun setup-test-column (size offset block-size &optional (e-scale 1) (mp-scale 1)
                          &key
                            (angle 0d0)
                            (friction 0d0)
                            (surcharge-load 72.5d3)
                            (piston-scale 1d0)
                            (piston-mps 2)
                            (init-stress 50d3))
  (let* ((multigrid-refinement 1)
         (sim (cl-mpm/setup:make-simple-sim
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm/dynamic-relaxation::mpm-sim-dr-damage-ul
               :args-list (list :enable-aggregate t
                                :enable-fbar nil
                                :gravity 0d0
                                :max-split-depth 8)))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 1.7d3)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size)))
    (declare (double-float h density))
    (progn
      (let* ((E *elastic-constant*)
             (angle 42d0)
             (angle-rad (* angle (/ pi 180)))
             (angle-r 30d0)
             (gf *gf*)
             (length-scale (* 1d-3))
             (ductility
               (cl-mpm/damage::estimate-ductility-jirsek2004
                gf
                length-scale
                init-stress E))
             (pd-inflection 0d0)
             )
        (format t "Estimated ductility ~E~%" ductility)
        (format t "Init stress ~E~%" init-stress)
        (format t "PD inflection point ~E~%" pd-inflection)
        (make-mps-plastic-damage)) 
      (cl-mpm:iterate-over-mps 
        (cl-mpm:sim-mps sim)
        (lambda (mp)
          (cl-mpm/damage::set-mp-damage mp *damage*)))
      (let* ((sur-height h-x)
             (sur-size (list 0.06d0 sur-height))
             ;(load surcharge-load)
             (load (/ surcharge-load (- 1d0 *damage*)))
             )
        (cl-mpm:iterate-over-mps
         (cl-mpm:sim-mps sim)
         (lambda (mp)
           (with-accessors ((stress cl-mpm/particle:mp-stress)
                            (strain cl-mpm/particle:mp-strain)
                            (strain-n cl-mpm/particle:mp-strain-n)
                            (E cl-mpm/particle::mp-E)
                            (nu cl-mpm/particle::mp-nu)
                            (de cl-mpm/particle::mp-elastic-matrix))
               mp
             (let* (;(k-ratio 0d0)
                    (k-ratio (/ nu (- 1d0 nu)))
                    ;;k0
                    ;; (k-ratio (- 1d0 (sin (* 42d0 (/ pi 180)))));;k0
                    (stresses (cl-mpm/utils:voigt-from-list (list
                                                             (- (* surcharge-load k-ratio))
                                                             (- surcharge-load)
                                                             (- (* surcharge-load k-ratio))
                                                             0d0
                                                             0d0
                                                             0d0)))
                    (strains (magicl:linear-solve de stresses)))
               (setf stress stresses
                     strain   strains
                     strain-n (cl-mpm/utils:voigt-copy strains)))))))
      (defparameter *mesh-resolution* h-x)
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-velocity-algorithm sim) :QUASI-STATIC)
      (when (typep *sim* 'cl-mpm/damage::mpm-sim-damage)
        (setf (cl-mpm::sim-nonlocal-damage sim) t)
        (setf (cl-mpm/damage::sim-enable-length-localisation sim) t))
      (cl-mpm/setup::set-mass-filter sim density :proportion 1d-15)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (cl-mpm/setup::setup-bcs
       sim :left '(0 nil nil))
      sim)))




(defun setup (&key
                (refine 1d0)
                (mps 4)
                (friction 0.0d0)
                (surcharge-load 72.5d3)
                (epsilon-scale 1d2)
                (piston-scale 1d0)
                (piston-mps 0)
                (init-stress 90d3)
                (mp-refine 2))
  (defparameter *displacement-increment* 0d0)
  (let* ((mps-per-dim mps)
         (mesh-size (/ 0.03d0 refine))
         (sunk-size 0.03d0)
         (box-size (* 2d0 sunk-size))
         (domain-size (* 3d0 box-size))
         (box-offset box-size)
         (rank (cl-mpi:mpi-comm-rank))
         (offset (list box-size box-offset)))
    (setf *box-size* box-size)
    (defparameter *sim* (setup-test-column
                         (list domain-size (+ (* 2 box-size) box-offset))
                         offset
                         (list box-size box-size)
                         (/ 1d0 mesh-size)
                         mps-per-dim
                         :piston-scale piston-scale
                         :piston-mps piston-mps
                         :surcharge-load surcharge-load
                         :init-stress init-stress))
    (make-penalty-box *sim* box-size (* 2d0 box-size) sunk-size friction box-offset
                      :epsilon-scale epsilon-scale
                      :corner-size (* 0.25d0 mesh-size)
                      :smoothness 2)
    (make-piston box-size box-offset surcharge-load epsilon-scale piston-scale)
    (dotimes (i mp-refine)
      (dolist (dir (list :y))
        (cl-mpm::split-mps-criteria
         *sim*
         (lambda (mp h)
           (when
               (and
                (>= (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 1)
                    (+ box-offset sunk-size (- mesh-size)))
                (<= (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 1)
                    (+ box-offset sunk-size mesh-size))
                (or
                 (<= (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 0)
                     (+ box-size (* 0.5d0 mesh-size)))
                 (>= (cl-mpm/utils:varef (cl-mpm/particle:mp-position mp) 0)
                     (- (* 2 box-size) (* 0.5d0 mesh-size))))
                )
             dir)))))
    (domain-decompose *sim*)
    (defparameter *true-load-bc* *shear-box-left-dynamic*)
    (when (= rank 0)
      (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
      (format t "Mesh-size: ~E~%" (cl-mpm/mesh::mesh-resolution (cl-mpm:sim-mesh *sim*)))))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *sim-step* 0))

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
      (time
        (progn
          ,@body))
      (progn
        ,@body)))

(defun get-load ()
  (let ((normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))))
    (cl-mpm/mpi:mpi-sum
     (+
      (cl-mpm/penalty::resolve-load-direction *shear-box-struct-left* normal)
      (cl-mpm/penalty::resolve-load-direction *shear-box-struct-right* normal)
      ))))

(defun get-load-left ()
  (let ((normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))))
    (cl-mpm/mpi:mpi-sum
      (cl-mpm/penalty::resolve-load-direction *shear-box-struct-left* normal))))

(defun get-load-right ()
  (let ((normal (cl-mpm/utils:vector-from-list (list 1d0 0d0 0d0))))
    (cl-mpm/mpi:mpi-sum 
      (cl-mpm/penalty::resolve-load-direction *shear-box-struct-right* normal))))



(defparameter *data-disp* nil)
(defparameter *data-v* nil)
(defparameter *data-damage* nil)
(defun run-adaptive
    (&key (output-dir "./output/") 
       (refine 1)
       (displacement 0.1d-3)
       (load-steps 10)
       (dt-scale 0.5d0)
       (enable-plasticity t)
       (enable-damage nil)
       (surcharge-load 0d0)
       )

  (let (;; (total-disp 0.12d-3)
        )
    (setf *displacement-increment* 0d0)
    (defparameter *data-disp* nil)
    (defparameter *data-v* nil)
    (defparameter *data-damage* nil)
    (cl-mpm/dynamic-relaxation::run-adaptive-load-control
     *sim*
     :output-dir output-dir
     :plotter (lambda (sim))
     :load-steps load-steps
     :substeps 50
     :criteria 1d-3
     :enable-damage enable-damage
     :enable-plastic enable-plasticity
     :max-adaptive-steps 0
     :min-adaptive-steps 0
     :max-damage-inc 0.5d0
     :max-plastic-inc nil;1d4
     :dt-scale 0.9d0
     :save-vtk-dr nil
     :save-vtk-loadstep t
     :post-iter-step (lambda (i o e)
                       (format t "Penalty load ~E - aim ~E~%"
                               (/ (get-piston-load *sim*) 0.06d0)
                               surcharge-load))
     :post-conv-step
     (lambda (sim)
       (save-disp sim output-dir)
       (push *displacement-increment* *data-disp*)
       (push (get-load) *data-v*)
       (push (cl-mpm/dynamic-relaxation::get-damage *sim*) *data-damage*))
     :loading-function (lambda (percent)
                         (setf *displacement-increment* (* displacement percent))))))


(defparameter *damage* 0d0)
(defun mpi-loop ()
  (let* ((refine (if (uiop:getenv "REFINE") (parse-integer (uiop:getenv "REFINE")) 2))
         (load (float (if (uiop:getenv "LOAD") (parse-float:parse-float (uiop:getenv "LOAD")) 72.5d3) 0d0))
         (damage (if (uiop:getenv "DAMAGE") (parse-float:parse-float (uiop:getenv "DAMAGE")) 0d0))
         (overscale (if (uiop:getenv "OVER") (parse-float:parse-float (uiop:getenv "OVER")) 0d0))
         (mps 4)
         (scale 1d0)
         (epsilon-scale (* 1d3 (- 1d0 damage)))
         ;(epsilon-scale 1d3)
         (piston-scale 1d0)
         ;(output-dir (format nil "./data/output-~F_0.5_~D_~f_~f_~F-~f/" refine mps scale damage overscale load))
         (output-dir (format nil "/nobackup/rmvn14/paper-1/plastic-damage-residual/output-~F_0.5_~D_~f_~f_~F-~f/" refine mps scale damage overscale load))
         )
    (setf *damage* damage)
    (format t "Refine: ~A~%" refine)
    (format t "Load: ~A~%" load)
    (format t "Damage: ~A~%" damage)
    (ensure-directories-exist (merge-pathnames output-dir))
    (setup
      :refine refine
      :mps mps
      :surcharge-load load
      :friction 0d0
      :epsilon-scale epsilon-scale
      :piston-scale piston-scale
      :init-stress
      (cl-mpm/damage::mohr-coloumb-coheasion-to-tensile
       131d3
       42d0))
    (cl-mpm::domain-sort-mps *sim*)


    (push (list :SCALAR "damage-tcs-c" #'cl-mpm/particle::mp-damage-compression) (cl-mpm::sim-output-list *sim*))
    (push (list :SCALAR "damage-tcs-s" #'cl-mpm/particle::mp-damage-shear) (cl-mpm::sim-output-list *sim*))
    (push (list :SCALAR "damage-tcs-t" #'cl-mpm/particle::mp-damage-tension) (cl-mpm::sim-output-list *sim*))
    ;(push (list :VOIGT "sig_u" (lambda (mp) (cl-mpm/fastmaths:fast-scale
    ;                                         (cl-mpm/particle::mp-undamaged-stress mp)
    ;                                         (/ 1d0 (cl-mpm/particle::mp-deformation-jacobian-strain mp))))) (cl-mpm::sim-output-list *sim*))
    (run-adaptive :output-dir output-dir
                  :displacement 3d-3
                  :load-steps 20
                  :refine refine
                  :enable-plasticity t
                  :enable-damage nil
                  :surcharge-load load)))

(let ((threads (parse-integer (if (uiop:getenv "OMP_NUM_THREADS") (uiop:getenv "OMP_NUM_THREADS") "1"))))
  ;(setf lparallel:*kernel* (lparallel:make-kernel threads :name "custom-kernel"))
  (cl-mpm/utils:set-workers threads)
  (format t "Thread count ~D~%" threads))
(defparameter *run-sim* nil)
(mpi-loop)
(lparallel:end-kernel)
