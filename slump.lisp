;; (sb-ext:restrict-compiler-policy 'speed 3 3)
;; (sb-ext:restrict-compiler-policy 'debug 0 0)
;; (sb-ext:restrict-compiler-policy 'safety 0 0)
(in-package :cl-mpm/examples/slump)
;; (ql:quickload :cl-mpm/buoyancy)
;; (ql:quickload :cl-mpm/mpi)
;; (ql:quickload :cl-mpm/output)

(defmethod cl-mpm::update-stress-mp (mesh (mp cl-mpm/particle::particle-chalk-delayed) dt fbar)
  (cl-mpm::update-stress-kirchoff-damaged mesh mp dt fbar)
  ;(cl-mpm::update-stress-kirchoff mesh mp dt fbar)
  )

(defun cl-mpm/damage::length-localisation (local-length local-length-damaged damage)
  ;; (+ (* local-length (- 1d0 damage)) (* local-length-damaged damage))
  ;(* local-length (max (sqrt (- 1d0 damage)) 1d-10))
  local-length
  )

(defmethod cl-mpm/damage::damage-model-calculate-y ((mp cl-mpm/particle::particle-visco-elasto-plastic-damage) dt)
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
          (setf damage-increment (cl-mpm/damage::tensile-energy-norm-pressure strain E nu de
                                                                              (* pressure damage))))
        (when (>= damage 1d0)
          (setf damage-increment 0d0))
        ;;Delocalisation switch
        (setf (cl-mpm/particle::mp-damage-y-local mp) damage-increment)
        (setf (cl-mpm/particle::mp-local-damage-increment mp) damage-increment)
        ))))

(defun setup-test-column (size block-size offset &optional (e-scale 1) (mp-scale 1))
  (format t "Setup-test-column ~%")
  (format t "~F ~F ~A mp count"
          e-scale
          mp-scale
          block-size
          (mapcar (lambda (e) (floor (* e e-scale mp-scale))) block-size)
          )
  (let* ((sim (cl-mpm/setup::make-block
               (/ 1d0 e-scale)
               (mapcar (lambda (x) (* x e-scale)) size)
               :sim-type 'cl-mpm/mpi::mpm-sim-mpi-nodes-damage
               ))
         (h (cl-mpm/mesh:mesh-resolution (cl-mpm:sim-mesh sim)))
         (h-x (/ h 1d0))
         (h-y (/ h 1d0))
         (density 900d0)
         (angle 0d0)
         (elements (mapcar (lambda (s) (* e-scale (/ s 2))) size))
         )
    (declare (double-float h density))
    (progn
      (let* ((length-scale h)
             (init-stress 0.3d6)
             ;; (gf 100d0)
             ;; (ductility (estimate-ductility-jirsek2004 gf length-scale init-stress 1d9))
             (ductility (cl-mpm/damage::estimate-ductility-jirsek2004 5000d0 length-scale init-stress 1d9))
             ;; (ductility 8d0)
             )
        (format t "~F ~F ~A mp count"
                e-scale
                mp-scale
                block-size
                (mapcar (lambda (e) (floor (* e e-scale mp-scale))) block-size)
                )
		    (when (= (cl-mpi:mpi-comm-rank) 0) 
			    (format t "Estimated ductility ~E~%" ductility)
			    (format t "Estimated gf ~E~%" (cl-mpm/damage::gf-from-ductility ductility length-scale init-stress 1d9))
			    (when (< ductility 1d0)
			      (error "Ductility too low ~A" ductility)))

        (setf (cl-mpm:sim-mps sim)
              (cl-mpm/setup::make-mps-from-list
               (cl-mpm/setup::make-block-mps-list
                offset
                block-size
                (mapcar (lambda (e) (* e e-scale mp-scale)) block-size)
                density
                'cl-mpm/particle::particle-visco-elasto-plastic-damage
                :visc-factor 1.11d6
                ;; :visc-factor 11d6
                :visc-power 3d0
                :E 1d9
                :nu 0.325d0
                :enable-plasticity t
                :friction-angle 80.0d0

                :kt-res-ratio 1d-15
                :kc-res-ratio 1d-2
                :g-res-ratio 5d-3

                :fracture-energy 3000d0
                :initiation-stress init-stress

                :delay-time 1d1
                :delay-exponent 2d0
                :ductility ductility
                :critical-damage 1d0;(- 1.0d0 1d-3)
                :damage-domain-rate 0.9d0;This slider changes how GIMP update turns to uGIMP under damage
                :local-length length-scale;(* 0.20d0 (sqrt 7))
                :local-length-damaged 10d-10

                :psi (* 00d0 (/ pi 180))
                :phi (* 40d0 (/ pi 180))
                :c 1000d3

                :gravity -9.8d0
                ;; :gravity-axis (magicl:from-list '(0.5d0 0.5d0) '(2 1))
                :index 0
                ;:angle angle
                ;; :slope slope
                ))))
      (let* ((mp-0 (aref (cl-mpm:sim-mps sim) 0))
             (fc (cl-mpm/particle::mp-fc mp-0))
             (ft (cl-mpm/particle::mp-ft mp-0))
             (angle-d (* (/ 180 pi) (atan (* 3 (/ (- fc ft) (+ fc ft))))))
             (rc (cl-mpm/particle::mp-k-compressive-residual-ratio mp-0))
             (rs (cl-mpm/particle::mp-shear-residual-ratio mp-0))
             (angle-plastic (cl-mpm/particle::mp-phi mp-0))
             (angle-plastic-damaged (atan (* (/ rs rc) (tan angle-plastic))))
             )
        (when (= (cl-mpi:mpi-comm-rank) 0) 
          (format t "Chalk damage growth angle: ~F~%"
                  angle-d)
          (format t "Chalk plastic virgin angle: ~F~%"
                  (* (/ 180 pi) angle-plastic))
          (format t "Chalk plastic residual angle: ~F~%"
                  (* (/ 180 pi) angle-plastic-damaged))))
      (setf (cl-mpm:sim-allow-mp-split sim) nil)
      (setf (cl-mpm::sim-enable-damage sim) nil)
      (setf (cl-mpm::sim-nonlocal-damage sim) t)
      (setf (cl-mpm::sim-enable-fbar sim) nil)
      (setf (cl-mpm::sim-allow-mp-damage-removal sim) nil)
      (setf (cl-mpm::sim-mp-damage-removal-instant sim) nil)
      (setf (cl-mpm::sim-mass-filter sim) 1d-10)
      (let ((mass-scale 1d4))
        (setf (cl-mpm::sim-mass-scale sim) mass-scale)
        (setf (cl-mpm:sim-damping-factor sim)
              ;; 0.7d0
              (* 100d0
                 (cl-mpm::sim-mass-scale sim)
                 (cl-mpm/setup::estimate-critical-damping sim))
              ))
      (let ((dt-scale 0.5d0))
        (setf
         (cl-mpm:sim-dt sim)
         (* dt-scale h
            (sqrt (cl-mpm::sim-mass-scale sim))
            (sqrt (/ density (cl-mpm/particle::mp-p-modulus (aref (cl-mpm:sim-mps sim) 0)))))))
      (setf (cl-mpm:sim-bcs sim)
            (cl-mpm/bc::make-outside-bc-var
             (cl-mpm:sim-mesh sim)
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 nil nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil 0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(0 0 nil)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))
             (lambda (i) (cl-mpm/bc:make-bc-fixed i '(nil nil 0)))))
      (let* ((terminus-size (+ (second block-size) (* 0 (first block-size))))
             (ocean-x 1000)
             ;; (ocean-y (+ h-y (* 0.0d0 terminus-size)))
             (ocean-y (+ h-y (- terminus-size 100d0)))
             (ocean-y (* (round ocean-y h-y) h-y))
            )

        (format t "Ocean level ~a~%" ocean-y)
        (defparameter *water-height* ocean-y)
        (defparameter *floor-bc*
          (cl-mpm/penalty::make-bc-penalty-point-normal
           sim
           (cl-mpm/utils:vector-from-list '(0d0 1d0 0d0))
           (cl-mpm/utils:vector-from-list (list 00d0 (* 2d0 h-y) 0d0))
           (* 1d9 0.1d0)
           0.0d0))
        (setf (cl-mpm::sim-bcs-force-list sim)
              (list
               (cl-mpm/bc:make-bcs-from-list
                (list
                 (cl-mpm/buoyancy::make-bc-buoyancy-clip
                  sim
                  ocean-y
                  *water-density*
                  (lambda (pos datum)
                    t
                    ))))
               (cl-mpm/bc:make-bcs-from-list
                (list *floor-bc*)
                )
               ))
        )
      (format t "End of setup-test-column~%")
      sim)))




;; (defmethod cl-mpm::post-stress-step (mesh (mp cl-mpm/particle::particle-mc) dt)
;;   (with-accessors ((ps cl-mpm/particle::mp-strain-plastic-vm)
;;                    (c cl-mpm/particle::mp-c))
;;       mp
;;     ;; (let ((rho_0 200d3)
;;     ;;       (rho_1 0.002d3)
;;     ;;       (soft 1d0))
;;     ;;   (setf c (max rho_1
;;     ;;                  (* rho_0 (exp (- (* soft ps)))))))
;;     )
;;   )

(defparameter *sim* nil)
(defparameter *run-sim* t)
(defparameter *t* 0)
(defparameter *sim-step* 0)
(defparameter *refine* (/ 1d0 2d0))
(let ((refine (uiop:getenv "REFINE")))
  (when refine
    (setf *refine* (parse-integer (uiop:getenv "REFINE")))  
    ))
;; (defmethod cl-mpm::update-node-forces ((sim cl-mpm::mpm-sim))
;;   (with-accessors ((damping cl-mpm::sim-damping-factor)
;;                    (mass-scale cl-mpm::sim-mass-scale)
;;                    (mesh cl-mpm::sim-mesh)
;;                    (dt cl-mpm::sim-dt))
;;       sim
;;     (cl-mpm::iterate-over-nodes
;;      mesh
;;      (lambda (node)
;;        (cl-mpm::calculate-forces-cundall node damping dt mass-scale)))))

(defun setup ()
  (let* ((mesh-size 10)
         (mps-per-cell 2)
         (slope 0d0)
         (shelf-height 400)
         (shelf-aspect 2)
         (runout-aspect 1)
         (shelf-length (* shelf-height shelf-aspect))
         (shelf-end-height (+ shelf-height (* (- slope) shelf-length)))
         (shelf-height-terminus shelf-height)
         (shelf-height shelf-end-height)
         (offset (list 0
                       (* 2 mesh-size))))
    (defparameter *removal-point* (- (+ shelf-length (* runout-aspect shelf-height)) (* 2 mesh-size)))
    (defparameter *sim*
      (setup-test-column (list (+ shelf-length (* runout-aspect shelf-height))
                               (* shelf-height 2)
                               )
                         (list shelf-length shelf-height)
                         offset
                         (/ 1 mesh-size) mps-per-cell))

    (print "Removing mps")
    (let* ((undercut-angle 00d0)
           (normal (magicl:from-list (list
                                      (cos (- (* pi (/ undercut-angle 180d0))))
                                      (sin (- (* pi (/ undercut-angle 180d0))))
                                      0d0) '(3 1))))
      (cl-mpm/setup::remove-sdf *sim* (lambda (p) (plane-point-sdf
                                     p
                                     normal
                                     (magicl:from-list (list shelf-length shelf-height 0)
                                                       '(3 1) :type 'double-float)
                                     ))))
    )
  (format t "MPs: ~D~%" (length (cl-mpm:sim-mps *sim*)))
                                        ;(loop for f in (uiop:directory-files (uiop:merge-pathnames* "./outframes/")) do (uiop:delete-file-if-exists f))
  (loop for f in (uiop:directory-files (uiop:merge-pathnames* *output-directory*)) do (uiop:delete-file-if-exists f))
  (defparameter *run-sim* t)
  (defparameter *t* 0)
  (defparameter *oobf* 0)
  (defparameter *energy* 0)
  (defparameter *sim-step* 0))

(defun domain-decompose (sim)
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (let ((dsize (floor (cl-mpi:mpi-comm-size))))
      (setf (cl-mpm/mpi::mpm-sim-mpi-domain-count sim) (list dsize 1 1)))
    (when (= rank 0)
      (format t "Sim MPs: ~a~%" (length (cl-mpm:sim-mps sim)))
      (format t "Decompose~%"))
    (cl-mpm/mpi::domain-decompose
     sim
     :domain-scaler
     (lambda (domain)
       (destructuring-bind (x y z) domain
         (let ((dnew (list (mapcar (lambda (p)
                                     (if (< p 1d0)
                                         (* p (/ 2 3))
                                         p)) x)
                           y z)))
           (format t "Domain ~A ~A~%" (list x y z) dnew)
           dnew))))
    (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps sim)))))

(defun mpi-loop ()
  (format t "Starting mpi~%")
  (let ((rank (cl-mpi:mpi-comm-rank)))
    (format t "Starting setup~%")
    (setup)
    (format t "Finished setup~%")
    (domain-decompose *sim*)
    ;(when (= rank 0))
    (format t "Rank ~D - Sim MPs: ~a~%" rank (length (cl-mpm:sim-mps *sim*)))
    ;; (let ((mp (aref (cl-mpm:sim-mps *sim*) 0)))
    ;;   (when (slot-exists-p mp 'cl-mpm/particle::local-length)
    ;;     (let ((dhalo-size (* 1 (cl-mpm/particle::mp-local-length (aref (cl-mpm:sim-mps *sim*) 0)))))
    ;;       (when (= rank 0)
    ;;         (format t "Min size ~A length scale ~F~%" (mapcar (lambda (x) (abs (reduce #'- x)))  (cl-mpm/mpi::mpm-sim-mpi-domain-bounds *sim*)) dhalo-size) )
    ;;       (setf (cl-mpm/mpi::mpm-sim-mpi-halo-damage-size *sim*) dhalo-size))))
    (when (= rank 0)
      (format t "Run mpi~%"))
    (run-mpi)
    (when (= rank 0)
      (format t "Done mpi~%"))
    )
  )

(defmacro rank-0-time (rank &rest body)
  `(if (= ,rank 0)
      (time
        (progn
          ,@body))
      (progn
        ,@body)))

(defun run-mpi ()
  ;; (change-class *sim* 'cl-mpm/damage::mpm-sim-damage)
  ;; (cl-mpm/output:save-vtk-mesh (merge-pathnames *output-directory* "mesh.vtk") *sim*)
  (let* ((rank (cl-mpi:mpi-comm-rank))
         (target-time 1d1)
         (target-time-original target-time)
         (mass-scale (cl-mpm::sim-mass-scale *sim*))
         (collapse-target-time 1d0)
         (collapse-mass-scale 1d0)
         (dt (cl-mpm:sim-dt *sim*))
         (substeps (floor target-time dt))
         (settle-steps 25)
         (damp-steps 20)
         (dt-scale 0.8d0)
         (dt-0 0d0)
         (damping-0
            (* 0.1d0
                ;(cl-mpm::sim-mass-scale *sim*)
                (cl-mpm/setup::estimate-critical-damping *sim*))
          )
         ;(damping-0 (cl-mpm:sim-damping-factor *sim*))
		 (damage-0
		   (cl-mpm/mpi:mpi-sum
		      (lparallel:pmap-reduce (lambda (mp)
			                             (*
				                            1d0
				                            (cl-mpm/particle::mp-mass mp)
				                            (cl-mpm/particle::mp-damage mp)))
			                           #'+ (cl-mpm:sim-mps *sim*)
			                           :initial-value 0d0))))
    (cl-mpm/output::save-vtk-nodes (merge-pathnames *output-directory* (format nil "sim_nodes_~2,'0d.vtk" rank)) *sim*)
    (setf (cl-mpm:sim-dt *sim*) (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
    (setf dt-0 (/ (cl-mpm:sim-dt *sim*) (sqrt (cl-mpm::sim-mass-scale *sim*))))

    (defparameter *data-damage* 0d0)
    (defparameter *data-energy* 0d0)
    (rank-0-time
      rank
      (cl-mpm::update-sim *sim*))

    (setf substeps (floor target-time (cl-mpm:sim-dt *sim*)))



    (when (= rank 0)
      (format t "Substeps ~D~%" substeps)
 	    (cl-mpm/output::save-simulation-parameters (merge-pathnames *output-directory* "settings.json")
                                             *sim*
                                             (list :dt target-time)) 
		(with-open-file (stream (merge-pathnames *output-directory* "timesteps.csv") :direction :output :if-exists :supersede)
			  (format stream "steps,time,damage,energy,oobf~%")))
    (time (loop for steps from 0 to 400
                while *run-sim*
                do
                   (progn
                     (when (= rank 0)
                       (format t "Step ~d ~%" steps))
                     (when (= (mod steps 1) 0)
                       (cl-mpm/output:save-vtk (merge-pathnames *output-directory* (format nil "sim_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                       (cl-mpm/output::save-vtk-nodes (merge-pathnames *output-directory* (format nil "sim_nodes_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)) 
                     ;; (cl-mpm/output::save-vtk-cells (merge-pathnames *output-directory* (format nil "sim_cells_~2,'0d_~5,'0d.vtk" rank *sim-step*)) *sim*)
                     (when (= rank 0)
                      (with-open-file (stream (merge-pathnames *output-directory* "timesteps.csv") :direction :output :if-exists :append)
                                                (format stream "~D,~f,~f,~f,~f~%"
                                                                             steps
                                                                             *t*
                                                        *data-damage*
                                                        *data-energy*
                                                        *oobf*)))
                     (let ((energy-estimate 0d0)
                           (work 0d0))
                       (rank-0-time
                        rank
                        (dotimes (i substeps)
                          (cl-mpm::update-sim *sim*)
                          ;; (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))

                          (incf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                          (incf *oobf* (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                          (incf energy-estimate (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*))
                          (setf *t* (+ *t* (cl-mpm::sim-dt *sim*)))))


                       (setf *oobf* (/ *oobf* substeps)
                             energy-estimate (/ energy-estimate substeps))
                       (setf energy-estimate (abs (/ energy-estimate work)))
                       ;; (setf work (cl-mpm/dynamic-relaxation::estimate-power-norm *sim*))
                       ;; (setf *oobf* (cl-mpm/dynamic-relaxation::estimate-oobf *sim*))
                       ;; (when (> work 0d0)
                       ;;   (setf energy-estimate (/ (cl-mpm/dynamic-relaxation::estimate-energy-norm *sim*) work)))

                       ;; (setf energy-estimate (/ energy-estimate substeps) 
                       ;;       *oobf* (/ *oobf* substeps) )

                       (setf *data-energy* energy-estimate)
                       (let ((damage-est
                               (- (cl-mpm/mpi:mpi-sum
                                    (lparallel:pmap-reduce (lambda (mp)
                                                             (* 
                                                               1d0
					                                           (cl-mpm/particle::mp-mass mp)
                                                               (cl-mpm/particle::mp-damage mp)))
                                                           #'+ (cl-mpm:sim-mps *sim*)
                                                           :initial-value 0d0))
                                  damage-0)))
                         (setf *data-damage* damage-est)
                         (when (= rank 0)
                           (format t "Total damage: ~E~%" damage-est)))

                       (when (= rank 0)
                         (format t "Energy estimate: ~E~%" energy-estimate)
                         (format t "OOBF estimate: ~E~%" *oobf*) )


                       (when (>= steps settle-steps)
                         (setf (cl-mpm::sim-enable-damage *sim*) t
                               ;dt-scale 0.7d0
                               )
                         (if (or 
                               (> energy-estimate 1d-3)
                               (> *oobf* 1d-2)
                               ;(> steps 200)
                               )
                           (progn
                             (when (= rank 0)
                               (format t "Collapse timestep~%"))
                             (setf
                               ;dt-scale 0.5d0
                              target-time collapse-target-time
                              (cl-mpm::sim-mass-scale *sim*) collapse-mass-scale))
                           (progn
                             (when (= rank 0)
                               (format t "Accelerate timestep~%"))
                             (setf
                              target-time 1d1
                              (cl-mpm::sim-mass-scale *sim*) 1d4)
                             ))))
                     (when (>= steps damp-steps)
                       (setf (cl-mpm:sim-damping-factor *sim*)
                             damping-0
                             ;; (* 0.01d0 damping-0)
                             ))
                     (print "calc adaptive")
                     (multiple-value-bind (dt-e substeps-e)
                         (cl-mpm:calculate-adaptive-time *sim* target-time :dt-scale dt-scale)
                      (when (= rank 0)
                        (format t "CFL dt estimate: ~f~%" dt-e)
                        (format t "CFL step count estimate: ~D~%" substeps-e))
                      (setf substeps substeps-e))
					           (let* (;(dt-est (cl-mpm/setup::estimate-elastic-dt *sim* :dt-scale dt-scale))
							              (dt-est (* dt-0 (sqrt (cl-mpm::sim-mass-scale *sim*))))
					 	                (substeps-est (floor target-time dt-est)))
					             (when (< substeps-est substeps)
					 	             (setf (cl-mpm:sim-dt *sim*) dt-est)
					 	             (setf substeps substeps-est)))
                     (when (= rank 0)
                       (format t "CFL dt estimate: ~f~%" (cl-mpm:sim-dt *sim*))
                       (format t "CFL step count estimate: ~D~%" substeps))
					 (cl-mpm/setup:remove-sdf
                       *sim*
                       (lambda (p)
                         (if (> (magicl:tref p 0 0) *removal-point*)
                             -1d0
                             1d0)))
                     (incf *sim-step*)
                     ;; (swank.live:update-swank)
                     )))))

(defparameter *output-directory* (merge-pathnames
                                  ;"./output/"
                                  "/nobackup/rmvn14/ham-slump/"
                                  ;; "/mnt/d/Temp/ham-slump/"
                                  ))
(ensure-directories-exist *output-directory*)
(let* ((omp-get (uiop:getenv "OMP_NUM_THREADS"))
       (omp-threads (if omp-get
                        (parse-integer omp-get)
                        8)))
  (setf lparallel:*kernel* (lparallel:make-kernel omp-threads :name "custom-kernel")))
;; ;(defparameter *run-sim* nil)
;; ;(setup)
;; ;(format t "MP count:~D~%" (length (cl-mpm:sim-mps *sim*)))
;; ;(run)

;; (format t "Running~%")
;; ;(setf lparallel:*kernel* (lparallel:make-kernel 32 :name "custom-kernel"))
(mpi-loop)

(print "Test")
