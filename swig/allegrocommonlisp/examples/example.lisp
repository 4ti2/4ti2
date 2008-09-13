(use-package :4ti2)

(let* ((state (-4ti2-rays-create-state *-4ti2-prec-int-32*))
       (mat (-4ti2-state-create-matrix state 
				       #|numrows:|# 2 
						    #|numcols:|# 3
								 "mat"))
       (rel (-4ti2-state-create-matrix state 1 2 "rel")))
  (-4ti2-matrix-set-entry-int32-t mat 0 2 47)
  (-4ti2-matrix-write-to-stdout mat)
  (-4ti2-matrix-write-to-stdout rel)
  (-4ti2-state-set-options state '("IGNORED-ARG" "--support"))
  (-4ti2-state-compute state)
  (let ((hom (-4ti2-state-get-matrix state "qhom")))
    ;;(-4ti2-matrix-write-to-stdout hom)
    (format t "Numrows: ~D Numcols: ~D~%" 
	    (-4ti2-matrix-get-num-rows hom)
	    (-4ti2-matrix-get-num-cols hom))
    )
  (-4ti2-state-delete state))

(let* ((state (-4ti2-rays-create-state *-4ti2-prec-int-arb*))
       ;; FIXME: Non-implemented precision does not signal an error!
       (mat (-4ti2-state-create-matrix state 
				       #|numrows:|# 2 
						    #|numcols:|# 3
								 "mat"))
       (rel (-4ti2-state-create-matrix state 1 2 "rel")))
  (-4ti2-matrix-set-entry-int32-t mat 0 2 47)
  (-4ti2-matrix-write-to-stdout mat)
  (-4ti2-matrix-write-to-stdout rel)
  (-4ti2-state-compute state)
  (let ((hom (-4ti2-state-get-matrix state "qhom")))
    (-4ti2-matrix-write-to-stdout hom))
  (-4ti2-state-delete state))

;;; Error handling:

(let* ((state (-4ti2-rays-create-state *-4ti2-prec-int-32*))
       (xyz (-4ti2-state-create-matrix state 1 2 "xyz")))
  ;; FIXME: Does not signal an error!!
  nil)

