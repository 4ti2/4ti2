;; originally from ACL gywopt 
(defpackage :swig-macros
  (:use :cl)
  (:export #:with-c-vector
	   #:with-int-vector
	   #:with-double-vector
	   #:with-char*-vector))

(in-package :swig-macros)

(defmacro with-c-vector ((out in ff-type lisp-type &key no-coercion-for-type) &body body)
  ;; The COERCE calls are necessary, but checking whether an integer is of type (SIGNED-BYTE 32)
  ;; is very expensive on Allegro CL.  A fast path is to check for FIXNUM -- then no coercion is needed. 
  ;; See profiling below.
  (if (not no-coercion-for-type)
      (setq no-coercion-for-type
	    (case ff-type
	      (:int 'fixnum)
	      (:double 'double-float))))
  `(macrolet ((coerce-to-type (value-form)
		`(let ((value ,value-form))
		   (if (typep value ',',no-coercion-for-type)
		       value
		       (coerce value ',',lisp-type)))))
     (etypecase ,in
       ((simple-array ,lisp-type (*))
	(let ((,out ,in))
	  ,@body))
       ((simple-array t (*))
	(let* ((length (length ,in))
	       (type (list :array ',ff-type (if (zerop length) 1 length))))
	  (ff:with-stack-fobject (,out type)
	    (loop for element across ,in
	       for index of-type fixnum from 0
	       do (setf (ff:fslot-value-typed '(:array ,ff-type) :foreign ,out index)
			(coerce-to-type element)))
	    ,@body)))
       (vector
	(let* ((length (length ,in))
	       (type (list :array ',ff-type (if (zerop length) 1 length))))
	  (ff:with-stack-fobject (,out type)
	    (loop for element across ,in
	       for index of-type fixnum from 0
	       do (setf (ff:fslot-value-typed '(:array ,ff-type) :foreign ,out index)
			(coerce-to-type element)))
	    ,@body)))
       ;;      (null
       ;;       (let ((,out nil))
       ;; 	,@body))
       (list
	(let* ((length (length ,in))
	       (,out (make-array length :element-type ',lisp-type 
				 :adjustable nil :fill-pointer nil)))
	  (declare (type (simple-array ,lisp-type (*)) ,out)
		   (optimize (safety 0) (speed 3)))
	  (loop for element in ,in
	     for index of-type fixnum from 0
	     do (setf (aref ,out index)
		      (coerce-to-type element)))
	  ,@body)))))

(defmacro with-int-vector ((out in) &body body)
  `(with-c-vector (,out ,in :int (signed-byte 32))
     ,@body))

(defmacro with-double-vector ((out in) &body body)
  `(with-c-vector (,out ,in :double double-float)
     ,@body))

(defmacro with-char*-vector ((out in) &body body)
  `(let* ((length (length ,in))
	  (type (list :array '(* :char) (if (zerop length) 1 length))))
     (ff:with-stack-fobject (,out type)
       (etypecase ,in
	 (vector
	  (loop for element across ,in
	     for index from 0
	     do (setf (ff:fslot-value-typed '(:array (* :char)) :foreign ,out index)
		      (if (null element)
			  0
			  (ff:string-to-char* element)))))
	 (list
	  (loop for element in ,in
	     for index from 0
	     do (setf (ff:fslot-value-typed '(:array (* :char)) :foreign ,out index)
		      (if (null element)
			  0
			  (ff:string-to-char* element))))))
       ,@body)))



#||
;; Profiling:

(let ((fun (compile nil (lambda () (swig-macros:with-c-vector (out (make-list 1000000 :initial-element 1) :double double-float :no-coercion-for-type double-float)
	   (ff:fslot-value-typed '(:array :double) :foreign out 0)))))) (prof:with-profiling (:start-sampling-p t) (dotimes (i 1000) (funcall fun))))

(let ((fun (compile nil (lambda () (swig-macros:with-c-vector (out (make-list 1000000 :initial-element 1.0) :double double-float :no-coercion-for-type double-float)
	   (ff:fslot-value-typed '(:array :double) :foreign out 0)))))) (prof:with-profiling (:start-sampling-p t) (dotimes (i 172) (funcall fun))))

(let ((fun (compile nil (lambda () (swig-macros:with-c-vector (out (make-list 1000000 :initial-element 1) :int (signed-byte 32) :no-coercion-for-type fixnum)
	   (ff:fslot-value-typed '(:array :int) :foreign out 0)))))) (prof:with-profiling (:start-sampling-p t) (dotimes (i 172) (funcall fun))))

(let ((fun (compile nil (lambda () (swig-macros:with-c-vector (out (make-list 1000000 :initial-element 1) :int (signed-byte 32) :no-coercion-for-type nil)
	   (ff:fslot-value-typed '(:array :int) :foreign out 0)))))) (prof:with-profiling (:start-sampling-p t) (dotimes (i 172) (funcall fun))))


||#
