;;; Basic definitions for the 4ti2 interface

(in-package :4ti2)

(define-condition 4ti2-error (error)
  ((status :reader 4ti2-error-status :initarg :status :type integer))
  (:documentation "The return values of the 4ti2 API are translated
into conditions of type 4TI2-ERROR."))

(defun check-4ti2-status (status)
  (unless (= status 0 #||*-4ti2-ok*||#)
    (error '4ti2-error :status status)))

