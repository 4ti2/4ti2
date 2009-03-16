/* -----------------------------------------------------------------------------
 * See the LICENSE file for information on copyright, usage and redistribution
 * of SWIG, and the README file for authors - http://www.swig.org/release.html.
 *
 * list-vector.i
 *
 * Allegro CL typemaps for converting between C arrays and Lisp lists or vectors  
 * ----------------------------------------------------------------------------- */

/* Here is a macro that will define typemaps for converting between C
   arrays and Lisp lists or vectors when passing arguments to the C
   function.

   TYPEMAP_LIST_VECTOR_INPUT_OUTPUT(...)
   
   Supported calling conventions:

   func(int VECTORLENINPUT, [const] C_TYPE *VECTORINPUT)

       Lisp wrapper will take one argument, a vector.  A temporary C
       array of elements of type C_TYPE will be allocated and filled
       with the elements of the vectors, converted to C with the
       SCM_TO_C function.  Length and address of the array are passed
       to the C function.

       SCM_TYPE is used to describe the Lisp type of the elements in
       the Guile procedure documentation.
   
   func(int LISTLENINPUT, [const] C_TYPE *LISTINPUT)

       Likewise, but the Lisp wrapper will take one argument, a list.

   func(int *VECTORLENOUTPUT, C_TYPE **VECTOROUTPUT)

       Lisp wrapper will take no arguments.  Addresses of an integer
       and a C_TYPE * variable will be passed to the C function.  The
       C function is expected to return address and length of a
       freshly allocated array of elements of type C_TYPE through
       these pointers.  The elements of this array are converted to
       Lisp with the C_TO_SCM function and returned as a Lisp
       vector. 

       If the function has a void return value, the vector constructed
       by this typemap becomes the return value of the Lisp wrapper.
       Otherwise, the function returns multiple values.  (See
       the documentation on how to deal with multiple values.)

   func(int *LISTLENOUTPUT, C_TYPE **LISTOUTPUT)

       Likewise, but the Lisp wrapper will return a list instead of
       a vector.

   It is also allowed to use "size_t LISTLENINPUT" rather than "int
   LISTLENINPUT".  */

%define TYPEMAP_LIST_VECTOR_INPUT_OUTPUT(C_TYPE, FF_TYPE, LISP_TYPE)

%typemap(lin, numinputs=0) int VECTORLENINPUT, size_t VECTORLENINPUT
%{  (let (($out nil))
      (labels ((swig-update-parallel-vector-length (sequence)
		 (let ((length (length sequence)))
		   (if $out
		       (assert (= $out length))
		       (setf $out length)))))
        (declare (ignorable swig-update-parallel-vector-length))
	$body)) %}

%typemap(lin) C_TYPE *VECTORINPUT, const C_TYPE *VECTORINPUT
%{  (swig-update-parallel-vector-length $in)
    (with-c-vector ($out $in FF_TYPE LISP_TYPE)
      $body)  %}

%typemap(lin, numinputs=0) int LISTLENINPUT, size_t LISTLENINPUT
%{  (let (($out nil))
      (labels ((swig-update-parallel-list-length (sequence)
		 (let ((length (length sequence)))
		   (if $out
		       (assert (= $out length))
		       (setf $out length)))))
        (declare (ignorable swig-update-parallel-list-length))
	$body)) %}

%typemap(lin) C_TYPE *LISTINPUT, const C_TYPE *LISTINPUT
%{  (swig-update-parallel-list-length $in)
    (with-c-vector ($out $in FF_TYPE LISP_TYPE)
      $body)  %}

%typemap(lin, numinputs=0) int *LISTLENOUTPUT, size_t *LISTLENOUTPUT
%{  (ff:with-stack-fobject ($out :int)
      (labels ((swig-array-length ()
		 (ff:fslot-value $out)))
	$body)) %}

%typemap(lin, numinputs=0) C_TYPE **LISTOUTPUT
%{  (ff:with-stack-fobject ($out (quote (* (* FF_TYPE))))
      $body
      (let ((array (fslot-value $out))
	    (length (swig-array-length)))
	(loop for index from 0 below length
	   collect (ff:fslot-value-typed (quote (:array FF_TYPE)) :c array index) into result
	   finally (push result ACL_result))))  %}

%typemap(lin, numinputs=0) int *VECTORLENOUTPUT, size_t *VECTORLENOUTPUT
%{  (ff:with-stack-fobject ($out :int)
      (labels ((swig-array-length ()
		 (ff:fslot-value $out)))
	$body)) %}

%typemap(lin, numinputs=0) C_TYPE **VECTOROUTPUT
%{  (ff:with-stack-fobject ($out (quote (* (* FF_TYPE))))
      $body
      (let* ((array (fslot-value $out))
	     (length (swig-array-length))
	     (result (make-array length :fill-pointer nil :adjustable nil)))
	(loop for index from 0 below length	      
	      do (setf (aref result index)
	      	       (ff:fslot-value-typed (quote (:array FF_TYPE)) :c array index))
	      finally (push result ACL_result))))  %} 

%enddef

TYPEMAP_LIST_VECTOR_INPUT_OUTPUT(int, :int, (signed-byte 32))
TYPEMAP_LIST_VECTOR_INPUT_OUTPUT(double, :double, double-float)
TYPEMAP_LIST_VECTOR_INPUT_OUTPUT(char, :char, character)

/* Following is a macro that emits typemaps that are much more
   flexible.  (They are also messier.)  It supports multiple parallel
   lists and vectors (sharing one length argument each).

   TYPEMAP_PARALLEL_LIST_VECTOR_INPUT_OUTPUT(...)
   
   Supported calling conventions:

   func(int PARALLEL_VECTORLENINPUT, [const] C_TYPE *PARALLEL_VECTORINPUT, ...)  or
   func([const] C_TYPE *PARALLEL_VECTORINPUT, ..., int PARALLEL_VECTORLENINPUT)

   func(int PARALLEL_LISTLENINPUT, [const] C_TYPE *PARALLEL_LISTINPUT, ...) or
   func([const] C_TYPE *PARALLEL_LISTINPUT, ..., int PARALLEL_LISTLENINPUT)

   func(int *PARALLEL_VECTORLENOUTPUT, C_TYPE **PARALLEL_VECTOROUTPUT, ...) or
   func(C_TYPE **PARALLEL_VECTOROUTPUT, int *PARALLEL_VECTORLENOUTPUT, ...)

   func(int *PARALLEL_LISTLENOUTPUT, C_TYPE **PARALLEL_LISTOUTPUT) or
   func(C_TYPE **PARALLEL_LISTOUTPUT, int *PARALLEL_LISTLENOUTPUT)

   It is also allowed to use "size_t PARALLEL_LISTLENINPUT" rather than "int
   PARALLEL_LISTLENINPUT".  */

%define TYPEMAP_PARALLEL_LIST_VECTOR_INPUT_OUTPUT(C_TYPE, FF_TYPE, LISP_TYPE)

%typemap(lin, numinputs=0) int PARALLEL_VECTORLENINPUT, size_t PARALLEL_VECTORLENINPUT
%{  (let (($out nil))
      (labels ((swig-update-parallel-vector-length (sequence)
		 (let ((length (length sequence)))
		   (if $out
		       (assert (= $out length))
		       (setf $out length)))))
	$body)) %}

%typemap(lin) C_TYPE *PARALLEL_VECTORINPUT, const C_TYPE *PARALLEL_VECTORINPUT
%{  (swig-update-parallel-vector-length $in)
    (with-c-vector ($out $in FF_TYPE LISP_TYPE)
      $body)  %}

%typemap(lin, numinputs=0) int PARALLEL_LISTLENINPUT, size_t PARALLEL_LISTLENINPUT
%{  (let (($out nil))
      (labels ((swig-update-parallel-list-length (sequence)
		 (let ((length (length sequence)))
		   (if $out
		       (assert (= $out length))
		       (setf $out length)))))
        (declare (ignorable swig-update-parallel-list-length))
	$body)) %}

%typemap(lin) C_TYPE *PARALLEL_LISTINPUT, const C_TYPE *PARALLEL_LISTINPUT
%{  (swig-update-parallel-list-length $in)
    (with-c-vector ($out $in FF_TYPE LISP_TYPE)
      $body)  %}

%typemap(lin, numinputs=0) int *PARALLEL_LISTLENOUTPUT, size_t *PARALLEL_LISTLENOUTPUT, int *LISTLENOUTPUT_NOFREE, size_t *LISTLENOUTPUT_NOFREE
%{  (ff:with-stack-fobject ($out :int)
      (labels ((swig-array-length ()
		 (ff:fslot-value $out)))
	$body)) %}

%typemap(lin, numinputs=0) C_TYPE **PARALLEL_LISTOUTPUT,
  C_TYPE **LISTOUTPUT_NOFREE
%{  (ff:with-stack-fobject ($out (quote (* (* FF_TYPE))))
      $body
      (let ((array (fslot-value $out))
	    (length (swig-array-length)))
	(loop for index from 0 below length
	   collect (ff:fslot-value-typed (quote (:array FF_TYPE)) :c array index) into result
	   finally (push result ACL_result))))  %}

%enddef

TYPEMAP_PARALLEL_LIST_VECTOR_INPUT_OUTPUT(int, :int, (signed-byte 32))
TYPEMAP_PARALLEL_LIST_VECTOR_INPUT_OUTPUT(double, :double, double-float)
TYPEMAP_PARALLEL_LIST_VECTOR_INPUT_OUTPUT(char, :char, character)
TYPEMAP_PARALLEL_LIST_VECTOR_INPUT_OUTPUT(char *, (* :char), string)

/* Allegro CL insists on converting character arrays itself from (array character *),
   even if we pass :strings-convert nil to def-foreign-call. */

%typemap(lin) char *PARALLEL_VECTORINPUT, const char *PARALLEL_VECTORINPUT
%{  (swig-update-parallel-vector-length $in)
    (let (($out $in))
      $body)  %}

%typemap(lin) char *PARALLEL_LISTINPUT, const char *PARALLEL_LISTINPUT
%{  (swig-update-parallel-list-length $in)
    (let (($out (coerce $in (quote string))))
      $body)  %}

/* For char *, we need to do everything ourselves anyway. */
%typemap(lin) char **PARALLEL_VECTORINPUT, const char **PARALLEL_VECTORINPUT
%{  (swig-update-parallel-vector-length $in)
    (with-char*-vector ($out $in)
      $body)  %}

%typemap(lin) char **PARALLEL_LISTINPUT, const char **PARALLEL_LISTINPUT
%{  (swig-update-parallel-list-length $in)
    (with-char*-vector ($out $in)
      $body)  %}

