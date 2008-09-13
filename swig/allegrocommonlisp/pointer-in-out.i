/* -----------------------------------------------------------------------------
 * See the LICENSE file for information on copyright, usage and redistribution
 * of SWIG, and the README file for authors - http://www.swig.org/release.html.
 *
 * pointer-in-out.i
 *
 * Allegro CL typemaps for passing pointers indirectly 
 * ----------------------------------------------------------------------------- */

/* Here is a macro that will define typemaps for passing C pointers indirectly.
  
   TYPEMAP_POINTER_INPUT_OUTPUT(PTR_TYPE, PRETTY_TYPE)

   PTR_TYPE is a pointer types.
   PRETTY_TYPE is only used to described the type in the generated docstrings.
   Supported calling conventions (in this example, PTR_TYPE is int *):

   func(int **INPUT)

       Scheme wrapper will take one argument, a wrapped C pointer.
       The address of a variable containing this pointer will be
       passed to the function.

   func(int **INPUT_CONSUMED)			NOT YET IMPLEMENTED

       Likewise, but mark the pointer object as not garbage
       collectable.

   func(int **INPUT_DESTROYED)			NOT YET IMPLEMENTED

       Likewise, but mark the pointer object as destroyed.
       
   func(int **OUTPUT)

       Scheme wrapper will take no arguments.  The address of an int *
       variable will be passed to the function.  The function is
       expected to modify the variable; its value is wrapped and
       becomes an extra return value.  (See the documentation on how
       to deal with multiple values.)
   
   func(int **OUTPUT_NONCOLLECTABLE)

       Likewise, but make the pointer object not garbage collectable.
   
   func(int **BOTH)				NOT YET IMPLEMENTED
   func(int **INOUT)				NOT YET IMPLEMENTED

       This annotation combines INPUT and OUTPUT.

*/

%define TYPEMAP_POINTER_INPUT_OUTPUT(PTR_TYPE, LISP_CLASS)

%typemap(lin, numinputs=0) PTR_TYPE *OUTPUT
%{ (let (($out (ff:allocate-fobject '$*in_fftype :c)))
     $body
     (let* ((address (ff:fslot-value-typed (quote $*in_fftype) :c $out))
	    (object (make-instance (quote LISP_CLASS)
				   :foreign-address address)))
       (unless (zerop address)
	 (excl:schedule-finalization object (function $ldestructor)))
       (push object ACL_result)
       (ff:free-fobject $out))) %}

%typemap(lin, numinputs=0) PTR_TYPE *OUTPUT_NONCOLLECTABLE
%{ (let (($out (ff:allocate-fobject '$*in_fftype :c)))
     $body
     (let* ((address (ff:fslot-value-typed (quote $*in_fftype) :c $out))
	    (object (make-instance (quote LISP_CLASS)
				   :foreign-address address)))
       (push object ACL_result)
       (ff:free-fobject $out))) %}

%typemap(lin) PTR_TYPE *INPUT
%{ (with-stack-fobject ($out '$*in_fftype)
     (setf (fslot-value $out) (foreign-pointer-address $in))
     $body) %}

/* As long as SWIG does not get the return class of typedefs right: */
%typemap(lout) PTR_TYPE
%{ (let* ((address $body)
	  (new-inst (make-instance (quote LISP_CLASS)
				   :foreign-address address)))
    (when (and $owner (not (zerop address)))
      (excl:schedule-finalization new-inst (function $ldestructor)))
    (setq ACL_ffresult new-inst)) %}

%enddef

TYPEMAP_POINTER_INPUT_OUTPUT(int *, ff:foreign-pointer)
TYPEMAP_POINTER_INPUT_OUTPUT(double *, ff:foreign-pointer)
TYPEMAP_POINTER_INPUT_OUTPUT(void *, ff:foreign-pointer)
