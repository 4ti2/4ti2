/* SWIG interface file for Allegro CL */

%module "4ti2";

%include "pointer-in-out.i"
%include "list-vector.i"

typedef int int32_t;
typedef long int int64_t;

%{
#include "4ti2/4ti2.h"
%}

%inline %{

struct _4ti2_state {};
struct _4ti2_matrix {};

%}

%typemap(lout) _4ti2_status
%{ (let ((status $body))
     (check-4ti2-status status)) %}

TYPEMAP_POINTER_INPUT_OUTPUT(_4ti2_state *, -4ti2-state)
%apply _4ti2_state** OUTPUT { _4ti2_state** state };

TYPEMAP_POINTER_INPUT_OUTPUT(_4ti2_matrix *, -4ti2-matrix)
%apply _4ti2_matrix** OUTPUT { _4ti2_matrix** matrix };

%apply int *OUTPUT { int32_t *value };

%apply int PARALLEL_LISTLENINPUT { int argc };
%apply char **PARALLEL_LISTINPUT { char **argv };

%include "4ti2/4ti2.h"

