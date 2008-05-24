/*
4ti2 -- A software package for algebraic, geometric and combinatorial
problems on linear spaces.

Copyright (C) 2006 4ti2 team.
Main author(s): Matthias Walter.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. 
*/

#include "zsolve.h"
#include "Algorithm.hpp"

using namespace _4ti2_zsolve_;




class LibraryMatrix
{
private:
    int width_;
    int height_;
    bool read_only_;
    int internal_precision_;

protected:
    virtual int32_t get32_ (int r, int c) = 0;
    virtual int64_t get64_ (int r, int c) = 0;
    virtual void set32_ (int r, int c, int32_t value) = 0;
    virtual void set64_ (int r, int c, int64_t value) = 0;

#ifdef _4ti2_GMP_
    virtual void getGMP_ (mpz_ptr result, int r, int c) = 0;
    virtual void setGMP_ (int r, int c, mpz_srcptr value) = 0;
#endif

public:
    LibraryMatrix (int h, int w, bool ro, int prec)
    {
        width_ = w;
        height_ = h;
        read_only_ = ro;
        internal_precision_ = prec;
    }

    virtual ~LibraryMatrix ()
    {
	
    }

    int getWidth ()
    {
        return width_;
    }

    int getHeight ()
    {
        return height_;
    }

    bool isReadOnly ()
    {
        return read_only_;
    }

    int getPrecision ()
    {
        return internal_precision_;
    }

    int32_t get32 (int r, int c)
    {
        if (internal_precision_ == 32)
            return get32_ (r, c);
        else
            return 0;
    }

    int64_t get64 (int r, int c)
    {
        if (internal_precision_ == 32)
            return get32_ (r, c);
        else if (internal_precision_ == 64)
            return get64_ (r, c);
        else
            return 0;
    }

    bool set32 (int r, int c, int32_t value)
    {
        if (read_only_)
            return false;
        if (internal_precision_ == 32)
            set32_ (r, c, value);
        else if (internal_precision_ == 64)
            set64_ (r, c, value);
        else
	{
#ifndef _4ti2_GMP_
	    return false;
#else
	    mpz_t mp;
	    mpz_init_set_si (mp, value);
	    setGMP_ (r, c, mp);
	    mpz_clear (mp);
#endif
	}
        return true;
    }

    bool set64 (int r, int c, int64_t value)
    {
        if (read_only_)
            return false;
        if (internal_precision_ == 32)
            return false;
        else if (internal_precision_ == 64)
            set64_ (r, c, value);
        else
        {
#ifndef _4ti2_GMP_
	    return false;
#else
	    mpz_t mp;
	    mpz_init_set_si (mp, value);
	    setGMP_ (r, c, mp);
	    mpz_clear (mp);
#endif
	}
        return true;
    }

#ifdef _4ti2_GMP_
    void getGMP (mpz_ptr result, int r, int c)
    {
        if (internal_precision_ == 32)
            mpz_set_si (result, get32_ (r,c));
        else if (internal_precision_ == 64)
            mpz_set_si (result, get64_ (r,c));
        else
            getGMP_ (result, r, c);
    }

    bool setGMP (int r, int c, mpz_srcptr value)
    {
        if (read_only_)
            return false;
        if (internal_precision_ == 32 || internal_precision_ == 64)
            return false;
        setGMP_(r, c, value);
        return true;
    }
#endif
};








class LibraryMatMatrix : LibraryMatrix
{
private:
    LinearSystem <int32_t>* system32_;
    LinearSystem <int64_t>* system64_;
#ifdef _4ti2_GMP_
    LinearSystem <mpz_class>* systemGMP_;
#endif

public:
    LibraryMatMatrix (LinearSystem <int32_t>* system) : LibraryMatrix (system->relations(), system->variables (), false, 32)
    {
	system32_ = system;
	system64_ = NULL;
#ifdef _4ti2_GMP_
	systemGMP_ = NULL;
#endif
    }

    LibraryMatMatrix (LinearSystem <int64_t>* system) : LibraryMatrix (system->relations(), system->variables (), false, 64)
    {
	system32_ = NULL;
        system64_ = system;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
    }

    LibraryMatMatrix (LinearSystem <mpz_class>* system) : LibraryMatrix (system->relations(), system->variables (), false, 0)
    {
	system32_ = NULL;
        system64_ = NULL;
        systemGMP_ = system;
#endif
    }

protected:
    int32_t get32_ (int r, int c)
    {
	return system32_->matrix ()[r][c];
    }

    int64_t get64_ (int r, int c)
    {
	return system64_->matrix ()[r][c];
    }

    void set32_ (int r, int c, int32_t value)
    {
	system32_->matrix ()[r][c] = value;
    }

    void set64_ (int r, int c, int64_t value)
    {
	system64_->matrix ()[r][c] = value;
    }

#ifdef _4ti2_GMP_

    void getGMP_ (mpz_ptr result, int r, int c)
    {
	mpz_set (result, systemGMP_->matrix ()[r][c].get_mpz_t ());
    }

    void setGMP_ (int r, int c, mpz_srcptr value)
    {
	mpz_class v (value);
	systemGMP_->matrix ()[r][c] = v;
    }
#endif
};





class LibraryRhsMatrix : LibraryMatrix
{
private:
    LinearSystem <int32_t>* system32_;
    LinearSystem <int64_t>* system64_;
#ifdef _4ti2_GMP_
    LinearSystem <mpz_class>* systemGMP_;
#endif

public:
    LibraryRhsMatrix (LinearSystem <int32_t>* system) : LibraryMatrix (system->relations(), 1, false, 32)
    {
        system32_ = system;
        system64_ = NULL;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
#endif
    }

    LibraryRhsMatrix (LinearSystem <int64_t>* system) : LibraryMatrix (system->relations(), 1, false, 64)
    {
        system32_ = NULL;
        system64_ = system;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
    }

    LibraryRhsMatrix (LinearSystem <mpz_class>* system) : LibraryMatrix (system->relations(), 1, false, 0)
    {
        system32_ = NULL;
        system64_ = NULL;
        systemGMP_ = system;
#endif
    }

protected:
    int32_t get32_ (int r, int c)
    {
	return system32_->rhs ()[r];
    }

    int64_t get64_ (int r, int c)
    {
        return system64_->rhs ()[r];
    }

    void set32_ (int r, int c, int32_t value)
    {
        system32_->rhs ()[r] = value;
    }

    void set64_ (int r, int c, int64_t value)
    {
        system64_->rhs ()[r] = value;
    }

#ifdef _4ti2_GMP_

    void getGMP_ (mpz_ptr result, int r, int c)
    {
        mpz_set (result, systemGMP_->rhs ()[r].get_mpz_t ());
    }

    void setGMP_ (int r, int c, mpz_srcptr value)
    {
        mpz_class v (value);
        systemGMP_->rhs ()[r] = v;
    }
#endif
};




class LibraryRelMatrix : LibraryMatrix
{
private:
    LinearSystem <int32_t>* system32_;
    LinearSystem <int64_t>* system64_;
#ifdef _4ti2_GMP_
    LinearSystem <mpz_class>* systemGMP_;
#endif

public:
    LibraryRelMatrix (LinearSystem <int32_t>* system) : LibraryMatrix (system->relations(), 1, false, 32)
    {
        system32_ = system;
        system64_ = NULL;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
#endif
    }

    LibraryRelMatrix (LinearSystem <int64_t>* system) : LibraryMatrix (system->relations(), 1, false, 32)
    {
        system32_ = NULL;
        system64_ = system;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
    }

    LibraryRelMatrix (LinearSystem <mpz_class>* system) : LibraryMatrix (system->relations(), 1, false, 32)
    {
        system32_ = NULL;
        system64_ = NULL;
        systemGMP_ = system;
#endif
    }

protected:
    int32_t get32_ (int r, int c)
    {
	if (system32_ != NULL)
	    return (int32_t) system32_->get_relation (r).get ();
	else if (system64_ != NULL)
	    return (int32_t) system64_->get_relation (r).get ();
#ifdef _4ti2_GMP_
	else if (systemGMP_ != NULL)
            return (int32_t) systemGMP_->get_relation (r).get ();
#endif
	return 0;
    }

    int64_t get64_ (int r, int c)
    {
	return 0;
    }

    void set32_ (int r, int c, int32_t value)
    {
	if (system32_ != NULL)
	    system32_->get_relation (r).set ( (Relation <int32_t>::RelationType) value );
	else if (system64_ != NULL)
            system64_->get_relation (r).set ( (Relation <int64_t>::RelationType) value );
#ifdef _4ti2_GMP_
	else if (systemGMP_ != NULL)
            systemGMP_->get_relation (r).set ( (Relation <mpz_class>::RelationType) value );
#endif
    }

    void set64_ (int r, int c, int64_t value)
    {

    }

#ifdef _4ti2_GMP_

    void getGMP_ (mpz_ptr result, int r, int c)
    {

    }

    void setGMP_ (int r, int c, mpz_srcptr value)
    {

    }
#endif
};






class LibrarySignMatrix : LibraryMatrix
{
private:
    LinearSystem <int32_t>* system32_;
    LinearSystem <int64_t>* system64_;
#ifdef _4ti2_GMP_
    LinearSystem <mpz_class>* systemGMP_;
#endif

public:
    LibrarySignMatrix (LinearSystem <int32_t>* system) : LibraryMatrix (1, system->variables(), false, 32)
    {
        system32_ = system;
        system64_ = NULL;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
#endif
    }

    LibrarySignMatrix (LinearSystem <int64_t>* system) : LibraryMatrix (1, system->variables(), false, 32)
    {
        system32_ = NULL;
        system64_ = system;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
    }

    LibrarySignMatrix (LinearSystem <mpz_class>* system) : LibraryMatrix (1, system->variables(), false, 32)
    {
        system32_ = NULL;
        system64_ = NULL;
        systemGMP_ = system;
#endif
    }

protected:
    int32_t get32_ (int r, int c)
    {	
        if (system32_ != NULL)
	    return (int32_t) system32_->get_variable (c).sign ();
        else if (system64_ != NULL)
	    return (int32_t) system64_->get_variable (c).sign ();
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
	    return systemGMP_->get_variable (c).sign ();
#endif
        return 0;
    }

    int64_t get64_ (int r, int c)
    {
        return 0;
    }

    void set32_ (int r, int c, int32_t value)
    {
        if (system32_ != NULL)
	    system32_->get_variable (c).setSign (value);
        else if (system64_ != NULL)
            system64_->get_variable (c).setSign (value);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            systemGMP_->get_variable (c).setSign (value);
#endif
    }

    void set64_ (int r, int c, int64_t value)
    {

    }

#ifdef _4ti2_GMP_

    void getGMP_ (mpz_ptr result, int r, int c)
    {

    }

    void setGMP_ (int r, int c, mpz_srcptr value)
    {

    }
#endif
};


class LibraryBoundMatrix : LibraryMatrix
{
private:
    LinearSystem <int32_t>* system32_;
    LinearSystem <int64_t>* system64_;
#ifdef _4ti2_GMP_
    LinearSystem <mpz_class>* systemGMP_;
#endif
    bool is_lower_;

public:
    LibraryBoundMatrix (LinearSystem <int32_t>* system, bool is_lower) : LibraryMatrix (1, system->variables(), false, 32)
    {
        system32_ = system;
        system64_ = NULL;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
#endif
	is_lower_ = is_lower;
    }

    LibraryBoundMatrix (LinearSystem <int64_t>* system, bool is_lower) : LibraryMatrix (1, system->variables(), false, 64)
    {
        system32_ = NULL;
        system64_ = system;
#ifdef _4ti2_GMP_
        systemGMP_ = NULL;
	is_lower_ = is_lower;
    }

    LibraryBoundMatrix (LinearSystem <mpz_class>* system, bool is_lower) : LibraryMatrix (1, system->variables(), false, 0)
    {
        system32_ = NULL;
        system64_ = NULL;
        systemGMP_ = system;
#endif
	is_lower_ = is_lower;
    }

protected:
    int32_t get32_ (int r, int c)
    {
	if (is_lower_)
	    return system32_->get_variable (c).lower ();
	else
	    return system32_->get_variable (c).upper ();
    }

    int64_t get64_ (int r, int c)
    {
	if (is_lower_)
            return system64_->get_variable (c).lower ();
        else
            return system64_->get_variable (c).upper ();
    }

    void set32_ (int r, int c, int32_t value)
    {
	if (is_lower_)
            system32_->get_variable (c).setLower (value);
	else
	    system32_->get_variable (c).setUpper (value);
    }

    void set64_ (int r, int c, int64_t value)
    {
	if (is_lower_)
            system64_->get_variable (c).setLower (value);
        else
            system64_->get_variable (c).setUpper (value);

    }

#ifdef _4ti2_GMP_

    void getGMP_ (mpz_ptr result, int r, int c)
    {
	mpz_class value = is_lower_ ? systemGMP_->get_variable (c).lower () : systemGMP_->get_variable (c).upper ();
	mpz_set (result, value.get_mpz_t ());
    }

    void setGMP_ (int r, int c, mpz_srcptr value)
    {
	mpz_class v (value);
	if (is_lower_)
	    systemGMP_->get_variable (c).setLower (v);
	else
	    systemGMP_->get_variable (c).setUpper (v);
    }
#endif
};





class LibraryResultMatrix : LibraryMatrix
{
private:
    VectorArray <int32_t>* array32_;
    VectorArray <int64_t>* array64_;
#ifdef _4ti2_GMP_
    VectorArray <mpz_class>* arrayGMP_;
#endif

public:
    LibraryResultMatrix (VectorArray <int32_t>* array) : LibraryMatrix (array->height(), array->width (), true, 32)
    {
	array32_ = array;
	array64_ = NULL;
#ifdef _4ti2_GMP_
	arrayGMP_ = NULL;
#endif
    }

    LibraryResultMatrix (VectorArray <int64_t>* array) : LibraryMatrix (array->height(), array->width (), true, 64)
    {
	array32_ = NULL;
        array64_ = array;
#ifdef _4ti2_GMP_
        arrayGMP_ = NULL;
    }

    LibraryResultMatrix (VectorArray <mpz_class>* array) : LibraryMatrix (array->height(), array->width (), true, 0)
    {
	array32_ = NULL;
        array64_ = NULL;
        arrayGMP_ = array;
#endif
    }

protected:
    int32_t get32_ (int r, int c)
    {
        return (*array32_)[r][c];
    }

    int64_t get64_ (int r, int c)
    {
	return (*array64_)[r][c];
    }

    void set32_ (int r, int c, int32_t value)
    {

    }

    void set64_ (int r, int c, int64_t value)
    {

    }

#ifdef _4ti2_GMP_

    void getGMP_ (mpz_ptr result, int r, int c)
    {
        mpz_set (result, (*arrayGMP_)[r][c].get_mpz_t ());
    }

    void setGMP_ (int r, int c, mpz_srcptr value)
    {
	
    }
#endif
};



class LibraryState
{
private:
    int height_;
    int width_;
    LinearSystem <int32_t>* system32_;
    LinearSystem <int64_t>* system64_;
    Algorithm <int32_t>* algorithm32_;
    Algorithm <int64_t>* algorithm64_;
    VectorArray <int32_t>* inhoms32_;
    VectorArray <int32_t>* homs32_;
    VectorArray <int32_t>* frees32_;
    VectorArray <int64_t>* inhoms64_;
    VectorArray <int64_t>* homs64_;
    VectorArray <int64_t>* frees64_;

#ifdef _4ti2_GMP_
    LinearSystem <mpz_class>* systemGMP_;
    Algorithm <mpz_class>* algorithmGMP_;
    VectorArray <mpz_class>* inhomsGMP_;
    VectorArray <mpz_class>* homsGMP_;
    VectorArray <mpz_class>* freesGMP_;
#endif

public:
    LibraryState (int height, int width, int precision)
    {
        height_ = height;
        width_ = width;
        system32_ = NULL; system64_ = NULL;
	algorithm32_ = NULL; algorithm64_ = NULL;
	inhoms32_ = homs32_ = frees32_ = NULL;
        inhoms64_ = homs64_ = frees64_ = NULL;
#ifdef _4ti2_GMP_
	systemGMP_ = NULL;
	algorithmGMP_ = NULL;
	inhomsGMP_ = homsGMP_ = freesGMP_ = NULL;
#endif

        if (precision == 32)
	{
	    VectorArray <int32_t> array (height, width);
	    int32_t* rhs = create_vector <int32_t> (height);
	    system32_ = new LinearSystem <int32_t> (array, rhs, true, 1, -1);
	    delete_vector <int32_t> (rhs);
	}
        else if (precision == 64)
	{
	    VectorArray <int64_t> array (height, width);
	    int64_t* rhs = create_vector <int64_t> (height);
	    system64_ = new LinearSystem <int64_t> (array, rhs, true, 1, -1);
	    delete_vector <int64_t> (rhs);
	}
        else
	{
#ifdef _4ti2_GMP_
	    VectorArray <mpz_class> array (height, width);
	    mpz_class* rhs = create_vector <mpz_class> (height);
	    systemGMP_ = new LinearSystem <mpz_class> (array, rhs, true, 1, -1);
	    delete_vector <mpz_class> (rhs);
#endif
	}
    }

    ~LibraryState ()
    {
        if (system32_ != NULL) delete system32_;
        if (system64_ != NULL) delete system64_;
        if (algorithm32_ != NULL) delete algorithm32_;
        if (algorithm64_ != NULL) delete algorithm64_;
	if (inhoms32_ != NULL) delete inhoms32_;
	if (homs32_ != NULL) delete homs32_;
	if (frees32_ != NULL) delete frees32_;
	if (inhoms64_ != NULL) delete inhoms64_;
        if (homs64_ != NULL) delete homs64_;
        if (frees64_ != NULL) delete frees64_;
#ifdef _4ti2_GMP_
	if (systemGMP_ != NULL) delete systemGMP_;
	if (algorithmGMP_ != NULL) delete algorithmGMP_;
	if (inhomsGMP_ != NULL) delete inhomsGMP_;
        if (homsGMP_ != NULL) delete homsGMP_;
        if (freesGMP_ != NULL) delete freesGMP_;
#endif
    }

    void compute ()
    {
	int variables;
	if (system32_ != NULL)
	{
	    algorithm32_ = new Algorithm <int32_t> (system32_, NULL);
	    algorithm32_->compute ();
	    variables = algorithm32_->get_result_variables ();
	    inhoms32_ = new VectorArray <int32_t> (variables);
	    homs32_ = new VectorArray <int32_t> (variables);
	    frees32_ = new VectorArray <int32_t> (variables);
	    algorithm32_->extract_zsolve_results (*inhoms32_, *homs32_, *frees32_);
	}
	else if (system64_ != NULL)
	{
	    algorithm64_ = new Algorithm <int64_t> (system64_, NULL);
	    algorithm64_->compute ();
	    variables = algorithm64_->get_result_variables ();
            inhoms64_ = new VectorArray <int64_t> (variables);
            homs64_ = new VectorArray <int64_t> (variables);
            frees64_ = new VectorArray <int64_t> (variables);
            algorithm64_->extract_zsolve_results (*inhoms64_, *homs64_, *frees64_);
	}
#ifdef _4ti2_GMP_
	else if (systemGMP_ != NULL)
	{    
	    algorithmGMP_ = new Algorithm <mpz_class> (systemGMP_, NULL);
	    algorithmGMP_->compute ();
	    variables = algorithmGMP_->get_result_variables ();
            inhomsGMP_ = new VectorArray <mpz_class> (variables);
            homsGMP_ = new VectorArray <mpz_class> (variables);
            freesGMP_ = new VectorArray <mpz_class> (variables);
            algorithmGMP_->extract_zsolve_results (*inhomsGMP_, *homsGMP_, *freesGMP_);
	}
#endif
    }

    LibraryMatMatrix* mat ()
    {
        if (system32_ != NULL)
            return new LibraryMatMatrix (system32_);
        else if (system64_ != NULL)
            return new LibraryMatMatrix (system64_);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            return new LibraryMatMatrix (systemGMP_);
#endif
        else
            return NULL;
    }

    LibraryRhsMatrix* rhs ()
    {
	if (system32_ != NULL)
            return new LibraryRhsMatrix (system32_);
        else if (system64_ != NULL)
            return new LibraryRhsMatrix (system64_);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            return new LibraryRhsMatrix (systemGMP_);
#endif
        else
            return NULL;
    }

    LibraryRelMatrix* rel ()
    {
        if (system32_ != NULL)
            return new LibraryRelMatrix (system32_);
        else if (system64_ != NULL)
            return new LibraryRelMatrix (system64_);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            return new LibraryRelMatrix (systemGMP_);
#endif
        else
            return NULL;
    }

    LibrarySignMatrix* sign ()
    {
	if (system32_ != NULL)
            return new LibrarySignMatrix (system32_);
        else if (system64_ != NULL)
            return new LibrarySignMatrix (system64_);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            return new LibrarySignMatrix (systemGMP_);
#endif
        else
            return NULL;
    }

    LibraryBoundMatrix* lb ()
    {
	if (system32_ != NULL)
            return new LibraryBoundMatrix (system32_, true);
        else if (system64_ != NULL)
            return new LibraryBoundMatrix (system64_, true);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            return new LibraryBoundMatrix (systemGMP_, true);
#endif
        else
            return NULL;
    }

    LibraryBoundMatrix* ub ()
    {
	if (system32_ != NULL)
            return new LibraryBoundMatrix (system32_, false);
        else if (system64_ != NULL)
            return new LibraryBoundMatrix (system64_, false);
#ifdef _4ti2_GMP_
        else if (systemGMP_ != NULL)
            return new LibraryBoundMatrix (systemGMP_, false);
#endif
        else
	    return NULL;
    }

    LibraryResultMatrix* zfree ()
    {
	if (frees32_ != NULL)
	    return new LibraryResultMatrix (frees32_);
	else if (frees64_ != NULL)
	    return new LibraryResultMatrix (frees64_);
#ifdef _4ti2_GMP_
	else if (freesGMP_ != NULL)
	    return new LibraryResultMatrix (freesGMP_);
#endif
	return NULL;
    }

    LibraryResultMatrix* zhom ()
    {
	if (homs32_ != NULL)
            return new LibraryResultMatrix (homs32_);
        else if (homs64_ != NULL)
            return new LibraryResultMatrix (homs64_);
#ifdef _4ti2_GMP_
        else if (homsGMP_ != NULL)
            return new LibraryResultMatrix (homsGMP_);
#endif
        return NULL;
    }

    LibraryResultMatrix* zinhom ()
    {
	if (inhoms32_ != NULL)
            return new LibraryResultMatrix (inhoms32_);
        else if (inhoms64_ != NULL)
            return new LibraryResultMatrix (inhoms64_);
#ifdef _4ti2_GMP_
        else if (inhomsGMP_ != NULL)
            return new LibraryResultMatrix (inhomsGMP_);
#endif
        return NULL;
    }

};













extern "C"
{

    ZSolveState zsolve_state_create (int height, int width, int precision)
    {
	LibraryState* state = new LibraryState (height, width, precision);
	return (ZSolveState) state;
    }

    void zsolve_state_delete (ZSolveState state)
    {
	delete ((LibraryState*) state);
    }

    void zsolve_state_compute (ZSolveState state)
    {
	((LibraryState*) state)->compute ();
    }

    ZSolveMatrix zsolve_state_mat (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->mat ();
    }

    ZSolveMatrix zsolve_state_rhs (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->rhs ();
    }

    ZSolveMatrix zsolve_state_rel (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->rel ();
    }

    ZSolveMatrix zsolve_state_sign (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->sign ();
    }

    ZSolveMatrix zsolve_state_lb (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->lb ();
    }

    ZSolveMatrix zsolve_state_ub (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->ub ();
    }

    ZSolveMatrix zsolve_state_zfree (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->zfree ();
    }

    ZSolveMatrix zsolve_state_zhom (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->zhom ();
    }

    ZSolveMatrix zsolve_state_zinhom (ZSolveState state)
    {
	return (ZSolveMatrix) ((LibraryState*) state)->zinhom ();
    }

    ZSolveMatrix zsolve_state_matrix (ZSolveState state, char* name)
    {
        if (!strcmp (name, "mat"))
            return zsolve_state_mat (state);
        else if (!strcmp (name, "rhs"))
            return zsolve_state_rhs (state);
        else if (!strcmp (name, "rel"))
            return zsolve_state_rel (state);
        else if (!strcmp (name, "sign"))
            return zsolve_state_sign (state);
        else if (!strcmp (name, "lb"))
            return zsolve_state_lb (state);
        else if (!strcmp (name, "ub"))
            return zsolve_state_ub (state);
        else if (!strcmp (name, "zinhom"))
            return zsolve_state_zinhom (state);
        else if (!strcmp (name, "zhom"))
            return zsolve_state_zhom (state);
        else if (!strcmp (name, "zfree"))
            return zsolve_state_zfree (state);
        else
            return NULL;
    }



    

    int zsolve_matrix_width (ZSolveMatrix matrix)
    {
	return ((LibraryMatrix*) matrix)->getWidth ();
    }

    int zsolve_matrix_height (ZSolveMatrix matrix)
    {
        return ((LibraryMatrix*) matrix)->getHeight ();
    }

    int zsolve_matrix_read_only (ZSolveMatrix matrix)
    {
	return ((LibraryMatrix*) matrix)->isReadOnly () ? 1 : 0;
    }

    void zsolve_matrix_delete (ZSolveMatrix matrix)
    {
	delete ((LibraryMatrix*) matrix);
    }

    int zsolve_matrix_set_32 (ZSolveMatrix matrix, int r, int c, int32_t value)
    {
	return ((LibraryMatrix*) matrix)->set32 (r,c,value);
    }

    int zsolve_matrix_set_64 (ZSolveMatrix matrix, int r, int c, int64_t value)
    {
	return ((LibraryMatrix*) matrix)->set64 (r,c, value);
    }

    int32_t zsolve_matrix_get_32 (ZSolveMatrix matrix, int r, int c)
    {
	return ((LibraryMatrix*) matrix)->get32 (r,c);
    }

    int64_t zsolve_matrix_get_64 (ZSolveMatrix matrix, int r, int c)
    {
	return ((LibraryMatrix*) matrix)->get64 (r,c);
    }

    void zsolve_matrix_print_32 (ZSolveMatrix matrix)
    {
	int h = zsolve_matrix_height (matrix);
	int w = zsolve_matrix_width (matrix);
	printf ("%d %d\n", h, w);
	for (int r = 0; r < h; r++)
	{
	    for (int c = 0; c < w; c++)
	    {
		printf ("%d ", zsolve_matrix_get_32 (matrix, r, c));
	    }
	    printf ("\n");
	}
    }

    void zsolve_matrix_print_64 (ZSolveMatrix matrix)
    {
        int h = zsolve_matrix_height (matrix);
        int w = zsolve_matrix_width (matrix);
        printf ("%d %d\n", h, w);
        for (int r = 0; r < h; r++)
	{
	    for (int c = 0; c < w; c++)
	    {
		printf ("%lld ", zsolve_matrix_get_64 (matrix, r, c));
	    }
	    printf ("\n");
	}
    }

#ifdef _4ti2_GMP_
    int zsolve_matrix_set_gmp (ZSolveMatrix matrix, int r, int c, mpz_srcptr value)
    {
	return ((LibraryMatrix*) matrix)->setGMP (r,c,value);
    }

    void zsolve_matrix_get_gmp (ZSolveMatrix matrix, mpz_ptr result, int r, int c)
    {
	((LibraryMatrix*) matrix)->getGMP (result, r, c);
    }

    void zsolve_matrix_print_gmp (ZSolveMatrix matrix)
    {
	int h = zsolve_matrix_height (matrix);
        int w = zsolve_matrix_width (matrix);
        printf ("%d %d\n", h, w);
        for (int r = 0; r < h; r++)
	{
	    for (int c = 0; c < w; c++)
	    {
		mpz_t value;
		mpz_init (value);
		zsolve_matrix_get_gmp (matrix, value, r, c);
		gmp_printf ("%Zd ", value);
		mpz_clear (value);
	    }
	    printf ("\n");
	}
    }
#endif

}

