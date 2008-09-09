#ifndef _4ti2_groebner__CircuitsAPI_
#define _4ti2_groebner__CircuitsAPI_

#include "groebner/QSolveAPI.h"

namespace _4ti2_ {

class CircuitsAPI : public QSolveAPI {
public:
    CircuitsAPI();
    virtual ~CircuitsAPI();

    virtual void compute();
    virtual void write(const char* basename);

protected:
    virtual void write_usage();
    virtual void write_output_files();
};

} // namespace _4ti2_

#endif
