#ifndef _4ti2_groebner__RaysAPI_
#define _4ti2_groebner__RaysAPI_

#include "groebner/QSolveAPI.h"

namespace _4ti2_ {

class RaysAPI : public QSolveAPI {
public:
    RaysAPI();
    virtual ~RaysAPI();

    virtual void compute();
    virtual void write(const char* basename);

protected:
    virtual void write_usage();
    virtual void write_output_files();
};

} // namespace _4ti2_

#endif
