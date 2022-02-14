

#ifndef WWgg_SPLIT
#define WWgg_SPLIT

#include <stdlib.h>
#include <math.h> 

namespace split {
    bool Ph1_phi(double phi){
        phi = abs(phi);
        phi *= 10000;
        phi -= floor(phi);
        int phi_int = int((phi*10));
        return phi_int%2 == 0;
    }

}

#endif
