#ifndef _CAPD_CONFIG1
#define _CAPD_CONFIG1

#include "capd/auxil/OutputStream.h"
//================================global variables/configuration====================================

//either use stack or heap, when using stack consider its limited size
#define _HEAP
#if defined(_STACK)
  #define _D 2*_m+2
#endif
#if defined(_HEAP)
  #define _D 0
#endif

#define __FILIB__ 0 //whether to use the filib library for interval arithmetic


#define __SMALL__ 1e-25 ///used in rough enclosure algorithm, small value used to increase in diameter set at the begining

#define __MAX_VALIDATE_STEPS__ 1000 ///maximal number of steps possible when validating T

///DEBUGGING LOGGERS
capd::auxil::OutputStream generalDebug(std::cout, false, true); //OutputStream for general data
capd::auxil::OutputStream enclosureDebug(std::cout, false, true); //OutputStream for enclosure data
capd::auxil::OutputStream tailDebug(std::cout, false, true); //OutputStream for tail data
capd::auxil::OutputStream NDebug(std::cout, false, true);
capd::auxil::OutputStream generalDebug2(std::cout, false, true);
capd::auxil::OutputStream fftDebug(std::cout, false, true);
extern inline void setLoggers(){
  generalDebug.logfile("general.txt", true);
  generalDebug.log = false;
  generalDebug2.logfile("general2.txt", true);
  generalDebug2.log = false;
  enclosureDebug.logfile("enclosure.txt", true);
  enclosureDebug.log = false;
  tailDebug.logfile("tail.txt", true);
  tailDebug.log = false;
  NDebug.logfile("N.txt", true);
  NDebug.log = false;
  fftDebug.logfile("fft.txt", true);
  fftDebug.log = false;
}

///verification flags
#define __VERIFY_ENCLOSURE__ 1 //whether to check if validated rough-enclosure satisfies theorem assumptions (is in fact enclosure)
                               //slower calculations, but with guaranteed enclosures.
#define __VERIFY_TAIL__ 1 //whether to check if validated tail satisfies theorem assumptions (is in fact enclosure and T([0,h])\subset T),
                          //slower calculations, but with guaranteed tail enclosures.


///tail validation flags
#define __DECREASE_L__ 1.02 ///used in the automatic choice of M procedure

#define __INCREASE_L__ 1.02 ///used in the automatic choice of M procedure

#define __L_THRESHOLD__ 1.05 ///used in the automatic choice of M procedure

///the tail validation constants

#define __DECREASE_CT_EPS__ 1.01

#define __D_STEP__ 0.01

#define __CT_EPS__ 1e+010

#define __L_CONST__ 2

#define __L_CONST_DECREASE__ 2

///decimal places upto which L is truncated, to determine if its value have stabilised
#define __L_TRUNCATION_DP__ 2

#endif //_CAPD_CONFIG1
