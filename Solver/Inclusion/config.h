#ifndef _CAPD_CONFIG
#define _CAPD_CONFIG

//================================global variables/configuration====================================

//=======================================more important=============================================
//either use stack or heap, when using stack consider its limited size
#define _HEAP
#if defined(_STACK)
  #define _D 2*_m+2
#endif
#if defined(_HEAP)
  #define _D 0
#endif

#define __CONSTANT_M__ 0 //decides whether M is constant or is changing throughout integration

#define __M_MAX__ 500 //default optimal values of M

#define __S_DISSIPATIVE__ -1 //threshold value, starting at which lambda_k is considered as dissipative,
               //enclosure calculating algorithm is being switched from C^1 algorithm into
               //dissipative algorithm based on isolation

///files where debug data is written
#define __FILE_NAME_DATA__ "data.txt"
#define __FILE_NAME_DEBUG__ "debug.txt"
#define __FILE_NAME_SOLUTIONS__ "solutions.txt"
#define __FILE_NAME_n__ "n.txt"
#define __FILE_NAME_ENCLOSURES__ "temp.txt"
#define __FILE_NAME_SM__ "sm.txt"
#define __FILE_NAME_BASIN__ "basin.txt"
#define __FILE_NAME_BASIN_IN__ "basin_in.txt"

#define __FILE_NAME_EXACT_FIXED_POINT_LOCATION__ "proof_part2_exact_fixed_point_location.txt"
#define __FILE_NAME_TRAPPING_REGION__ "proof_part1_the_trapping_region.txt"
#define __FILE_NAME_L__ "proof_part3_the_l.txt"
#define __FILE_NAME_BASIN_OF_ATTRACTION__ "proof_part4_basin_of_attraction.txt"
#define __FILE_NAME_FORCING__ "forcing.txt"

#define __PRINT_DATA_IN_LATEX_FORMAT__ 0 ///if data is outputted in latex format

#define __SMALL__ 1e-25 ///used in rough enclosure algorithm, small value used to increase in diameter set at the begining

#define __REFINEMENT_STEPS__ 2 ///perform how many steps to refine enclosure and tail AT LEAST ONE STEP

#define __INITIAL_RADIUS__ 1e-12 ///tries to find a trapping region in a box that is product of [-__INITIAL_RADIUS__, __INITIAL_RADIUS__]
                                 ///if not succeed then the box is increased

#define __MAX_VALIDATE_STEPS__ 1000 ///maximal number of steps possible when validating T

///debugging flags
#define __GENERAL_DEBUG__ 0 //debugging flag, when on writing W_1, W_2, \delta, \Delta, x
#define __INFLATE_DEBUG__ 0 //debugging flag of inflate coefficients function, writes out the array with coefficients used to inflate the elements in the near tail
#define __DEBUG_ENCLOSURE__ 0 ///rough enclosure calculation debug loging all steps
#define __DEBUG_TAIL__ 0 ///tail debug, logging T, T[0], T[0, h], T[h)

///debugging flags of function "enclosureWithTail" (all steps used in order to validate T ( such that T[0, h]\subset T)
#define __ENCLOSURE_WITH_TAIL_DEBUG__ 0
#define __ENCLOSURE_WITH_TAIL_DEBUG_FAR_TAIL__ 0

#define __N_DEBUG__ 0 //debugging flag of function "N", calculating non-linear part for the projection (k\leq m)
#define __n_DEBUG__ 0 //debugging falg of function "n", calculating non-linear part for the tail (k>m)
#define __nk_DEBUG__ 0
#define __nDF_DEBUG__ 0
#define __DELTA_DEBUG__ 0 //debugging flag of function calculating \delta
#define __CALCULATE_TH_DEBUG__ 0 //debugging flag of function calculating T(h)
#define __D12inf_DEBUG__ 0 //debugging flag of function calculating constants D_F(k\leq 2M), D_F(k>2M), D_I
#define __BOX_DEBUG__ 0 //debugging flag of function calculating enclosing box of fixed point
#define __CALCULATE_L_DEBUG__ 0 //debugging flag of function calculating the Lipshitz constant of the flow

///verification flags
#define __VERIFY_ENCLOSURE__ 1 //whether to check if validated rough-enclosure satisfies theorem assumptions (is in fact enclosure)
                               //slower calculations, but with guaranteed enclosures.
#define __VERIFY_TAIL__ 1 //whether to check if validated tail satisfies theorem assumptions (is in fact enclosure and T([0,h])\subset T),
                          //slower calculations, but with guaranteed tail enclosures.


///tail validation flags
#define __SM__ 1 ///whether to log changes in the far tail exponent s and the dimension M

#define __DECREASE_L__ 1.02

#define __INCREASE_L__ 1.02

#define __L_THRESHOLD__ 1.05 ///choosing M procedure

///the tail validation constants
#define __D_G__ 0

#define __D_2__ 1.001

#define __DECREASE_CT_EPS__ 1.01

#define __D_STEP__ 0.01

#define __INFLATE_C__ 1.01

#define __INFLATE_RADIUS_1__ 5

#define __INFLATE_RADIUS_2__ 5

#define __CT_EPS__ 1e+010

#define __CT_MIN__ 1e-100 ///minimal value of C(T) that is considered, below this the tail T is not validated

#define __L_CONST__ 2

#define __L_CONST_DECREASE__ 2

///decimal places upto which L is truncated, to determine if its value have stabilised
#define __L_TRUNCATION_DP__ 2

//box finding flags
#define __BOX_FIND_INFLATE_C__ 1.01

#define __SPLIT_PIECES__ 2//defines how many smaller boxes are result of splitting after set blow up
#define __SET_BLOWUP_L_THRESHOLD__ 2000 //threshold value, when value of L reaches this limit we consider this as a set blow up

//======================================DOUBLE FUNCTIONS========================================
//Define double versions of the basic math functions, architecture dependent.
#define POW power
#define POW_DOUBLE powl
#define SIN sin
#define COS cos
#define EXP exp
#define EXP_DOUBLE expl
#define SQRT sqrtl
#define ROUND roundl
#define CEIL ceill
#define FLOOR floorl

//=======================================less important=============================================

///various default values
#define _m 5
#define _M 10
#define _ni 0.1
#define _STEP 0.001
#define __ORDER__ 5
#define __END_TIME_INTERVAL__ ap::pi()

const int _p = 2;
const int _s = 3;

//#define _p 2
//#define _s 3
#define _C 0
#define _R 1



#endif
