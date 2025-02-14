#include <math.h>
#include <string.h>
#include "pm_std.h"
#include "sm_std.h"
#include "ne_std.h"
#include "ne_dae.h"
#include "sm_ssci_run_time_errors.h"
#include "sm_RuntimeDerivedValuesBundle.h"
#include "sm_CTarget.h"
void MagneticBasket_Simscape_Optimizer_dda62cd9_1_setTargets ( const
RuntimeDerivedValuesBundle * rtdv , CTarget * targets ) { ( void ) rtdv ; (
void ) targets ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_resetAsmStateVector ( const void
* mech , double * state ) { double xx [ 1 ] ; ( void ) mech ; xx [ 0 ] = 0.0
; state [ 0 ] = xx [ 0 ] ; state [ 1 ] = xx [ 0 ] ; state [ 2 ] = xx [ 0 ] ;
state [ 3 ] = xx [ 0 ] ; state [ 4 ] = xx [ 0 ] ; state [ 5 ] = xx [ 0 ] ;
state [ 6 ] = xx [ 0 ] ; state [ 7 ] = xx [ 0 ] ; state [ 8 ] = xx [ 0 ] ;
state [ 9 ] = xx [ 0 ] ; state [ 10 ] = xx [ 0 ] ; state [ 11 ] = xx [ 0 ] ;
state [ 12 ] = xx [ 0 ] ; state [ 13 ] = xx [ 0 ] ; state [ 14 ] = xx [ 0 ] ;
state [ 15 ] = xx [ 0 ] ; state [ 16 ] = xx [ 0 ] ; state [ 17 ] = xx [ 0 ] ;
} void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_initializeTrackedAngleState (
const void * mech , const RuntimeDerivedValuesBundle * rtdv , const int *
modeVector , const double * motionData , double * state ) { const double *
rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts .
mValues ; double xx [ 12 ] ; ( void ) mech ; ( void ) rtdvd ; ( void ) rtdvi
; ( void ) modeVector ; xx [ 0 ] = motionData [ 70 ] ; xx [ 1 ] = motionData
[ 71 ] ; xx [ 2 ] = motionData [ 72 ] ; xx [ 3 ] = motionData [ 73 ] ; xx [ 4
] = motionData [ 91 ] ; xx [ 5 ] = motionData [ 92 ] ; xx [ 6 ] = motionData
[ 93 ] ; xx [ 7 ] = motionData [ 94 ] ; pm_math_Quaternion_inverseCompose_ra
( xx + 0 , xx + 4 , xx + 8 ) ; xx [ 0 ] = motionData [ 168 ] ; xx [ 1 ] =
motionData [ 169 ] ; xx [ 2 ] = motionData [ 170 ] ;
pm_math_Quaternion_inverseXform_ra ( xx + 8 , xx + 0 , xx + 3 ) ; state [ 16
] = pm_math_canonicalAngle ( 2.0 * atan2 ( sqrt ( xx [ 9 ] * xx [ 9 ] + xx [
10 ] * xx [ 10 ] + xx [ 11 ] * xx [ 11 ] ) , fabs ( xx [ 8 ] ) ) * ( ( - ( xx
[ 8 ] * xx [ 10 ] ) ) < 0.0 ? - 1.0 : + 1.0 ) ) ; state [ 17 ] = xx [ 4 ] -
motionData [ 145 ] ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computeDiscreteState ( const
void * mech , const RuntimeDerivedValuesBundle * rtdv , double * state ) {
const double * rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv
-> mInts . mValues ; ( void ) mech ; ( void ) rtdvd ; ( void ) rtdvi ; ( void
) state ; } void MagneticBasket_Simscape_Optimizer_dda62cd9_1_adjustPosition
( const void * mech , const double * dofDeltas , double * state ) { ( void )
mech ; state [ 0 ] = state [ 0 ] + dofDeltas [ 0 ] ; state [ 2 ] = state [ 2
] + dofDeltas [ 1 ] ; state [ 4 ] = state [ 4 ] + dofDeltas [ 2 ] ; state [ 6
] = state [ 6 ] + dofDeltas [ 3 ] ; state [ 8 ] = state [ 8 ] + dofDeltas [ 4
] ; state [ 10 ] = state [ 10 ] + dofDeltas [ 5 ] ; state [ 12 ] = state [ 12
] + dofDeltas [ 6 ] ; state [ 14 ] = state [ 14 ] + dofDeltas [ 7 ] ; }
static void perturbAsmJointPrimitiveState_0_0 ( double mag , double * state )
{ state [ 0 ] = state [ 0 ] + mag ; } static void
perturbAsmJointPrimitiveState_0_0v ( double mag , double * state ) { state [
0 ] = state [ 0 ] + mag ; state [ 1 ] = state [ 1 ] - 0.875 * mag ; } static
void perturbAsmJointPrimitiveState_1_0 ( double mag , double * state ) {
state [ 2 ] = state [ 2 ] + mag ; } static void
perturbAsmJointPrimitiveState_1_0v ( double mag , double * state ) { state [
2 ] = state [ 2 ] + mag ; state [ 3 ] = state [ 3 ] - 0.875 * mag ; } static
void perturbAsmJointPrimitiveState_2_0 ( double mag , double * state ) {
state [ 4 ] = state [ 4 ] + mag ; } static void
perturbAsmJointPrimitiveState_2_0v ( double mag , double * state ) { state [
4 ] = state [ 4 ] + mag ; state [ 5 ] = state [ 5 ] - 0.875 * mag ; } static
void perturbAsmJointPrimitiveState_3_0 ( double mag , double * state ) {
state [ 6 ] = state [ 6 ] + mag ; } static void
perturbAsmJointPrimitiveState_3_0v ( double mag , double * state ) { state [
6 ] = state [ 6 ] + mag ; state [ 7 ] = state [ 7 ] - 0.875 * mag ; } static
void perturbAsmJointPrimitiveState_4_0 ( double mag , double * state ) {
state [ 8 ] = state [ 8 ] + mag ; } static void
perturbAsmJointPrimitiveState_4_0v ( double mag , double * state ) { state [
8 ] = state [ 8 ] + mag ; state [ 9 ] = state [ 9 ] - 0.875 * mag ; } static
void perturbAsmJointPrimitiveState_5_0 ( double mag , double * state ) {
state [ 10 ] = state [ 10 ] + mag ; } static void
perturbAsmJointPrimitiveState_5_0v ( double mag , double * state ) { state [
10 ] = state [ 10 ] + mag ; state [ 11 ] = state [ 11 ] - 0.875 * mag ; }
static void perturbAsmJointPrimitiveState_6_0 ( double mag , double * state )
{ state [ 12 ] = state [ 12 ] + mag ; } static void
perturbAsmJointPrimitiveState_6_0v ( double mag , double * state ) { state [
12 ] = state [ 12 ] + mag ; state [ 13 ] = state [ 13 ] - 0.875 * mag ; }
static void perturbAsmJointPrimitiveState_7_0 ( double mag , double * state )
{ state [ 14 ] = state [ 14 ] + mag ; } static void
perturbAsmJointPrimitiveState_7_0v ( double mag , double * state ) { state [
14 ] = state [ 14 ] + mag ; state [ 15 ] = state [ 15 ] - 0.875 * mag ; }
void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_perturbAsmJointPrimitiveState (
const void * mech , size_t stageIdx , size_t primIdx , double mag , boolean_T
doPerturbVelocity , double * state ) { ( void ) mech ; ( void ) stageIdx ; (
void ) primIdx ; ( void ) mag ; ( void ) doPerturbVelocity ; ( void ) state ;
switch ( ( stageIdx * 6 + primIdx ) * 2 + ( doPerturbVelocity ? 1 : 0 ) ) {
case 0 : perturbAsmJointPrimitiveState_0_0 ( mag , state ) ; break ; case 1 :
perturbAsmJointPrimitiveState_0_0v ( mag , state ) ; break ; case 12 :
perturbAsmJointPrimitiveState_1_0 ( mag , state ) ; break ; case 13 :
perturbAsmJointPrimitiveState_1_0v ( mag , state ) ; break ; case 24 :
perturbAsmJointPrimitiveState_2_0 ( mag , state ) ; break ; case 25 :
perturbAsmJointPrimitiveState_2_0v ( mag , state ) ; break ; case 36 :
perturbAsmJointPrimitiveState_3_0 ( mag , state ) ; break ; case 37 :
perturbAsmJointPrimitiveState_3_0v ( mag , state ) ; break ; case 48 :
perturbAsmJointPrimitiveState_4_0 ( mag , state ) ; break ; case 49 :
perturbAsmJointPrimitiveState_4_0v ( mag , state ) ; break ; case 60 :
perturbAsmJointPrimitiveState_5_0 ( mag , state ) ; break ; case 61 :
perturbAsmJointPrimitiveState_5_0v ( mag , state ) ; break ; case 72 :
perturbAsmJointPrimitiveState_6_0 ( mag , state ) ; break ; case 73 :
perturbAsmJointPrimitiveState_6_0v ( mag , state ) ; break ; case 84 :
perturbAsmJointPrimitiveState_7_0 ( mag , state ) ; break ; case 85 :
perturbAsmJointPrimitiveState_7_0v ( mag , state ) ; break ; } } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computePosDofBlendMatrix ( const
void * mech , size_t stageIdx , size_t primIdx , const double * state , int
partialType , double * matrix ) { ( void ) mech ; ( void ) stageIdx ; ( void
) primIdx ; ( void ) state ; ( void ) partialType ; ( void ) matrix ; switch
( ( stageIdx * 6 + primIdx ) ) { } } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computeVelDofBlendMatrix ( const
void * mech , size_t stageIdx , size_t primIdx , const double * state , int
partialType , double * matrix ) { ( void ) mech ; ( void ) stageIdx ; ( void
) primIdx ; ( void ) state ; ( void ) partialType ; ( void ) matrix ; switch
( ( stageIdx * 6 + primIdx ) ) { } } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_projectPartiallyTargetedPos (
const void * mech , size_t stageIdx , size_t primIdx , const double *
origState , int partialType , double * state ) { ( void ) mech ; ( void )
stageIdx ; ( void ) primIdx ; ( void ) origState ; ( void ) partialType ; (
void ) state ; switch ( ( stageIdx * 6 + primIdx ) ) { } } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_propagateMotion ( const void *
mech , const RuntimeDerivedValuesBundle * rtdv , const double * state ,
double * motionData ) { const double * rtdvd = rtdv -> mDoubles . mValues ;
const int * rtdvi = rtdv -> mInts . mValues ; double xx [ 85 ] ; ( void )
mech ; ( void ) rtdvd ; ( void ) rtdvi ; xx [ 0 ] = 0.7071067811865476 ; xx [
1 ] = - xx [ 0 ] ; xx [ 2 ] = 0.0 ; xx [ 3 ] = 0.014 + state [ 0 ] ; xx [ 4 ]
= 0.5 ; xx [ 5 ] = xx [ 4 ] * state [ 2 ] ; xx [ 6 ] = xx [ 0 ] * cos ( xx [
5 ] ) ; xx [ 7 ] = xx [ 0 ] * sin ( xx [ 5 ] ) ; xx [ 5 ] = xx [ 6 ] + xx [ 7
] ; xx [ 8 ] = xx [ 6 ] - xx [ 7 ] ; xx [ 6 ] = 1.0e-3 ; xx [ 7 ] = 2.0 ; xx
[ 9 ] = xx [ 8 ] * xx [ 6 ] ; xx [ 10 ] = 1.500000000000001e-3 + xx [ 6 ] -
xx [ 7 ] * xx [ 8 ] * xx [ 9 ] ; xx [ 11 ] = xx [ 7 ] * xx [ 9 ] * xx [ 5 ] ;
xx [ 9 ] = xx [ 4 ] * state [ 4 ] ; xx [ 12 ] = cos ( xx [ 9 ] ) ; xx [ 13 ]
= sin ( xx [ 9 ] ) ; xx [ 9 ] = xx [ 6 ] * xx [ 13 ] ; xx [ 14 ] = xx [ 6 ] +
xx [ 6 ] - xx [ 7 ] * xx [ 9 ] * xx [ 13 ] ; xx [ 15 ] = xx [ 7 ] * xx [ 12 ]
* xx [ 9 ] ; xx [ 9 ] = xx [ 4 ] * state [ 6 ] ; xx [ 16 ] = cos ( xx [ 9 ] )
; xx [ 17 ] = sin ( xx [ 9 ] ) ; xx [ 9 ] = xx [ 6 ] * xx [ 17 ] ; xx [ 18 ]
= xx [ 6 ] + xx [ 6 ] - xx [ 7 ] * xx [ 9 ] * xx [ 17 ] ; xx [ 19 ] = xx [ 7
] * xx [ 16 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 4 ] * state [ 8 ] ; xx [ 20 ] =
cos ( xx [ 9 ] ) ; xx [ 21 ] = sin ( xx [ 9 ] ) ; xx [ 9 ] = xx [ 6 ] * xx [
21 ] ; xx [ 22 ] = xx [ 7 ] * xx [ 9 ] * xx [ 21 ] - xx [ 6 ] ; xx [ 23 ] =
1.500000000000001e-3 + xx [ 7 ] * xx [ 20 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 4 ]
* state [ 10 ] ; xx [ 24 ] = cos ( xx [ 9 ] ) ; xx [ 25 ] = sin ( xx [ 9 ] )
; xx [ 9 ] = xx [ 6 ] * xx [ 25 ] ; xx [ 26 ] = xx [ 6 ] - ( xx [ 7 ] * xx [
9 ] * xx [ 25 ] - xx [ 6 ] ) ; xx [ 27 ] = xx [ 7 ] * xx [ 24 ] * xx [ 9 ] ;
xx [ 9 ] = xx [ 4 ] * state [ 12 ] ; xx [ 28 ] = cos ( xx [ 9 ] ) ; xx [ 29 ]
= sin ( xx [ 9 ] ) ; xx [ 9 ] = xx [ 6 ] * xx [ 29 ] ; xx [ 30 ] = xx [ 6 ] -
( xx [ 7 ] * xx [ 9 ] * xx [ 29 ] - xx [ 6 ] ) ; xx [ 31 ] = xx [ 7 ] * xx [
28 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 4 ] * state [ 14 ] ; xx [ 4 ] = cos ( xx [
9 ] ) ; xx [ 32 ] = sin ( xx [ 9 ] ) ; xx [ 9 ] = xx [ 6 ] * xx [ 32 ] ; xx [
33 ] = xx [ 6 ] - ( xx [ 7 ] * xx [ 9 ] * xx [ 32 ] - xx [ 6 ] ) ; xx [ 6 ] =
xx [ 7 ] * xx [ 4 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 20 ] * xx [ 24 ] - xx [ 21 ]
* xx [ 25 ] ; xx [ 34 ] = xx [ 20 ] * xx [ 25 ] + xx [ 24 ] * xx [ 21 ] ; xx
[ 35 ] = xx [ 27 ] * xx [ 21 ] ; xx [ 36 ] = xx [ 21 ] * xx [ 26 ] ; xx [ 37
] = xx [ 26 ] - xx [ 7 ] * ( xx [ 20 ] * xx [ 35 ] + xx [ 36 ] * xx [ 21 ] )
- xx [ 22 ] ; xx [ 38 ] = ( xx [ 20 ] * xx [ 36 ] - xx [ 35 ] * xx [ 21 ] ) *
xx [ 7 ] + xx [ 27 ] + xx [ 23 ] ; xx [ 35 ] = xx [ 34 ] * xx [ 29 ] - xx [
28 ] * xx [ 9 ] ; xx [ 36 ] = xx [ 29 ] * xx [ 9 ] + xx [ 34 ] * xx [ 28 ] ;
xx [ 39 ] = xx [ 34 ] * xx [ 31 ] ; xx [ 40 ] = xx [ 34 ] * xx [ 30 ] ; xx [
41 ] = xx [ 30 ] - xx [ 7 ] * ( xx [ 39 ] * xx [ 9 ] + xx [ 34 ] * xx [ 40 ]
) + xx [ 37 ] ; xx [ 42 ] = ( xx [ 40 ] * xx [ 9 ] - xx [ 34 ] * xx [ 39 ] )
* xx [ 7 ] + xx [ 31 ] + xx [ 38 ] ; xx [ 39 ] = xx [ 36 ] * xx [ 33 ] ; xx [
40 ] = xx [ 36 ] * xx [ 6 ] ; xx [ 43 ] = xx [ 0 ] * xx [ 5 ] ; xx [ 44 ] =
xx [ 8 ] * xx [ 0 ] ; xx [ 45 ] = xx [ 43 ] + xx [ 44 ] ; xx [ 46 ] = xx [ 43
] - xx [ 44 ] ; xx [ 43 ] = xx [ 0 ] * xx [ 0 ] * xx [ 10 ] ; xx [ 44 ] = xx
[ 0 ] * xx [ 0 ] * xx [ 11 ] ; xx [ 0 ] = xx [ 7 ] * ( xx [ 43 ] - xx [ 44 ]
) - xx [ 10 ] + xx [ 3 ] ; xx [ 47 ] = ( xx [ 43 ] + xx [ 44 ] ) * xx [ 7 ] -
xx [ 11 ] ; xx [ 43 ] = xx [ 46 ] * xx [ 13 ] - xx [ 12 ] * xx [ 45 ] ; xx [
44 ] = xx [ 13 ] * xx [ 45 ] + xx [ 46 ] * xx [ 12 ] ; xx [ 48 ] = xx [ 46 ]
* xx [ 15 ] ; xx [ 49 ] = xx [ 46 ] * xx [ 14 ] ; xx [ 50 ] = xx [ 7 ] * ( xx
[ 48 ] * xx [ 45 ] + xx [ 46 ] * xx [ 49 ] ) - xx [ 14 ] + xx [ 0 ] ; xx [ 51
] = xx [ 15 ] - ( xx [ 46 ] * xx [ 48 ] - xx [ 49 ] * xx [ 45 ] ) * xx [ 7 ]
+ xx [ 47 ] ; xx [ 48 ] = xx [ 44 ] * xx [ 19 ] ; xx [ 49 ] = xx [ 44 ] * xx
[ 18 ] ; xx [ 52 ] = xx [ 28 ] * xx [ 4 ] - xx [ 29 ] * xx [ 32 ] ; xx [ 53 ]
= xx [ 28 ] * xx [ 32 ] + xx [ 4 ] * xx [ 29 ] ; xx [ 54 ] = xx [ 6 ] * xx [
29 ] ; xx [ 55 ] = xx [ 29 ] * xx [ 33 ] ; xx [ 56 ] = xx [ 33 ] - xx [ 7 ] *
( xx [ 28 ] * xx [ 54 ] + xx [ 55 ] * xx [ 29 ] ) + xx [ 30 ] ; xx [ 57 ] = (
xx [ 28 ] * xx [ 55 ] - xx [ 54 ] * xx [ 29 ] ) * xx [ 7 ] + xx [ 6 ] + xx [
31 ] ; xx [ 54 ] = xx [ 25 ] * xx [ 57 ] ; xx [ 55 ] = xx [ 56 ] * xx [ 25 ]
; xx [ 58 ] = xx [ 12 ] * xx [ 16 ] - xx [ 13 ] * xx [ 17 ] ; xx [ 59 ] = xx
[ 12 ] * xx [ 17 ] + xx [ 16 ] * xx [ 13 ] ; xx [ 60 ] = xx [ 19 ] * xx [ 13
] ; xx [ 61 ] = xx [ 13 ] * xx [ 18 ] ; xx [ 62 ] = xx [ 7 ] * ( xx [ 12 ] *
xx [ 60 ] + xx [ 61 ] * xx [ 13 ] ) - xx [ 18 ] - xx [ 14 ] ; xx [ 63 ] = xx
[ 19 ] - ( xx [ 60 ] * xx [ 13 ] - xx [ 12 ] * xx [ 61 ] ) * xx [ 7 ] + xx [
15 ] ; xx [ 60 ] = xx [ 8 ] * xx [ 63 ] ; xx [ 61 ] = xx [ 62 ] * xx [ 8 ] ;
xx [ 64 ] = xx [ 8 ] * state [ 1 ] ; xx [ 65 ] = xx [ 7 ] * xx [ 64 ] * xx [
5 ] ; xx [ 66 ] = 1.0e-3 ; xx [ 67 ] = state [ 1 ] - xx [ 7 ] * xx [ 8 ] * xx
[ 64 ] + xx [ 66 ] * state [ 3 ] ; xx [ 64 ] = state [ 3 ] + state [ 5 ] ; xx
[ 68 ] = xx [ 15 ] * state [ 3 ] + xx [ 65 ] ; xx [ 69 ] = xx [ 67 ] + xx [
14 ] * state [ 3 ] ; xx [ 70 ] = xx [ 13 ] * xx [ 69 ] ; xx [ 71 ] = xx [ 13
] * xx [ 68 ] ; xx [ 72 ] = xx [ 68 ] - ( xx [ 12 ] * xx [ 70 ] + xx [ 71 ] *
xx [ 13 ] ) * xx [ 7 ] ; xx [ 68 ] = xx [ 69 ] + xx [ 7 ] * ( xx [ 12 ] * xx
[ 71 ] - xx [ 70 ] * xx [ 13 ] ) + xx [ 66 ] * state [ 5 ] ; xx [ 69 ] = xx [
64 ] * xx [ 19 ] + xx [ 72 ] ; xx [ 70 ] = xx [ 68 ] + xx [ 64 ] * xx [ 18 ]
; xx [ 71 ] = xx [ 17 ] * xx [ 70 ] ; xx [ 73 ] = xx [ 69 ] * xx [ 17 ] ; xx
[ 74 ] = xx [ 66 ] * state [ 9 ] ; xx [ 75 ] = state [ 9 ] + state [ 11 ] ;
xx [ 76 ] = xx [ 26 ] * state [ 9 ] + xx [ 74 ] ; xx [ 77 ] = xx [ 76 ] * xx
[ 25 ] ; xx [ 78 ] = xx [ 27 ] * state [ 9 ] ; xx [ 79 ] = xx [ 78 ] * xx [
25 ] ; xx [ 80 ] = xx [ 7 ] * ( xx [ 24 ] * xx [ 77 ] + xx [ 79 ] * xx [ 25 ]
) - xx [ 78 ] ; xx [ 78 ] = xx [ 76 ] - ( xx [ 77 ] * xx [ 25 ] - xx [ 24 ] *
xx [ 79 ] ) * xx [ 7 ] + xx [ 66 ] * state [ 11 ] ; xx [ 76 ] = xx [ 75 ] +
state [ 13 ] ; xx [ 77 ] = xx [ 80 ] - xx [ 75 ] * xx [ 31 ] ; xx [ 79 ] = xx
[ 75 ] * xx [ 30 ] + xx [ 78 ] ; xx [ 81 ] = xx [ 79 ] * xx [ 29 ] ; xx [ 82
] = xx [ 77 ] * xx [ 29 ] ; xx [ 83 ] = xx [ 77 ] + xx [ 7 ] * ( xx [ 28 ] *
xx [ 81 ] - xx [ 82 ] * xx [ 29 ] ) ; xx [ 77 ] = xx [ 79 ] - ( xx [ 28 ] *
xx [ 82 ] + xx [ 81 ] * xx [ 29 ] ) * xx [ 7 ] + xx [ 66 ] * state [ 13 ] ;
xx [ 79 ] = xx [ 83 ] - xx [ 76 ] * xx [ 6 ] ; xx [ 81 ] = xx [ 76 ] * xx [
33 ] + xx [ 77 ] ; xx [ 82 ] = xx [ 81 ] * xx [ 32 ] ; xx [ 84 ] = xx [ 79 ]
* xx [ 32 ] ; motionData [ 0 ] = xx [ 1 ] ; motionData [ 1 ] = xx [ 2 ] ;
motionData [ 2 ] = xx [ 1 ] ; motionData [ 3 ] = xx [ 2 ] ; motionData [ 4 ]
= xx [ 3 ] ; motionData [ 5 ] = xx [ 2 ] ; motionData [ 6 ] = xx [ 2 ] ;
motionData [ 7 ] = - xx [ 5 ] ; motionData [ 8 ] = xx [ 2 ] ; motionData [ 9
] = xx [ 8 ] ; motionData [ 10 ] = xx [ 2 ] ; motionData [ 11 ] = - xx [ 10 ]
; motionData [ 12 ] = xx [ 2 ] ; motionData [ 13 ] = - xx [ 11 ] ; motionData
[ 14 ] = - xx [ 12 ] ; motionData [ 15 ] = xx [ 2 ] ; motionData [ 16 ] = -
xx [ 13 ] ; motionData [ 17 ] = xx [ 2 ] ; motionData [ 18 ] = - xx [ 14 ] ;
motionData [ 19 ] = xx [ 2 ] ; motionData [ 20 ] = xx [ 15 ] ; motionData [
21 ] = - xx [ 16 ] ; motionData [ 22 ] = xx [ 2 ] ; motionData [ 23 ] = - xx
[ 17 ] ; motionData [ 24 ] = xx [ 2 ] ; motionData [ 25 ] = - xx [ 18 ] ;
motionData [ 26 ] = xx [ 2 ] ; motionData [ 27 ] = xx [ 19 ] ; motionData [
28 ] = - xx [ 20 ] ; motionData [ 29 ] = xx [ 2 ] ; motionData [ 30 ] = xx [
21 ] ; motionData [ 31 ] = xx [ 2 ] ; motionData [ 32 ] = - xx [ 22 ] ;
motionData [ 33 ] = xx [ 2 ] ; motionData [ 34 ] = xx [ 23 ] ; motionData [
35 ] = - xx [ 24 ] ; motionData [ 36 ] = xx [ 2 ] ; motionData [ 37 ] = xx [
25 ] ; motionData [ 38 ] = xx [ 2 ] ; motionData [ 39 ] = xx [ 26 ] ;
motionData [ 40 ] = xx [ 2 ] ; motionData [ 41 ] = xx [ 27 ] ; motionData [
42 ] = - xx [ 28 ] ; motionData [ 43 ] = xx [ 2 ] ; motionData [ 44 ] = xx [
29 ] ; motionData [ 45 ] = xx [ 2 ] ; motionData [ 46 ] = xx [ 30 ] ;
motionData [ 47 ] = xx [ 2 ] ; motionData [ 48 ] = xx [ 31 ] ; motionData [
49 ] = - xx [ 4 ] ; motionData [ 50 ] = xx [ 2 ] ; motionData [ 51 ] = xx [
32 ] ; motionData [ 52 ] = xx [ 2 ] ; motionData [ 53 ] = xx [ 33 ] ;
motionData [ 54 ] = xx [ 2 ] ; motionData [ 55 ] = xx [ 6 ] ; motionData [ 56
] = xx [ 9 ] ; motionData [ 57 ] = xx [ 2 ] ; motionData [ 58 ] = - xx [ 34 ]
; motionData [ 59 ] = xx [ 2 ] ; motionData [ 60 ] = xx [ 37 ] ; motionData [
61 ] = xx [ 2 ] ; motionData [ 62 ] = xx [ 38 ] ; motionData [ 63 ] = xx [ 35
] ; motionData [ 64 ] = xx [ 2 ] ; motionData [ 65 ] = xx [ 36 ] ; motionData
[ 66 ] = xx [ 2 ] ; motionData [ 67 ] = xx [ 41 ] ; motionData [ 68 ] = xx [
2 ] ; motionData [ 69 ] = xx [ 42 ] ; motionData [ 70 ] = - ( xx [ 4 ] * xx [
35 ] + xx [ 36 ] * xx [ 32 ] ) ; motionData [ 71 ] = xx [ 2 ] ; motionData [
72 ] = xx [ 32 ] * xx [ 35 ] - xx [ 36 ] * xx [ 4 ] ; motionData [ 73 ] = xx
[ 2 ] ; motionData [ 74 ] = xx [ 33 ] - ( xx [ 36 ] * xx [ 39 ] - xx [ 40 ] *
xx [ 35 ] ) * xx [ 7 ] + xx [ 41 ] ; motionData [ 75 ] = xx [ 2 ] ;
motionData [ 76 ] = xx [ 6 ] - xx [ 7 ] * ( xx [ 36 ] * xx [ 40 ] + xx [ 39 ]
* xx [ 35 ] ) + xx [ 42 ] ; motionData [ 77 ] = xx [ 45 ] ; motionData [ 78 ]
= xx [ 2 ] ; motionData [ 79 ] = xx [ 46 ] ; motionData [ 80 ] = xx [ 2 ] ;
motionData [ 81 ] = xx [ 0 ] ; motionData [ 82 ] = xx [ 2 ] ; motionData [ 83
] = xx [ 47 ] ; motionData [ 84 ] = xx [ 43 ] ; motionData [ 85 ] = xx [ 2 ]
; motionData [ 86 ] = - xx [ 44 ] ; motionData [ 87 ] = xx [ 2 ] ; motionData
[ 88 ] = xx [ 50 ] ; motionData [ 89 ] = xx [ 2 ] ; motionData [ 90 ] = xx [
51 ] ; motionData [ 91 ] = - ( xx [ 16 ] * xx [ 43 ] + xx [ 44 ] * xx [ 17 ]
) ; motionData [ 92 ] = xx [ 2 ] ; motionData [ 93 ] = xx [ 44 ] * xx [ 16 ]
- xx [ 17 ] * xx [ 43 ] ; motionData [ 94 ] = xx [ 2 ] ; motionData [ 95 ] =
xx [ 50 ] - ( xx [ 18 ] + ( xx [ 48 ] * xx [ 43 ] - xx [ 44 ] * xx [ 49 ] ) *
xx [ 7 ] ) ; motionData [ 96 ] = xx [ 2 ] ; motionData [ 97 ] = xx [ 19 ] -
xx [ 7 ] * ( xx [ 49 ] * xx [ 43 ] + xx [ 44 ] * xx [ 48 ] ) + xx [ 51 ] ;
motionData [ 98 ] = xx [ 52 ] ; motionData [ 99 ] = xx [ 2 ] ; motionData [
100 ] = - xx [ 53 ] ; motionData [ 101 ] = xx [ 2 ] ; motionData [ 102 ] = xx
[ 56 ] ; motionData [ 103 ] = xx [ 2 ] ; motionData [ 104 ] = xx [ 57 ] ;
motionData [ 105 ] = xx [ 53 ] * xx [ 25 ] - xx [ 24 ] * xx [ 52 ] ;
motionData [ 106 ] = xx [ 2 ] ; motionData [ 107 ] = xx [ 53 ] * xx [ 24 ] +
xx [ 25 ] * xx [ 52 ] ; motionData [ 108 ] = xx [ 2 ] ; motionData [ 109 ] =
xx [ 56 ] - ( xx [ 24 ] * xx [ 54 ] + xx [ 55 ] * xx [ 25 ] ) * xx [ 7 ] + xx
[ 26 ] ; motionData [ 110 ] = xx [ 2 ] ; motionData [ 111 ] = xx [ 57 ] + xx
[ 7 ] * ( xx [ 24 ] * xx [ 55 ] - xx [ 54 ] * xx [ 25 ] ) + xx [ 27 ] ;
motionData [ 112 ] = xx [ 58 ] ; motionData [ 113 ] = xx [ 2 ] ; motionData [
114 ] = xx [ 59 ] ; motionData [ 115 ] = xx [ 2 ] ; motionData [ 116 ] = xx [
62 ] ; motionData [ 117 ] = xx [ 2 ] ; motionData [ 118 ] = xx [ 63 ] ;
motionData [ 119 ] = - ( xx [ 58 ] * xx [ 5 ] + xx [ 59 ] * xx [ 8 ] ) ;
motionData [ 120 ] = xx [ 2 ] ; motionData [ 121 ] = xx [ 8 ] * xx [ 58 ] -
xx [ 59 ] * xx [ 5 ] ; motionData [ 122 ] = xx [ 2 ] ; motionData [ 123 ] =
xx [ 62 ] - xx [ 7 ] * ( xx [ 60 ] * xx [ 5 ] + xx [ 8 ] * xx [ 61 ] ) - xx [
10 ] ; motionData [ 124 ] = xx [ 2 ] ; motionData [ 125 ] = xx [ 63 ] - ( xx
[ 8 ] * xx [ 60 ] - xx [ 61 ] * xx [ 5 ] ) * xx [ 7 ] - xx [ 11 ] ;
motionData [ 126 ] = xx [ 2 ] ; motionData [ 127 ] = xx [ 2 ] ; motionData [
128 ] = xx [ 2 ] ; motionData [ 129 ] = xx [ 2 ] ; motionData [ 130 ] = xx [
2 ] ; motionData [ 131 ] = state [ 1 ] ; motionData [ 132 ] = xx [ 2 ] ;
motionData [ 133 ] = state [ 3 ] ; motionData [ 134 ] = xx [ 2 ] ; motionData
[ 135 ] = xx [ 65 ] ; motionData [ 136 ] = xx [ 2 ] ; motionData [ 137 ] = xx
[ 67 ] ; motionData [ 138 ] = xx [ 2 ] ; motionData [ 139 ] = xx [ 64 ] ;
motionData [ 140 ] = xx [ 2 ] ; motionData [ 141 ] = xx [ 72 ] ; motionData [
142 ] = xx [ 2 ] ; motionData [ 143 ] = xx [ 68 ] ; motionData [ 144 ] = xx [
2 ] ; motionData [ 145 ] = xx [ 64 ] + state [ 7 ] ; motionData [ 146 ] = xx
[ 2 ] ; motionData [ 147 ] = xx [ 69 ] - ( xx [ 16 ] * xx [ 71 ] + xx [ 73 ]
* xx [ 17 ] ) * xx [ 7 ] ; motionData [ 148 ] = xx [ 2 ] ; motionData [ 149 ]
= xx [ 70 ] + xx [ 7 ] * ( xx [ 16 ] * xx [ 73 ] - xx [ 71 ] * xx [ 17 ] ) +
xx [ 66 ] * state [ 7 ] ; motionData [ 150 ] = xx [ 2 ] ; motionData [ 151 ]
= - state [ 9 ] ; motionData [ 152 ] = xx [ 2 ] ; motionData [ 153 ] = xx [ 2
] ; motionData [ 154 ] = xx [ 2 ] ; motionData [ 155 ] = xx [ 74 ] ;
motionData [ 156 ] = xx [ 2 ] ; motionData [ 157 ] = - xx [ 75 ] ; motionData
[ 158 ] = xx [ 2 ] ; motionData [ 159 ] = xx [ 80 ] ; motionData [ 160 ] = xx
[ 2 ] ; motionData [ 161 ] = xx [ 78 ] ; motionData [ 162 ] = xx [ 2 ] ;
motionData [ 163 ] = - xx [ 76 ] ; motionData [ 164 ] = xx [ 2 ] ; motionData
[ 165 ] = xx [ 83 ] ; motionData [ 166 ] = xx [ 2 ] ; motionData [ 167 ] = xx
[ 77 ] ; motionData [ 168 ] = xx [ 2 ] ; motionData [ 169 ] = - ( xx [ 76 ] +
state [ 15 ] ) ; motionData [ 170 ] = xx [ 2 ] ; motionData [ 171 ] = xx [ 79
] + xx [ 7 ] * ( xx [ 4 ] * xx [ 82 ] - xx [ 84 ] * xx [ 32 ] ) ; motionData
[ 172 ] = xx [ 2 ] ; motionData [ 173 ] = xx [ 81 ] - ( xx [ 4 ] * xx [ 84 ]
+ xx [ 82 ] * xx [ 32 ] ) * xx [ 7 ] + xx [ 66 ] * state [ 15 ] ; } static
size_t computeAssemblyError_0 ( const RuntimeDerivedValuesBundle * rtdv ,
const int * modeVector , const double * motionData , double * error ) { const
double * rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv ->
mInts . mValues ; double xx [ 13 ] ; ( void ) rtdvd ; ( void ) rtdvi ; ( void
) modeVector ; xx [ 0 ] = 0.7071067811865476 ; xx [ 1 ] = xx [ 0 ] *
motionData [ 70 ] ; xx [ 2 ] = xx [ 0 ] * motionData [ 71 ] ; xx [ 3 ] = xx [
0 ] * motionData [ 72 ] ; xx [ 4 ] = xx [ 0 ] * motionData [ 73 ] ; xx [ 5 ]
= xx [ 1 ] - xx [ 2 ] ; xx [ 6 ] = xx [ 1 ] + xx [ 2 ] ; xx [ 7 ] = xx [ 3 ]
+ xx [ 4 ] ; xx [ 8 ] = xx [ 4 ] - xx [ 3 ] ; xx [ 1 ] = xx [ 0 ] *
motionData [ 91 ] ; xx [ 2 ] = xx [ 0 ] * motionData [ 92 ] ; xx [ 3 ] = xx [
0 ] * motionData [ 93 ] ; xx [ 4 ] = xx [ 0 ] * motionData [ 94 ] ; xx [ 9 ]
= xx [ 1 ] - xx [ 2 ] ; xx [ 10 ] = xx [ 1 ] + xx [ 2 ] ; xx [ 11 ] = xx [ 3
] + xx [ 4 ] ; xx [ 12 ] = xx [ 4 ] - xx [ 3 ] ;
pm_math_Quaternion_inverseCompose_ra ( xx + 5 , xx + 9 , xx + 0 ) ; xx [ 0 ]
= 1.0e-3 ; xx [ 3 ] = xx [ 0 ] * motionData [ 93 ] ; xx [ 4 ] = xx [ 0 ] *
motionData [ 94 ] ; xx [ 5 ] = 2.0 ; xx [ 6 ] = xx [ 0 ] * motionData [ 72 ]
; xx [ 7 ] = xx [ 0 ] * motionData [ 73 ] ; error [ 0 ] = xx [ 1 ] ; error [
1 ] = xx [ 2 ] ; error [ 2 ] = motionData [ 95 ] + ( xx [ 3 ] * motionData [
93 ] + xx [ 4 ] * motionData [ 94 ] ) * xx [ 5 ] - ( motionData [ 74 ] - ( xx
[ 6 ] * motionData [ 72 ] + xx [ 7 ] * motionData [ 73 ] ) * xx [ 5 ] ) -
2.0e-3 ; error [ 3 ] = motionData [ 96 ] - ( xx [ 4 ] * motionData [ 91 ] +
xx [ 3 ] * motionData [ 92 ] ) * xx [ 5 ] - ( ( xx [ 7 ] * motionData [ 70 ]
+ xx [ 6 ] * motionData [ 71 ] ) * xx [ 5 ] + motionData [ 75 ] ) ; error [ 4
] = xx [ 5 ] * ( xx [ 3 ] * motionData [ 91 ] - xx [ 4 ] * motionData [ 92 ]
) + motionData [ 97 ] - ( xx [ 5 ] * ( xx [ 7 ] * motionData [ 71 ] - xx [ 6
] * motionData [ 70 ] ) + motionData [ 76 ] ) ; return 5 ; } size_t
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computeAssemblyError ( const
void * mech , const RuntimeDerivedValuesBundle * rtdv , size_t constraintIdx
, const int * modeVector , const double * motionData , double * error ) { (
void ) mech ; ( void ) rtdv ; ( void ) modeVector ; ( void ) motionData ; (
void ) error ; switch ( constraintIdx ) { case 0 : return
computeAssemblyError_0 ( rtdv , modeVector , motionData , error ) ; } return
0 ; } static size_t computeAssemblyJacobian_0 ( const
RuntimeDerivedValuesBundle * rtdv , const double * state , const int *
modeVector , const double * motionData , double * J ) { const double * rtdvd
= rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts . mValues ;
double xx [ 137 ] ; ( void ) rtdvd ; ( void ) rtdvi ; ( void ) modeVector ;
xx [ 0 ] = 0.7071067811865476 ; xx [ 1 ] = xx [ 0 ] * motionData [ 70 ] ; xx
[ 2 ] = xx [ 0 ] * motionData [ 71 ] ; xx [ 3 ] = xx [ 0 ] * motionData [ 72
] ; xx [ 4 ] = xx [ 0 ] * motionData [ 73 ] ; xx [ 5 ] = xx [ 1 ] - xx [ 2 ]
; xx [ 6 ] = xx [ 1 ] + xx [ 2 ] ; xx [ 7 ] = xx [ 3 ] + xx [ 4 ] ; xx [ 8 ]
= xx [ 4 ] - xx [ 3 ] ; xx [ 1 ] = xx [ 0 ] * motionData [ 91 ] ; xx [ 2 ] =
xx [ 0 ] * motionData [ 92 ] ; xx [ 3 ] = xx [ 0 ] * motionData [ 93 ] ; xx [
4 ] = xx [ 0 ] * motionData [ 94 ] ; xx [ 9 ] = xx [ 1 ] - xx [ 2 ] ; xx [ 10
] = xx [ 1 ] + xx [ 2 ] ; xx [ 11 ] = xx [ 3 ] + xx [ 4 ] ; xx [ 12 ] = xx [
4 ] - xx [ 3 ] ; pm_math_Quaternion_inverseCompose_ra ( xx + 5 , xx + 9 , xx
+ 1 ) ; xx [ 5 ] = 0.0 ; xx [ 6 ] = xx [ 5 ] ; xx [ 7 ] = xx [ 5 ] ; xx [ 8 ]
= xx [ 5 ] ; pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 6 , xx + 13 ) ;
xx [ 6 ] = 2.0 ; xx [ 7 ] = 1.0 ; xx [ 8 ] = xx [ 7 ] - ( motionData [ 115 ]
* motionData [ 115 ] + motionData [ 113 ] * motionData [ 113 ] ) * xx [ 6 ] ;
xx [ 17 ] = xx [ 0 ] * xx [ 0 ] * xx [ 8 ] ; xx [ 18 ] = xx [ 6 ] * (
motionData [ 114 ] * motionData [ 115 ] - motionData [ 112 ] * motionData [
113 ] ) ; xx [ 19 ] = xx [ 0 ] * xx [ 0 ] * xx [ 18 ] ; xx [ 20 ] = (
motionData [ 112 ] * motionData [ 115 ] + motionData [ 113 ] * motionData [
114 ] ) * xx [ 6 ] ; xx [ 21 ] = xx [ 8 ] - ( xx [ 17 ] - xx [ 19 ] ) * xx [
6 ] ; xx [ 22 ] = xx [ 18 ] - xx [ 6 ] * ( xx [ 17 ] + xx [ 19 ] ) ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 20 , xx + 23 ) ; xx [ 17 ] =
xx [ 7 ] - ( motionData [ 24 ] * motionData [ 24 ] + motionData [ 22 ] *
motionData [ 22 ] ) * xx [ 6 ] ; xx [ 19 ] = xx [ 0 ] * xx [ 0 ] * xx [ 17 ]
; xx [ 20 ] = xx [ 6 ] * ( motionData [ 23 ] * motionData [ 24 ] - motionData
[ 21 ] * motionData [ 22 ] ) ; xx [ 21 ] = xx [ 0 ] * xx [ 0 ] * xx [ 20 ] ;
xx [ 27 ] = ( motionData [ 21 ] * motionData [ 24 ] + motionData [ 22 ] *
motionData [ 23 ] ) * xx [ 6 ] ; xx [ 28 ] = xx [ 17 ] - ( xx [ 19 ] - xx [
21 ] ) * xx [ 6 ] ; xx [ 29 ] = xx [ 20 ] - xx [ 6 ] * ( xx [ 19 ] + xx [ 21
] ) ; pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 27 , xx + 30 ) ; xx [
27 ] = xx [ 5 ] ; xx [ 28 ] = - 2.220446049250313e-16 ; xx [ 29 ] = - 1.0 ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 27 , xx + 34 ) ; xx [ 38 ] =
motionData [ 70 ] ; xx [ 39 ] = motionData [ 71 ] ; xx [ 40 ] = motionData [
72 ] ; xx [ 41 ] = motionData [ 73 ] ; xx [ 19 ] = ( motionData [ 108 ] *
motionData [ 108 ] + motionData [ 106 ] * motionData [ 106 ] ) * xx [ 6 ] -
xx [ 7 ] ; xx [ 21 ] = motionData [ 105 ] * motionData [ 106 ] - motionData [
107 ] * motionData [ 108 ] ; xx [ 27 ] = - ( ( motionData [ 105 ] *
motionData [ 108 ] + motionData [ 106 ] * motionData [ 107 ] ) * xx [ 6 ] ) ;
xx [ 28 ] = xx [ 19 ] ; xx [ 29 ] = xx [ 6 ] * xx [ 21 ] ;
pm_math_Quaternion_xform_ra ( xx + 38 , xx + 27 , xx + 42 ) ;
pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 42 , xx + 27 ) ; xx [ 42 ]
= - xx [ 27 ] ; xx [ 43 ] = - xx [ 28 ] ; xx [ 44 ] = - xx [ 29 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 42 , xx + 45 ) ; xx [ 22 ] =
( motionData [ 101 ] * motionData [ 101 ] + motionData [ 99 ] * motionData [
99 ] ) * xx [ 6 ] - xx [ 7 ] ; xx [ 27 ] = motionData [ 98 ] * motionData [
99 ] - motionData [ 100 ] * motionData [ 101 ] ; xx [ 42 ] = - ( ( motionData
[ 98 ] * motionData [ 101 ] + motionData [ 99 ] * motionData [ 100 ] ) * xx [
6 ] ) ; xx [ 43 ] = xx [ 22 ] ; xx [ 44 ] = xx [ 6 ] * xx [ 27 ] ;
pm_math_Quaternion_xform_ra ( xx + 38 , xx + 42 , xx + 49 ) ;
pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 49 , xx + 42 ) ; xx [ 49 ]
= - xx [ 42 ] ; xx [ 50 ] = - xx [ 43 ] ; xx [ 51 ] = - xx [ 44 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 49 , xx + 52 ) ; xx [ 28 ] =
( motionData [ 52 ] * motionData [ 52 ] + motionData [ 50 ] * motionData [ 50
] ) * xx [ 6 ] - xx [ 7 ] ; xx [ 29 ] = motionData [ 49 ] * motionData [ 50 ]
- motionData [ 51 ] * motionData [ 52 ] ; xx [ 42 ] = - ( ( motionData [ 49 ]
* motionData [ 52 ] + motionData [ 50 ] * motionData [ 51 ] ) * xx [ 6 ] ) ;
xx [ 43 ] = xx [ 28 ] ; xx [ 44 ] = xx [ 6 ] * xx [ 29 ] ;
pm_math_Quaternion_xform_ra ( xx + 38 , xx + 42 , xx + 49 ) ;
pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 49 , xx + 38 ) ; xx [ 41 ]
= - xx [ 38 ] ; xx [ 42 ] = - xx [ 39 ] ; xx [ 43 ] = - xx [ 40 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 41 , xx + 56 ) ; xx [ 38 ] =
xx [ 6 ] * ( motionData [ 70 ] * motionData [ 73 ] - motionData [ 71 ] *
motionData [ 72 ] ) ; xx [ 39 ] = ( motionData [ 73 ] * motionData [ 73 ] +
motionData [ 71 ] * motionData [ 71 ] ) * xx [ 6 ] - xx [ 7 ] ; xx [ 40 ] = -
( ( motionData [ 70 ] * motionData [ 71 ] + motionData [ 72 ] * motionData [
73 ] ) * xx [ 6 ] ) ; pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 38 ,
xx + 41 ) ; xx [ 9 ] = - xx [ 41 ] ; xx [ 10 ] = - xx [ 42 ] ; xx [ 11 ] = -
xx [ 43 ] ; pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 9 , xx + 38 ) ;
xx [ 1 ] = 0.5 ; xx [ 2 ] = xx [ 1 ] * state [ 2 ] ; xx [ 3 ] = xx [ 0 ] *
cos ( xx [ 2 ] ) ; xx [ 4 ] = xx [ 0 ] * sin ( xx [ 2 ] ) ; xx [ 0 ] = xx [ 3
] + xx [ 4 ] ; xx [ 2 ] = xx [ 3 ] - xx [ 4 ] ; xx [ 9 ] = - ( xx [ 0 ] *
motionData [ 0 ] + xx [ 2 ] * motionData [ 2 ] ) ; xx [ 10 ] = - ( xx [ 0 ] *
motionData [ 1 ] + xx [ 2 ] * motionData [ 3 ] ) ; xx [ 11 ] = xx [ 2 ] *
motionData [ 0 ] - xx [ 0 ] * motionData [ 2 ] ; xx [ 12 ] = xx [ 2 ] *
motionData [ 1 ] - xx [ 0 ] * motionData [ 3 ] ; xx [ 41 ] = motionData [ 112
] ; xx [ 42 ] = motionData [ 113 ] ; xx [ 43 ] = motionData [ 114 ] ; xx [ 44
] = motionData [ 115 ] ; pm_math_Quaternion_compose_ra ( xx + 9 , xx + 41 ,
xx + 48 ) ; xx [ 3 ] = 1.0e-3 ; xx [ 4 ] = xx [ 3 ] * xx [ 8 ] ; xx [ 8 ] =
xx [ 3 ] * xx [ 18 ] ; xx [ 9 ] = xx [ 50 ] * xx [ 4 ] + xx [ 51 ] * xx [ 8 ]
; xx [ 10 ] = xx [ 49 ] * xx [ 4 ] ; xx [ 11 ] = xx [ 49 ] * xx [ 8 ] ; xx [
41 ] = xx [ 9 ] ; xx [ 42 ] = - xx [ 10 ] ; xx [ 43 ] = - xx [ 11 ] ;
pm_math_Vector3_cross_ra ( xx + 49 , xx + 41 , xx + 59 ) ; xx [ 12 ] = xx [ 2
] * motionData [ 118 ] ; xx [ 13 ] = xx [ 2 ] * motionData [ 116 ] ; xx [ 16
] = 1.0e-3 ; xx [ 18 ] = xx [ 2 ] * xx [ 16 ] ; xx [ 23 ] = motionData [ 118
] - ( xx [ 2 ] * xx [ 12 ] - xx [ 13 ] * xx [ 0 ] ) * xx [ 6 ] - xx [ 6 ] *
xx [ 18 ] * xx [ 0 ] ; xx [ 26 ] = xx [ 6 ] * ( xx [ 2 ] * xx [ 13 ] + xx [
12 ] * xx [ 0 ] ) - ( motionData [ 116 ] + xx [ 6 ] * xx [ 2 ] * xx [ 18 ] )
+ xx [ 16 ] ; xx [ 0 ] = xx [ 26 ] * motionData [ 2 ] ; xx [ 41 ] =
motionData [ 1 ] ; xx [ 42 ] = motionData [ 2 ] ; xx [ 43 ] = motionData [ 3
] ; xx [ 2 ] = xx [ 23 ] * motionData [ 3 ] - xx [ 26 ] * motionData [ 1 ] ;
xx [ 12 ] = xx [ 23 ] * motionData [ 2 ] ; xx [ 49 ] = xx [ 0 ] ; xx [ 50 ] =
xx [ 2 ] ; xx [ 51 ] = - xx [ 12 ] ; pm_math_Vector3_cross_ra ( xx + 41 , xx
+ 49 , xx + 62 ) ; xx [ 13 ] = xx [ 1 ] * state [ 4 ] ; xx [ 18 ] = sin ( xx
[ 13 ] ) ; xx [ 30 ] = cos ( xx [ 13 ] ) ; xx [ 41 ] = xx [ 18 ] * motionData
[ 79 ] - xx [ 30 ] * motionData [ 77 ] ; xx [ 42 ] = xx [ 18 ] * motionData [
80 ] - xx [ 30 ] * motionData [ 78 ] ; xx [ 43 ] = - ( xx [ 18 ] * motionData
[ 77 ] + xx [ 30 ] * motionData [ 79 ] ) ; xx [ 44 ] = - ( xx [ 30 ] *
motionData [ 80 ] + xx [ 18 ] * motionData [ 78 ] ) ; xx [ 49 ] = motionData
[ 21 ] ; xx [ 50 ] = motionData [ 22 ] ; xx [ 51 ] = motionData [ 23 ] ; xx [
52 ] = motionData [ 24 ] ; pm_math_Quaternion_compose_ra ( xx + 41 , xx + 49
, xx + 65 ) ; xx [ 13 ] = xx [ 3 ] * xx [ 17 ] ; xx [ 17 ] = xx [ 3 ] * xx [
20 ] ; xx [ 20 ] = xx [ 67 ] * xx [ 13 ] + xx [ 68 ] * xx [ 17 ] ; xx [ 33 ]
= xx [ 66 ] * xx [ 13 ] ; xx [ 34 ] = xx [ 66 ] * xx [ 17 ] ; xx [ 41 ] = xx
[ 20 ] ; xx [ 42 ] = - xx [ 33 ] ; xx [ 43 ] = - xx [ 34 ] ;
pm_math_Vector3_cross_ra ( xx + 66 , xx + 41 , xx + 49 ) ; xx [ 37 ] = xx [
18 ] * motionData [ 25 ] ; xx [ 38 ] = xx [ 18 ] * motionData [ 27 ] ; xx [
41 ] = xx [ 16 ] * xx [ 18 ] ; xx [ 42 ] = motionData [ 27 ] - ( xx [ 30 ] *
xx [ 37 ] + xx [ 38 ] * xx [ 18 ] ) * xx [ 6 ] + xx [ 6 ] * xx [ 30 ] * xx [
41 ] ; xx [ 43 ] = xx [ 6 ] * ( xx [ 37 ] * xx [ 18 ] - xx [ 30 ] * xx [ 38 ]
) - ( motionData [ 25 ] + xx [ 6 ] * xx [ 41 ] * xx [ 18 ] ) + xx [ 16 ] ; xx
[ 18 ] = xx [ 43 ] * motionData [ 79 ] ; xx [ 66 ] = motionData [ 78 ] ; xx [
67 ] = motionData [ 79 ] ; xx [ 68 ] = motionData [ 80 ] ; xx [ 30 ] = xx [
42 ] * motionData [ 80 ] - xx [ 43 ] * motionData [ 78 ] ; xx [ 37 ] = xx [
42 ] * motionData [ 79 ] ; xx [ 69 ] = xx [ 18 ] ; xx [ 70 ] = xx [ 30 ] ; xx
[ 71 ] = - xx [ 37 ] ; pm_math_Vector3_cross_ra ( xx + 66 , xx + 69 , xx + 72
) ; xx [ 38 ] = xx [ 1 ] * state [ 6 ] ; xx [ 41 ] = cos ( xx [ 38 ] ) ; xx [
44 ] = sin ( xx [ 38 ] ) ; xx [ 38 ] = xx [ 16 ] * xx [ 44 ] ; xx [ 45 ] = xx
[ 6 ] * xx [ 41 ] * xx [ 38 ] ; xx [ 52 ] = xx [ 16 ] - xx [ 6 ] * xx [ 38 ]
* xx [ 44 ] ; xx [ 38 ] = xx [ 52 ] * motionData [ 86 ] ; xx [ 66 ] =
motionData [ 85 ] ; xx [ 67 ] = motionData [ 86 ] ; xx [ 68 ] = motionData [
87 ] ; xx [ 55 ] = xx [ 45 ] * motionData [ 87 ] - xx [ 52 ] * motionData [
85 ] ; xx [ 56 ] = xx [ 45 ] * motionData [ 86 ] ; xx [ 69 ] = xx [ 38 ] ; xx
[ 70 ] = xx [ 55 ] ; xx [ 71 ] = - xx [ 56 ] ; pm_math_Vector3_cross_ra ( xx
+ 66 , xx + 69 , xx + 75 ) ; xx [ 66 ] = xx [ 44 ] * motionData [ 84 ] + xx [
41 ] * motionData [ 86 ] ; xx [ 67 ] = xx [ 66 ] * xx [ 3 ] ; xx [ 68 ] = xx
[ 44 ] * motionData [ 86 ] - xx [ 41 ] * motionData [ 84 ] ; xx [ 69 ] = xx [
41 ] * motionData [ 87 ] + xx [ 44 ] * motionData [ 85 ] ; xx [ 70 ] = xx [
44 ] * motionData [ 87 ] - xx [ 41 ] * motionData [ 85 ] ; xx [ 41 ] = xx [ 3
] * xx [ 70 ] ; xx [ 44 ] = xx [ 1 ] * state [ 8 ] ; xx [ 71 ] = sin ( xx [
44 ] ) ; xx [ 78 ] = cos ( xx [ 44 ] ) ; xx [ 44 ] = xx [ 71 ] * motionData [
108 ] - xx [ 78 ] * motionData [ 106 ] ; xx [ 79 ] = xx [ 71 ] * motionData [
105 ] - xx [ 78 ] * motionData [ 107 ] ; xx [ 80 ] = xx [ 78 ] * motionData [
108 ] + xx [ 71 ] * motionData [ 106 ] ; xx [ 81 ] = xx [ 44 ] ; xx [ 82 ] =
xx [ 79 ] ; xx [ 83 ] = - xx [ 80 ] ; xx [ 84 ] = 2.0e-3 ; xx [ 85 ] = xx [
84 ] * xx [ 21 ] ; xx [ 21 ] = xx [ 19 ] * xx [ 3 ] ; xx [ 19 ] = xx [ 80 ] *
xx [ 85 ] - xx [ 21 ] * xx [ 79 ] ; xx [ 79 ] = xx [ 21 ] * xx [ 44 ] ; xx [
80 ] = xx [ 85 ] * xx [ 44 ] ; xx [ 86 ] = xx [ 19 ] ; xx [ 87 ] = xx [ 79 ]
; xx [ 88 ] = xx [ 80 ] ; pm_math_Vector3_cross_ra ( xx + 81 , xx + 86 , xx +
89 ) ; xx [ 44 ] = xx [ 78 ] * motionData [ 105 ] + xx [ 71 ] * motionData [
107 ] ; xx [ 81 ] = xx [ 71 ] * motionData [ 111 ] ; xx [ 82 ] = xx [ 71 ] *
motionData [ 109 ] ; xx [ 83 ] = xx [ 16 ] * xx [ 71 ] ; xx [ 86 ] = xx [ 1 ]
* state [ 10 ] ; xx [ 87 ] = cos ( xx [ 86 ] ) ; xx [ 88 ] = sin ( xx [ 86 ]
) ; xx [ 92 ] = - ( xx [ 87 ] * motionData [ 28 ] + xx [ 88 ] * motionData [
30 ] ) ; xx [ 93 ] = - ( xx [ 87 ] * motionData [ 29 ] + xx [ 88 ] *
motionData [ 31 ] ) ; xx [ 94 ] = xx [ 88 ] * motionData [ 28 ] - xx [ 87 ] *
motionData [ 30 ] ; xx [ 95 ] = xx [ 88 ] * motionData [ 29 ] - xx [ 87 ] *
motionData [ 31 ] ; xx [ 96 ] = motionData [ 98 ] ; xx [ 97 ] = motionData [
99 ] ; xx [ 98 ] = motionData [ 100 ] ; xx [ 99 ] = motionData [ 101 ] ;
pm_math_Quaternion_compose_ra ( xx + 92 , xx + 96 , xx + 100 ) ; xx [ 86 ] =
xx [ 22 ] * xx [ 3 ] ; xx [ 22 ] = xx [ 84 ] * xx [ 27 ] ; xx [ 27 ] = xx [
102 ] * xx [ 86 ] + xx [ 103 ] * xx [ 22 ] ; xx [ 92 ] = xx [ 101 ] * xx [ 86
] ; xx [ 93 ] = xx [ 101 ] * xx [ 22 ] ; xx [ 94 ] = - xx [ 27 ] ; xx [ 95 ]
= xx [ 92 ] ; xx [ 96 ] = xx [ 93 ] ; pm_math_Vector3_cross_ra ( xx + 101 ,
xx + 94 , xx + 97 ) ; xx [ 94 ] = xx [ 88 ] * motionData [ 104 ] ; xx [ 95 ]
= xx [ 88 ] * motionData [ 102 ] ; xx [ 96 ] = xx [ 16 ] * xx [ 88 ] ; xx [
101 ] = xx [ 6 ] * ( xx [ 94 ] * xx [ 88 ] - xx [ 87 ] * xx [ 95 ] ) -
motionData [ 104 ] - xx [ 6 ] * xx [ 87 ] * xx [ 96 ] ; xx [ 102 ] =
motionData [ 102 ] - ( ( xx [ 87 ] * xx [ 94 ] + xx [ 95 ] * xx [ 88 ] ) * xx
[ 6 ] + xx [ 6 ] * xx [ 96 ] * xx [ 88 ] ) + xx [ 16 ] ; xx [ 87 ] = xx [ 102
] * motionData [ 30 ] ; xx [ 94 ] = motionData [ 29 ] ; xx [ 95 ] =
motionData [ 30 ] ; xx [ 96 ] = motionData [ 31 ] ; xx [ 88 ] = xx [ 101 ] *
motionData [ 31 ] - xx [ 102 ] * motionData [ 29 ] ; xx [ 103 ] = xx [ 101 ]
* motionData [ 30 ] ; xx [ 104 ] = xx [ 87 ] ; xx [ 105 ] = xx [ 88 ] ; xx [
106 ] = - xx [ 103 ] ; pm_math_Vector3_cross_ra ( xx + 94 , xx + 104 , xx +
107 ) ; xx [ 94 ] = xx [ 1 ] * state [ 12 ] ; xx [ 95 ] = cos ( xx [ 94 ] ) ;
xx [ 96 ] = sin ( xx [ 94 ] ) ; xx [ 110 ] = - ( xx [ 95 ] * motionData [ 56
] + xx [ 96 ] * motionData [ 58 ] ) ; xx [ 111 ] = - ( xx [ 95 ] * motionData
[ 57 ] + xx [ 96 ] * motionData [ 59 ] ) ; xx [ 112 ] = xx [ 96 ] *
motionData [ 56 ] - xx [ 95 ] * motionData [ 58 ] ; xx [ 113 ] = xx [ 96 ] *
motionData [ 57 ] - xx [ 95 ] * motionData [ 59 ] ; xx [ 114 ] = motionData [
49 ] ; xx [ 115 ] = motionData [ 50 ] ; xx [ 116 ] = motionData [ 51 ] ; xx [
117 ] = motionData [ 52 ] ; pm_math_Quaternion_compose_ra ( xx + 110 , xx +
114 , xx + 118 ) ; xx [ 94 ] = xx [ 28 ] * xx [ 3 ] ; xx [ 28 ] = xx [ 84 ] *
xx [ 29 ] ; xx [ 29 ] = xx [ 120 ] * xx [ 94 ] + xx [ 121 ] * xx [ 28 ] ; xx
[ 84 ] = xx [ 119 ] * xx [ 94 ] ; xx [ 104 ] = xx [ 119 ] * xx [ 28 ] ; xx [
110 ] = - xx [ 29 ] ; xx [ 111 ] = xx [ 84 ] ; xx [ 112 ] = xx [ 104 ] ;
pm_math_Vector3_cross_ra ( xx + 119 , xx + 110 , xx + 113 ) ; xx [ 105 ] = xx
[ 96 ] * motionData [ 55 ] ; xx [ 106 ] = xx [ 96 ] * motionData [ 53 ] ; xx
[ 110 ] = xx [ 16 ] * xx [ 96 ] ; xx [ 111 ] = xx [ 6 ] * ( xx [ 105 ] * xx [
96 ] - xx [ 95 ] * xx [ 106 ] ) - motionData [ 55 ] - xx [ 6 ] * xx [ 95 ] *
xx [ 110 ] ; xx [ 112 ] = motionData [ 53 ] - ( ( xx [ 95 ] * xx [ 105 ] + xx
[ 106 ] * xx [ 96 ] ) * xx [ 6 ] + xx [ 6 ] * xx [ 110 ] * xx [ 96 ] ) + xx [
16 ] ; xx [ 95 ] = xx [ 112 ] * motionData [ 58 ] ; xx [ 119 ] = motionData [
57 ] ; xx [ 120 ] = motionData [ 58 ] ; xx [ 121 ] = motionData [ 59 ] ; xx [
96 ] = xx [ 111 ] * motionData [ 59 ] - xx [ 112 ] * motionData [ 57 ] ; xx [
105 ] = xx [ 111 ] * motionData [ 58 ] ; xx [ 122 ] = xx [ 95 ] ; xx [ 123 ]
= xx [ 96 ] ; xx [ 124 ] = - xx [ 105 ] ; pm_math_Vector3_cross_ra ( xx + 119
, xx + 122 , xx + 125 ) ; xx [ 106 ] = xx [ 1 ] * state [ 14 ] ; xx [ 1 ] =
cos ( xx [ 106 ] ) ; xx [ 110 ] = sin ( xx [ 106 ] ) ; xx [ 106 ] = xx [ 1 ]
* motionData [ 63 ] + xx [ 110 ] * motionData [ 65 ] ; xx [ 116 ] = xx [ 110
] * motionData [ 63 ] - xx [ 1 ] * motionData [ 65 ] ; xx [ 117 ] = xx [ 3 ]
* xx [ 116 ] ; xx [ 119 ] = xx [ 1 ] * motionData [ 64 ] + xx [ 110 ] *
motionData [ 66 ] ; xx [ 120 ] = xx [ 119 ] * xx [ 3 ] ; xx [ 121 ] = xx [
110 ] * motionData [ 64 ] - xx [ 1 ] * motionData [ 66 ] ; xx [ 122 ] = xx [
16 ] * xx [ 110 ] ; xx [ 123 ] = xx [ 16 ] - xx [ 6 ] * xx [ 122 ] * xx [ 110
] ; xx [ 110 ] = xx [ 123 ] * motionData [ 65 ] ; xx [ 128 ] = motionData [
64 ] ; xx [ 129 ] = motionData [ 65 ] ; xx [ 130 ] = motionData [ 66 ] ; xx [
124 ] = xx [ 6 ] * xx [ 1 ] * xx [ 122 ] ; xx [ 1 ] = xx [ 124 ] * motionData
[ 66 ] + xx [ 123 ] * motionData [ 64 ] ; xx [ 122 ] = xx [ 124 ] *
motionData [ 65 ] ; xx [ 131 ] = xx [ 110 ] ; xx [ 132 ] = - xx [ 1 ] ; xx [
133 ] = xx [ 122 ] ; pm_math_Vector3_cross_ra ( xx + 128 , xx + 131 , xx +
134 ) ; J [ 0 ] = xx [ 14 ] ; J [ 1 ] = xx [ 24 ] ; J [ 2 ] = xx [ 31 ] ; J [
3 ] = xx [ 35 ] ; J [ 4 ] = xx [ 46 ] ; J [ 5 ] = xx [ 53 ] ; J [ 6 ] = xx [
57 ] ; J [ 7 ] = xx [ 39 ] ; J [ 8 ] = xx [ 15 ] ; J [ 9 ] = xx [ 25 ] ; J [
10 ] = xx [ 32 ] ; J [ 11 ] = xx [ 36 ] ; J [ 12 ] = xx [ 47 ] ; J [ 13 ] =
xx [ 54 ] ; J [ 14 ] = xx [ 58 ] ; J [ 15 ] = xx [ 40 ] ; J [ 16 ] = xx [ 7 ]
; J [ 17 ] = xx [ 6 ] * ( xx [ 59 ] + xx [ 9 ] * xx [ 48 ] ) + xx [ 23 ] + (
xx [ 0 ] * motionData [ 0 ] + xx [ 62 ] ) * xx [ 6 ] ; J [ 18 ] = xx [ 6 ] *
( xx [ 49 ] + xx [ 20 ] * xx [ 65 ] ) + xx [ 42 ] + ( xx [ 18 ] * motionData
[ 77 ] + xx [ 72 ] ) * xx [ 6 ] ; J [ 19 ] = xx [ 45 ] + ( xx [ 38 ] *
motionData [ 84 ] + xx [ 75 ] ) * xx [ 6 ] - ( xx [ 67 ] * xx [ 68 ] + xx [
69 ] * xx [ 41 ] ) * xx [ 6 ] ; J [ 20 ] = - ( xx [ 6 ] * ( xx [ 89 ] - xx [
44 ] * xx [ 19 ] ) + xx [ 6 ] * ( xx [ 81 ] * xx [ 71 ] - xx [ 78 ] * xx [ 82
] ) - motionData [ 111 ] - xx [ 6 ] * xx [ 78 ] * xx [ 83 ] ) ; J [ 21 ] = -
( xx [ 6 ] * ( xx [ 97 ] - xx [ 27 ] * xx [ 100 ] ) + xx [ 101 ] + ( xx [ 87
] * motionData [ 28 ] + xx [ 107 ] ) * xx [ 6 ] ) ; J [ 22 ] = - ( xx [ 6 ] *
( xx [ 113 ] - xx [ 29 ] * xx [ 118 ] ) + xx [ 111 ] + ( xx [ 95 ] *
motionData [ 56 ] + xx [ 125 ] ) * xx [ 6 ] ) ; J [ 23 ] = ( xx [ 106 ] * xx
[ 117 ] + xx [ 120 ] * xx [ 121 ] ) * xx [ 6 ] - ( ( xx [ 110 ] * motionData
[ 63 ] + xx [ 134 ] ) * xx [ 6 ] - xx [ 124 ] ) ; J [ 24 ] = xx [ 5 ] ; J [
25 ] = ( xx [ 60 ] - xx [ 48 ] * xx [ 10 ] ) * xx [ 6 ] - xx [ 8 ] + ( xx [ 2
] * motionData [ 0 ] + xx [ 63 ] ) * xx [ 6 ] ; J [ 26 ] = ( xx [ 50 ] - xx [
65 ] * xx [ 33 ] ) * xx [ 6 ] - xx [ 17 ] + ( xx [ 30 ] * motionData [ 77 ] +
xx [ 73 ] ) * xx [ 6 ] ; J [ 27 ] = xx [ 6 ] * ( xx [ 69 ] * xx [ 67 ] - xx [
41 ] * xx [ 68 ] ) + ( xx [ 55 ] * motionData [ 84 ] + xx [ 76 ] ) * xx [ 6 ]
; J [ 28 ] = - ( xx [ 85 ] + xx [ 6 ] * ( xx [ 90 ] - xx [ 44 ] * xx [ 79 ] )
) ; J [ 29 ] = - ( xx [ 22 ] + ( xx [ 100 ] * xx [ 92 ] + xx [ 98 ] ) * xx [
6 ] + ( xx [ 88 ] * motionData [ 28 ] + xx [ 108 ] ) * xx [ 6 ] ) ; J [ 30 ]
= - ( xx [ 28 ] + ( xx [ 118 ] * xx [ 84 ] + xx [ 114 ] ) * xx [ 6 ] + ( xx [
96 ] * motionData [ 56 ] + xx [ 126 ] ) * xx [ 6 ] ) ; J [ 31 ] = - ( xx [ 6
] * ( xx [ 117 ] * xx [ 121 ] - xx [ 106 ] * xx [ 120 ] ) + xx [ 6 ] * ( xx [
135 ] - xx [ 1 ] * motionData [ 63 ] ) ) ; J [ 32 ] = xx [ 5 ] ; J [ 33 ] = (
xx [ 61 ] - xx [ 48 ] * xx [ 11 ] ) * xx [ 6 ] + xx [ 4 ] + xx [ 26 ] + xx [
6 ] * ( xx [ 64 ] - xx [ 12 ] * motionData [ 0 ] ) ; J [ 34 ] = ( xx [ 51 ] -
xx [ 65 ] * xx [ 34 ] ) * xx [ 6 ] + xx [ 13 ] + xx [ 43 ] + xx [ 6 ] * ( xx
[ 74 ] - xx [ 37 ] * motionData [ 77 ] ) ; J [ 35 ] = xx [ 52 ] + xx [ 6 ] *
( xx [ 77 ] - xx [ 56 ] * motionData [ 84 ] ) - ( xx [ 41 ] * xx [ 70 ] + xx
[ 66 ] * xx [ 67 ] ) * xx [ 6 ] + xx [ 3 ] ; J [ 36 ] = - ( xx [ 6 ] * ( xx [
91 ] - xx [ 44 ] * xx [ 80 ] ) + motionData [ 109 ] - ( ( xx [ 78 ] * xx [ 81
] + xx [ 82 ] * xx [ 71 ] ) * xx [ 6 ] + xx [ 6 ] * xx [ 83 ] * xx [ 71 ] ) -
xx [ 21 ] + xx [ 16 ] ) ; J [ 37 ] = - ( ( xx [ 100 ] * xx [ 93 ] + xx [ 99 ]
) * xx [ 6 ] - xx [ 86 ] + xx [ 102 ] + xx [ 6 ] * ( xx [ 109 ] - xx [ 103 ]
* motionData [ 28 ] ) ) ; J [ 38 ] = - ( ( xx [ 118 ] * xx [ 104 ] + xx [ 115
] ) * xx [ 6 ] - xx [ 94 ] + xx [ 112 ] + xx [ 6 ] * ( xx [ 127 ] - xx [ 105
] * motionData [ 56 ] ) ) ; J [ 39 ] = - ( xx [ 123 ] + ( xx [ 122 ] *
motionData [ 63 ] + xx [ 136 ] ) * xx [ 6 ] - ( xx [ 119 ] * xx [ 120 ] + xx
[ 117 ] * xx [ 116 ] ) * xx [ 6 ] + xx [ 3 ] ) ; return 5 ; } size_t
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computeAssemblyJacobian ( const
void * mech , const RuntimeDerivedValuesBundle * rtdv , size_t constraintIdx
, boolean_T forVelocitySatisfaction , const double * state , const int *
modeVector , const double * motionData , double * J ) { ( void ) mech ; (
void ) rtdv ; ( void ) state ; ( void ) modeVector ; ( void )
forVelocitySatisfaction ; ( void ) motionData ; ( void ) J ; switch (
constraintIdx ) { case 0 : return computeAssemblyJacobian_0 ( rtdv , state ,
modeVector , motionData , J ) ; } return 0 ; } size_t
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computeFullAssemblyJacobian (
const void * mech , const RuntimeDerivedValuesBundle * rtdv , const double *
state , const int * modeVector , const double * motionData , double * J ) {
const double * rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv
-> mInts . mValues ; double xx [ 137 ] ; ( void ) mech ; ( void ) rtdvd ; (
void ) rtdvi ; ( void ) modeVector ; xx [ 0 ] = 0.7071067811865476 ; xx [ 1 ]
= xx [ 0 ] * motionData [ 70 ] ; xx [ 2 ] = xx [ 0 ] * motionData [ 71 ] ; xx
[ 3 ] = xx [ 0 ] * motionData [ 72 ] ; xx [ 4 ] = xx [ 0 ] * motionData [ 73
] ; xx [ 5 ] = xx [ 1 ] - xx [ 2 ] ; xx [ 6 ] = xx [ 1 ] + xx [ 2 ] ; xx [ 7
] = xx [ 3 ] + xx [ 4 ] ; xx [ 8 ] = xx [ 4 ] - xx [ 3 ] ; xx [ 1 ] = xx [ 0
] * motionData [ 91 ] ; xx [ 2 ] = xx [ 0 ] * motionData [ 92 ] ; xx [ 3 ] =
xx [ 0 ] * motionData [ 93 ] ; xx [ 4 ] = xx [ 0 ] * motionData [ 94 ] ; xx [
9 ] = xx [ 1 ] - xx [ 2 ] ; xx [ 10 ] = xx [ 1 ] + xx [ 2 ] ; xx [ 11 ] = xx
[ 3 ] + xx [ 4 ] ; xx [ 12 ] = xx [ 4 ] - xx [ 3 ] ;
pm_math_Quaternion_inverseCompose_ra ( xx + 5 , xx + 9 , xx + 1 ) ; xx [ 5 ]
= 0.0 ; xx [ 6 ] = xx [ 5 ] ; xx [ 7 ] = xx [ 5 ] ; xx [ 8 ] = xx [ 5 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 6 , xx + 13 ) ; xx [ 6 ] =
2.0 ; xx [ 7 ] = 1.0 ; xx [ 8 ] = xx [ 7 ] - ( motionData [ 115 ] *
motionData [ 115 ] + motionData [ 113 ] * motionData [ 113 ] ) * xx [ 6 ] ;
xx [ 17 ] = xx [ 0 ] * xx [ 0 ] * xx [ 8 ] ; xx [ 18 ] = xx [ 6 ] * (
motionData [ 114 ] * motionData [ 115 ] - motionData [ 112 ] * motionData [
113 ] ) ; xx [ 19 ] = xx [ 0 ] * xx [ 0 ] * xx [ 18 ] ; xx [ 20 ] = (
motionData [ 112 ] * motionData [ 115 ] + motionData [ 113 ] * motionData [
114 ] ) * xx [ 6 ] ; xx [ 21 ] = xx [ 8 ] - ( xx [ 17 ] - xx [ 19 ] ) * xx [
6 ] ; xx [ 22 ] = xx [ 18 ] - xx [ 6 ] * ( xx [ 17 ] + xx [ 19 ] ) ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 20 , xx + 23 ) ; xx [ 17 ] =
xx [ 7 ] - ( motionData [ 24 ] * motionData [ 24 ] + motionData [ 22 ] *
motionData [ 22 ] ) * xx [ 6 ] ; xx [ 19 ] = xx [ 0 ] * xx [ 0 ] * xx [ 17 ]
; xx [ 20 ] = xx [ 6 ] * ( motionData [ 23 ] * motionData [ 24 ] - motionData
[ 21 ] * motionData [ 22 ] ) ; xx [ 21 ] = xx [ 0 ] * xx [ 0 ] * xx [ 20 ] ;
xx [ 27 ] = ( motionData [ 21 ] * motionData [ 24 ] + motionData [ 22 ] *
motionData [ 23 ] ) * xx [ 6 ] ; xx [ 28 ] = xx [ 17 ] - ( xx [ 19 ] - xx [
21 ] ) * xx [ 6 ] ; xx [ 29 ] = xx [ 20 ] - xx [ 6 ] * ( xx [ 19 ] + xx [ 21
] ) ; pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 27 , xx + 30 ) ; xx [
27 ] = xx [ 5 ] ; xx [ 28 ] = - 2.220446049250313e-16 ; xx [ 29 ] = - 1.0 ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 27 , xx + 34 ) ; xx [ 38 ] =
motionData [ 70 ] ; xx [ 39 ] = motionData [ 71 ] ; xx [ 40 ] = motionData [
72 ] ; xx [ 41 ] = motionData [ 73 ] ; xx [ 19 ] = ( motionData [ 108 ] *
motionData [ 108 ] + motionData [ 106 ] * motionData [ 106 ] ) * xx [ 6 ] -
xx [ 7 ] ; xx [ 21 ] = motionData [ 105 ] * motionData [ 106 ] - motionData [
107 ] * motionData [ 108 ] ; xx [ 27 ] = - ( ( motionData [ 105 ] *
motionData [ 108 ] + motionData [ 106 ] * motionData [ 107 ] ) * xx [ 6 ] ) ;
xx [ 28 ] = xx [ 19 ] ; xx [ 29 ] = xx [ 6 ] * xx [ 21 ] ;
pm_math_Quaternion_xform_ra ( xx + 38 , xx + 27 , xx + 42 ) ;
pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 42 , xx + 27 ) ; xx [ 42 ]
= - xx [ 27 ] ; xx [ 43 ] = - xx [ 28 ] ; xx [ 44 ] = - xx [ 29 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 42 , xx + 45 ) ; xx [ 22 ] =
( motionData [ 101 ] * motionData [ 101 ] + motionData [ 99 ] * motionData [
99 ] ) * xx [ 6 ] - xx [ 7 ] ; xx [ 27 ] = motionData [ 98 ] * motionData [
99 ] - motionData [ 100 ] * motionData [ 101 ] ; xx [ 42 ] = - ( ( motionData
[ 98 ] * motionData [ 101 ] + motionData [ 99 ] * motionData [ 100 ] ) * xx [
6 ] ) ; xx [ 43 ] = xx [ 22 ] ; xx [ 44 ] = xx [ 6 ] * xx [ 27 ] ;
pm_math_Quaternion_xform_ra ( xx + 38 , xx + 42 , xx + 49 ) ;
pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 49 , xx + 42 ) ; xx [ 49 ]
= - xx [ 42 ] ; xx [ 50 ] = - xx [ 43 ] ; xx [ 51 ] = - xx [ 44 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 49 , xx + 52 ) ; xx [ 28 ] =
( motionData [ 52 ] * motionData [ 52 ] + motionData [ 50 ] * motionData [ 50
] ) * xx [ 6 ] - xx [ 7 ] ; xx [ 29 ] = motionData [ 49 ] * motionData [ 50 ]
- motionData [ 51 ] * motionData [ 52 ] ; xx [ 42 ] = - ( ( motionData [ 49 ]
* motionData [ 52 ] + motionData [ 50 ] * motionData [ 51 ] ) * xx [ 6 ] ) ;
xx [ 43 ] = xx [ 28 ] ; xx [ 44 ] = xx [ 6 ] * xx [ 29 ] ;
pm_math_Quaternion_xform_ra ( xx + 38 , xx + 42 , xx + 49 ) ;
pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 49 , xx + 38 ) ; xx [ 41 ]
= - xx [ 38 ] ; xx [ 42 ] = - xx [ 39 ] ; xx [ 43 ] = - xx [ 40 ] ;
pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 41 , xx + 56 ) ; xx [ 38 ] =
xx [ 6 ] * ( motionData [ 70 ] * motionData [ 73 ] - motionData [ 71 ] *
motionData [ 72 ] ) ; xx [ 39 ] = ( motionData [ 73 ] * motionData [ 73 ] +
motionData [ 71 ] * motionData [ 71 ] ) * xx [ 6 ] - xx [ 7 ] ; xx [ 40 ] = -
( ( motionData [ 70 ] * motionData [ 71 ] + motionData [ 72 ] * motionData [
73 ] ) * xx [ 6 ] ) ; pm_math_Quaternion_inverseXform_ra ( xx + 9 , xx + 38 ,
xx + 41 ) ; xx [ 9 ] = - xx [ 41 ] ; xx [ 10 ] = - xx [ 42 ] ; xx [ 11 ] = -
xx [ 43 ] ; pm_math_Quaternion_compDeriv_ra ( xx + 1 , xx + 9 , xx + 38 ) ;
xx [ 1 ] = 0.5 ; xx [ 2 ] = xx [ 1 ] * state [ 2 ] ; xx [ 3 ] = xx [ 0 ] *
cos ( xx [ 2 ] ) ; xx [ 4 ] = xx [ 0 ] * sin ( xx [ 2 ] ) ; xx [ 0 ] = xx [ 3
] + xx [ 4 ] ; xx [ 2 ] = xx [ 3 ] - xx [ 4 ] ; xx [ 9 ] = - ( xx [ 0 ] *
motionData [ 0 ] + xx [ 2 ] * motionData [ 2 ] ) ; xx [ 10 ] = - ( xx [ 0 ] *
motionData [ 1 ] + xx [ 2 ] * motionData [ 3 ] ) ; xx [ 11 ] = xx [ 2 ] *
motionData [ 0 ] - xx [ 0 ] * motionData [ 2 ] ; xx [ 12 ] = xx [ 2 ] *
motionData [ 1 ] - xx [ 0 ] * motionData [ 3 ] ; xx [ 41 ] = motionData [ 112
] ; xx [ 42 ] = motionData [ 113 ] ; xx [ 43 ] = motionData [ 114 ] ; xx [ 44
] = motionData [ 115 ] ; pm_math_Quaternion_compose_ra ( xx + 9 , xx + 41 ,
xx + 48 ) ; xx [ 3 ] = 1.0e-3 ; xx [ 4 ] = xx [ 3 ] * xx [ 8 ] ; xx [ 8 ] =
xx [ 3 ] * xx [ 18 ] ; xx [ 9 ] = xx [ 50 ] * xx [ 4 ] + xx [ 51 ] * xx [ 8 ]
; xx [ 10 ] = xx [ 49 ] * xx [ 4 ] ; xx [ 11 ] = xx [ 49 ] * xx [ 8 ] ; xx [
41 ] = xx [ 9 ] ; xx [ 42 ] = - xx [ 10 ] ; xx [ 43 ] = - xx [ 11 ] ;
pm_math_Vector3_cross_ra ( xx + 49 , xx + 41 , xx + 59 ) ; xx [ 12 ] = xx [ 2
] * motionData [ 118 ] ; xx [ 13 ] = xx [ 2 ] * motionData [ 116 ] ; xx [ 16
] = 1.0e-3 ; xx [ 18 ] = xx [ 2 ] * xx [ 16 ] ; xx [ 23 ] = motionData [ 118
] - ( xx [ 2 ] * xx [ 12 ] - xx [ 13 ] * xx [ 0 ] ) * xx [ 6 ] - xx [ 6 ] *
xx [ 18 ] * xx [ 0 ] ; xx [ 26 ] = xx [ 6 ] * ( xx [ 2 ] * xx [ 13 ] + xx [
12 ] * xx [ 0 ] ) - ( motionData [ 116 ] + xx [ 6 ] * xx [ 2 ] * xx [ 18 ] )
+ xx [ 16 ] ; xx [ 0 ] = xx [ 26 ] * motionData [ 2 ] ; xx [ 41 ] =
motionData [ 1 ] ; xx [ 42 ] = motionData [ 2 ] ; xx [ 43 ] = motionData [ 3
] ; xx [ 2 ] = xx [ 23 ] * motionData [ 3 ] - xx [ 26 ] * motionData [ 1 ] ;
xx [ 12 ] = xx [ 23 ] * motionData [ 2 ] ; xx [ 49 ] = xx [ 0 ] ; xx [ 50 ] =
xx [ 2 ] ; xx [ 51 ] = - xx [ 12 ] ; pm_math_Vector3_cross_ra ( xx + 41 , xx
+ 49 , xx + 62 ) ; xx [ 13 ] = xx [ 1 ] * state [ 4 ] ; xx [ 18 ] = sin ( xx
[ 13 ] ) ; xx [ 30 ] = cos ( xx [ 13 ] ) ; xx [ 41 ] = xx [ 18 ] * motionData
[ 79 ] - xx [ 30 ] * motionData [ 77 ] ; xx [ 42 ] = xx [ 18 ] * motionData [
80 ] - xx [ 30 ] * motionData [ 78 ] ; xx [ 43 ] = - ( xx [ 18 ] * motionData
[ 77 ] + xx [ 30 ] * motionData [ 79 ] ) ; xx [ 44 ] = - ( xx [ 30 ] *
motionData [ 80 ] + xx [ 18 ] * motionData [ 78 ] ) ; xx [ 49 ] = motionData
[ 21 ] ; xx [ 50 ] = motionData [ 22 ] ; xx [ 51 ] = motionData [ 23 ] ; xx [
52 ] = motionData [ 24 ] ; pm_math_Quaternion_compose_ra ( xx + 41 , xx + 49
, xx + 65 ) ; xx [ 13 ] = xx [ 3 ] * xx [ 17 ] ; xx [ 17 ] = xx [ 3 ] * xx [
20 ] ; xx [ 20 ] = xx [ 67 ] * xx [ 13 ] + xx [ 68 ] * xx [ 17 ] ; xx [ 33 ]
= xx [ 66 ] * xx [ 13 ] ; xx [ 34 ] = xx [ 66 ] * xx [ 17 ] ; xx [ 41 ] = xx
[ 20 ] ; xx [ 42 ] = - xx [ 33 ] ; xx [ 43 ] = - xx [ 34 ] ;
pm_math_Vector3_cross_ra ( xx + 66 , xx + 41 , xx + 49 ) ; xx [ 37 ] = xx [
18 ] * motionData [ 25 ] ; xx [ 38 ] = xx [ 18 ] * motionData [ 27 ] ; xx [
41 ] = xx [ 16 ] * xx [ 18 ] ; xx [ 42 ] = motionData [ 27 ] - ( xx [ 30 ] *
xx [ 37 ] + xx [ 38 ] * xx [ 18 ] ) * xx [ 6 ] + xx [ 6 ] * xx [ 30 ] * xx [
41 ] ; xx [ 43 ] = xx [ 6 ] * ( xx [ 37 ] * xx [ 18 ] - xx [ 30 ] * xx [ 38 ]
) - ( motionData [ 25 ] + xx [ 6 ] * xx [ 41 ] * xx [ 18 ] ) + xx [ 16 ] ; xx
[ 18 ] = xx [ 43 ] * motionData [ 79 ] ; xx [ 66 ] = motionData [ 78 ] ; xx [
67 ] = motionData [ 79 ] ; xx [ 68 ] = motionData [ 80 ] ; xx [ 30 ] = xx [
42 ] * motionData [ 80 ] - xx [ 43 ] * motionData [ 78 ] ; xx [ 37 ] = xx [
42 ] * motionData [ 79 ] ; xx [ 69 ] = xx [ 18 ] ; xx [ 70 ] = xx [ 30 ] ; xx
[ 71 ] = - xx [ 37 ] ; pm_math_Vector3_cross_ra ( xx + 66 , xx + 69 , xx + 72
) ; xx [ 38 ] = xx [ 1 ] * state [ 6 ] ; xx [ 41 ] = cos ( xx [ 38 ] ) ; xx [
44 ] = sin ( xx [ 38 ] ) ; xx [ 38 ] = xx [ 16 ] * xx [ 44 ] ; xx [ 45 ] = xx
[ 6 ] * xx [ 41 ] * xx [ 38 ] ; xx [ 52 ] = xx [ 16 ] - xx [ 6 ] * xx [ 38 ]
* xx [ 44 ] ; xx [ 38 ] = xx [ 52 ] * motionData [ 86 ] ; xx [ 66 ] =
motionData [ 85 ] ; xx [ 67 ] = motionData [ 86 ] ; xx [ 68 ] = motionData [
87 ] ; xx [ 55 ] = xx [ 45 ] * motionData [ 87 ] - xx [ 52 ] * motionData [
85 ] ; xx [ 56 ] = xx [ 45 ] * motionData [ 86 ] ; xx [ 69 ] = xx [ 38 ] ; xx
[ 70 ] = xx [ 55 ] ; xx [ 71 ] = - xx [ 56 ] ; pm_math_Vector3_cross_ra ( xx
+ 66 , xx + 69 , xx + 75 ) ; xx [ 66 ] = xx [ 44 ] * motionData [ 84 ] + xx [
41 ] * motionData [ 86 ] ; xx [ 67 ] = xx [ 66 ] * xx [ 3 ] ; xx [ 68 ] = xx
[ 44 ] * motionData [ 86 ] - xx [ 41 ] * motionData [ 84 ] ; xx [ 69 ] = xx [
41 ] * motionData [ 87 ] + xx [ 44 ] * motionData [ 85 ] ; xx [ 70 ] = xx [
44 ] * motionData [ 87 ] - xx [ 41 ] * motionData [ 85 ] ; xx [ 41 ] = xx [ 3
] * xx [ 70 ] ; xx [ 44 ] = xx [ 1 ] * state [ 8 ] ; xx [ 71 ] = sin ( xx [
44 ] ) ; xx [ 78 ] = cos ( xx [ 44 ] ) ; xx [ 44 ] = xx [ 71 ] * motionData [
108 ] - xx [ 78 ] * motionData [ 106 ] ; xx [ 79 ] = xx [ 71 ] * motionData [
105 ] - xx [ 78 ] * motionData [ 107 ] ; xx [ 80 ] = xx [ 78 ] * motionData [
108 ] + xx [ 71 ] * motionData [ 106 ] ; xx [ 81 ] = xx [ 44 ] ; xx [ 82 ] =
xx [ 79 ] ; xx [ 83 ] = - xx [ 80 ] ; xx [ 84 ] = 2.0e-3 ; xx [ 85 ] = xx [
84 ] * xx [ 21 ] ; xx [ 21 ] = xx [ 19 ] * xx [ 3 ] ; xx [ 19 ] = xx [ 80 ] *
xx [ 85 ] - xx [ 21 ] * xx [ 79 ] ; xx [ 79 ] = xx [ 21 ] * xx [ 44 ] ; xx [
80 ] = xx [ 85 ] * xx [ 44 ] ; xx [ 86 ] = xx [ 19 ] ; xx [ 87 ] = xx [ 79 ]
; xx [ 88 ] = xx [ 80 ] ; pm_math_Vector3_cross_ra ( xx + 81 , xx + 86 , xx +
89 ) ; xx [ 44 ] = xx [ 78 ] * motionData [ 105 ] + xx [ 71 ] * motionData [
107 ] ; xx [ 81 ] = xx [ 71 ] * motionData [ 111 ] ; xx [ 82 ] = xx [ 71 ] *
motionData [ 109 ] ; xx [ 83 ] = xx [ 16 ] * xx [ 71 ] ; xx [ 86 ] = xx [ 1 ]
* state [ 10 ] ; xx [ 87 ] = cos ( xx [ 86 ] ) ; xx [ 88 ] = sin ( xx [ 86 ]
) ; xx [ 92 ] = - ( xx [ 87 ] * motionData [ 28 ] + xx [ 88 ] * motionData [
30 ] ) ; xx [ 93 ] = - ( xx [ 87 ] * motionData [ 29 ] + xx [ 88 ] *
motionData [ 31 ] ) ; xx [ 94 ] = xx [ 88 ] * motionData [ 28 ] - xx [ 87 ] *
motionData [ 30 ] ; xx [ 95 ] = xx [ 88 ] * motionData [ 29 ] - xx [ 87 ] *
motionData [ 31 ] ; xx [ 96 ] = motionData [ 98 ] ; xx [ 97 ] = motionData [
99 ] ; xx [ 98 ] = motionData [ 100 ] ; xx [ 99 ] = motionData [ 101 ] ;
pm_math_Quaternion_compose_ra ( xx + 92 , xx + 96 , xx + 100 ) ; xx [ 86 ] =
xx [ 22 ] * xx [ 3 ] ; xx [ 22 ] = xx [ 84 ] * xx [ 27 ] ; xx [ 27 ] = xx [
102 ] * xx [ 86 ] + xx [ 103 ] * xx [ 22 ] ; xx [ 92 ] = xx [ 101 ] * xx [ 86
] ; xx [ 93 ] = xx [ 101 ] * xx [ 22 ] ; xx [ 94 ] = - xx [ 27 ] ; xx [ 95 ]
= xx [ 92 ] ; xx [ 96 ] = xx [ 93 ] ; pm_math_Vector3_cross_ra ( xx + 101 ,
xx + 94 , xx + 97 ) ; xx [ 94 ] = xx [ 88 ] * motionData [ 104 ] ; xx [ 95 ]
= xx [ 88 ] * motionData [ 102 ] ; xx [ 96 ] = xx [ 16 ] * xx [ 88 ] ; xx [
101 ] = xx [ 6 ] * ( xx [ 94 ] * xx [ 88 ] - xx [ 87 ] * xx [ 95 ] ) -
motionData [ 104 ] - xx [ 6 ] * xx [ 87 ] * xx [ 96 ] ; xx [ 102 ] =
motionData [ 102 ] - ( ( xx [ 87 ] * xx [ 94 ] + xx [ 95 ] * xx [ 88 ] ) * xx
[ 6 ] + xx [ 6 ] * xx [ 96 ] * xx [ 88 ] ) + xx [ 16 ] ; xx [ 87 ] = xx [ 102
] * motionData [ 30 ] ; xx [ 94 ] = motionData [ 29 ] ; xx [ 95 ] =
motionData [ 30 ] ; xx [ 96 ] = motionData [ 31 ] ; xx [ 88 ] = xx [ 101 ] *
motionData [ 31 ] - xx [ 102 ] * motionData [ 29 ] ; xx [ 103 ] = xx [ 101 ]
* motionData [ 30 ] ; xx [ 104 ] = xx [ 87 ] ; xx [ 105 ] = xx [ 88 ] ; xx [
106 ] = - xx [ 103 ] ; pm_math_Vector3_cross_ra ( xx + 94 , xx + 104 , xx +
107 ) ; xx [ 94 ] = xx [ 1 ] * state [ 12 ] ; xx [ 95 ] = cos ( xx [ 94 ] ) ;
xx [ 96 ] = sin ( xx [ 94 ] ) ; xx [ 110 ] = - ( xx [ 95 ] * motionData [ 56
] + xx [ 96 ] * motionData [ 58 ] ) ; xx [ 111 ] = - ( xx [ 95 ] * motionData
[ 57 ] + xx [ 96 ] * motionData [ 59 ] ) ; xx [ 112 ] = xx [ 96 ] *
motionData [ 56 ] - xx [ 95 ] * motionData [ 58 ] ; xx [ 113 ] = xx [ 96 ] *
motionData [ 57 ] - xx [ 95 ] * motionData [ 59 ] ; xx [ 114 ] = motionData [
49 ] ; xx [ 115 ] = motionData [ 50 ] ; xx [ 116 ] = motionData [ 51 ] ; xx [
117 ] = motionData [ 52 ] ; pm_math_Quaternion_compose_ra ( xx + 110 , xx +
114 , xx + 118 ) ; xx [ 94 ] = xx [ 28 ] * xx [ 3 ] ; xx [ 28 ] = xx [ 84 ] *
xx [ 29 ] ; xx [ 29 ] = xx [ 120 ] * xx [ 94 ] + xx [ 121 ] * xx [ 28 ] ; xx
[ 84 ] = xx [ 119 ] * xx [ 94 ] ; xx [ 104 ] = xx [ 119 ] * xx [ 28 ] ; xx [
110 ] = - xx [ 29 ] ; xx [ 111 ] = xx [ 84 ] ; xx [ 112 ] = xx [ 104 ] ;
pm_math_Vector3_cross_ra ( xx + 119 , xx + 110 , xx + 113 ) ; xx [ 105 ] = xx
[ 96 ] * motionData [ 55 ] ; xx [ 106 ] = xx [ 96 ] * motionData [ 53 ] ; xx
[ 110 ] = xx [ 16 ] * xx [ 96 ] ; xx [ 111 ] = xx [ 6 ] * ( xx [ 105 ] * xx [
96 ] - xx [ 95 ] * xx [ 106 ] ) - motionData [ 55 ] - xx [ 6 ] * xx [ 95 ] *
xx [ 110 ] ; xx [ 112 ] = motionData [ 53 ] - ( ( xx [ 95 ] * xx [ 105 ] + xx
[ 106 ] * xx [ 96 ] ) * xx [ 6 ] + xx [ 6 ] * xx [ 110 ] * xx [ 96 ] ) + xx [
16 ] ; xx [ 95 ] = xx [ 112 ] * motionData [ 58 ] ; xx [ 119 ] = motionData [
57 ] ; xx [ 120 ] = motionData [ 58 ] ; xx [ 121 ] = motionData [ 59 ] ; xx [
96 ] = xx [ 111 ] * motionData [ 59 ] - xx [ 112 ] * motionData [ 57 ] ; xx [
105 ] = xx [ 111 ] * motionData [ 58 ] ; xx [ 122 ] = xx [ 95 ] ; xx [ 123 ]
= xx [ 96 ] ; xx [ 124 ] = - xx [ 105 ] ; pm_math_Vector3_cross_ra ( xx + 119
, xx + 122 , xx + 125 ) ; xx [ 106 ] = xx [ 1 ] * state [ 14 ] ; xx [ 1 ] =
cos ( xx [ 106 ] ) ; xx [ 110 ] = sin ( xx [ 106 ] ) ; xx [ 106 ] = xx [ 1 ]
* motionData [ 63 ] + xx [ 110 ] * motionData [ 65 ] ; xx [ 116 ] = xx [ 110
] * motionData [ 63 ] - xx [ 1 ] * motionData [ 65 ] ; xx [ 117 ] = xx [ 3 ]
* xx [ 116 ] ; xx [ 119 ] = xx [ 1 ] * motionData [ 64 ] + xx [ 110 ] *
motionData [ 66 ] ; xx [ 120 ] = xx [ 119 ] * xx [ 3 ] ; xx [ 121 ] = xx [
110 ] * motionData [ 64 ] - xx [ 1 ] * motionData [ 66 ] ; xx [ 122 ] = xx [
16 ] * xx [ 110 ] ; xx [ 123 ] = xx [ 16 ] - xx [ 6 ] * xx [ 122 ] * xx [ 110
] ; xx [ 110 ] = xx [ 123 ] * motionData [ 65 ] ; xx [ 128 ] = motionData [
64 ] ; xx [ 129 ] = motionData [ 65 ] ; xx [ 130 ] = motionData [ 66 ] ; xx [
124 ] = xx [ 6 ] * xx [ 1 ] * xx [ 122 ] ; xx [ 1 ] = xx [ 124 ] * motionData
[ 66 ] + xx [ 123 ] * motionData [ 64 ] ; xx [ 122 ] = xx [ 124 ] *
motionData [ 65 ] ; xx [ 131 ] = xx [ 110 ] ; xx [ 132 ] = - xx [ 1 ] ; xx [
133 ] = xx [ 122 ] ; pm_math_Vector3_cross_ra ( xx + 128 , xx + 131 , xx +
134 ) ; J [ 0 ] = xx [ 14 ] ; J [ 1 ] = xx [ 24 ] ; J [ 2 ] = xx [ 31 ] ; J [
3 ] = xx [ 35 ] ; J [ 4 ] = xx [ 46 ] ; J [ 5 ] = xx [ 53 ] ; J [ 6 ] = xx [
57 ] ; J [ 7 ] = xx [ 39 ] ; J [ 8 ] = xx [ 15 ] ; J [ 9 ] = xx [ 25 ] ; J [
10 ] = xx [ 32 ] ; J [ 11 ] = xx [ 36 ] ; J [ 12 ] = xx [ 47 ] ; J [ 13 ] =
xx [ 54 ] ; J [ 14 ] = xx [ 58 ] ; J [ 15 ] = xx [ 40 ] ; J [ 16 ] = xx [ 7 ]
; J [ 17 ] = xx [ 6 ] * ( xx [ 59 ] + xx [ 9 ] * xx [ 48 ] ) + xx [ 23 ] + (
xx [ 0 ] * motionData [ 0 ] + xx [ 62 ] ) * xx [ 6 ] ; J [ 18 ] = xx [ 6 ] *
( xx [ 49 ] + xx [ 20 ] * xx [ 65 ] ) + xx [ 42 ] + ( xx [ 18 ] * motionData
[ 77 ] + xx [ 72 ] ) * xx [ 6 ] ; J [ 19 ] = xx [ 45 ] + ( xx [ 38 ] *
motionData [ 84 ] + xx [ 75 ] ) * xx [ 6 ] - ( xx [ 67 ] * xx [ 68 ] + xx [
69 ] * xx [ 41 ] ) * xx [ 6 ] ; J [ 20 ] = - ( xx [ 6 ] * ( xx [ 89 ] - xx [
44 ] * xx [ 19 ] ) + xx [ 6 ] * ( xx [ 81 ] * xx [ 71 ] - xx [ 78 ] * xx [ 82
] ) - motionData [ 111 ] - xx [ 6 ] * xx [ 78 ] * xx [ 83 ] ) ; J [ 21 ] = -
( xx [ 6 ] * ( xx [ 97 ] - xx [ 27 ] * xx [ 100 ] ) + xx [ 101 ] + ( xx [ 87
] * motionData [ 28 ] + xx [ 107 ] ) * xx [ 6 ] ) ; J [ 22 ] = - ( xx [ 6 ] *
( xx [ 113 ] - xx [ 29 ] * xx [ 118 ] ) + xx [ 111 ] + ( xx [ 95 ] *
motionData [ 56 ] + xx [ 125 ] ) * xx [ 6 ] ) ; J [ 23 ] = ( xx [ 106 ] * xx
[ 117 ] + xx [ 120 ] * xx [ 121 ] ) * xx [ 6 ] - ( ( xx [ 110 ] * motionData
[ 63 ] + xx [ 134 ] ) * xx [ 6 ] - xx [ 124 ] ) ; J [ 24 ] = xx [ 5 ] ; J [
25 ] = ( xx [ 60 ] - xx [ 48 ] * xx [ 10 ] ) * xx [ 6 ] - xx [ 8 ] + ( xx [ 2
] * motionData [ 0 ] + xx [ 63 ] ) * xx [ 6 ] ; J [ 26 ] = ( xx [ 50 ] - xx [
65 ] * xx [ 33 ] ) * xx [ 6 ] - xx [ 17 ] + ( xx [ 30 ] * motionData [ 77 ] +
xx [ 73 ] ) * xx [ 6 ] ; J [ 27 ] = xx [ 6 ] * ( xx [ 69 ] * xx [ 67 ] - xx [
41 ] * xx [ 68 ] ) + ( xx [ 55 ] * motionData [ 84 ] + xx [ 76 ] ) * xx [ 6 ]
; J [ 28 ] = - ( xx [ 85 ] + xx [ 6 ] * ( xx [ 90 ] - xx [ 44 ] * xx [ 79 ] )
) ; J [ 29 ] = - ( xx [ 22 ] + ( xx [ 100 ] * xx [ 92 ] + xx [ 98 ] ) * xx [
6 ] + ( xx [ 88 ] * motionData [ 28 ] + xx [ 108 ] ) * xx [ 6 ] ) ; J [ 30 ]
= - ( xx [ 28 ] + ( xx [ 118 ] * xx [ 84 ] + xx [ 114 ] ) * xx [ 6 ] + ( xx [
96 ] * motionData [ 56 ] + xx [ 126 ] ) * xx [ 6 ] ) ; J [ 31 ] = - ( xx [ 6
] * ( xx [ 117 ] * xx [ 121 ] - xx [ 106 ] * xx [ 120 ] ) + xx [ 6 ] * ( xx [
135 ] - xx [ 1 ] * motionData [ 63 ] ) ) ; J [ 32 ] = xx [ 5 ] ; J [ 33 ] = (
xx [ 61 ] - xx [ 48 ] * xx [ 11 ] ) * xx [ 6 ] + xx [ 4 ] + xx [ 26 ] + xx [
6 ] * ( xx [ 64 ] - xx [ 12 ] * motionData [ 0 ] ) ; J [ 34 ] = ( xx [ 51 ] -
xx [ 65 ] * xx [ 34 ] ) * xx [ 6 ] + xx [ 13 ] + xx [ 43 ] + xx [ 6 ] * ( xx
[ 74 ] - xx [ 37 ] * motionData [ 77 ] ) ; J [ 35 ] = xx [ 52 ] + xx [ 6 ] *
( xx [ 77 ] - xx [ 56 ] * motionData [ 84 ] ) - ( xx [ 41 ] * xx [ 70 ] + xx
[ 66 ] * xx [ 67 ] ) * xx [ 6 ] + xx [ 3 ] ; J [ 36 ] = - ( xx [ 6 ] * ( xx [
91 ] - xx [ 44 ] * xx [ 80 ] ) + motionData [ 109 ] - ( ( xx [ 78 ] * xx [ 81
] + xx [ 82 ] * xx [ 71 ] ) * xx [ 6 ] + xx [ 6 ] * xx [ 83 ] * xx [ 71 ] ) -
xx [ 21 ] + xx [ 16 ] ) ; J [ 37 ] = - ( ( xx [ 100 ] * xx [ 93 ] + xx [ 99 ]
) * xx [ 6 ] - xx [ 86 ] + xx [ 102 ] + xx [ 6 ] * ( xx [ 109 ] - xx [ 103 ]
* motionData [ 28 ] ) ) ; J [ 38 ] = - ( ( xx [ 118 ] * xx [ 104 ] + xx [ 115
] ) * xx [ 6 ] - xx [ 94 ] + xx [ 112 ] + xx [ 6 ] * ( xx [ 127 ] - xx [ 105
] * motionData [ 56 ] ) ) ; J [ 39 ] = - ( xx [ 123 ] + ( xx [ 122 ] *
motionData [ 63 ] + xx [ 136 ] ) * xx [ 6 ] - ( xx [ 119 ] * xx [ 120 ] + xx
[ 117 ] * xx [ 116 ] ) * xx [ 6 ] + xx [ 3 ] ) ; return 5 ; } static
boolean_T isInKinematicSingularity_0 ( const RuntimeDerivedValuesBundle *
rtdv , const int * modeVector , const double * motionData ) { const double *
rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts .
mValues ; ( void ) rtdvd ; ( void ) rtdvi ; ( void ) modeVector ; ( void )
motionData ; return 0 ; } boolean_T
MagneticBasket_Simscape_Optimizer_dda62cd9_1_isInKinematicSingularity ( const
void * mech , const RuntimeDerivedValuesBundle * rtdv , size_t constraintIdx
, const int * modeVector , const double * motionData ) { ( void ) mech ; (
void ) rtdv ; ( void ) modeVector ; ( void ) motionData ; switch (
constraintIdx ) { case 0 : return isInKinematicSingularity_0 ( rtdv ,
modeVector , motionData ) ; } return 0 ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_convertStateVector ( const void
* asmMech , const RuntimeDerivedValuesBundle * rtdv , const void * simMech ,
const double * asmState , const int * asmModeVector , const int *
simModeVector , double * simState ) { const double * rtdvd = rtdv -> mDoubles
. mValues ; const int * rtdvi = rtdv -> mInts . mValues ; ( void ) asmMech ;
( void ) rtdvd ; ( void ) rtdvi ; ( void ) simMech ; ( void ) asmModeVector ;
( void ) simModeVector ; simState [ 0 ] = asmState [ 0 ] ; simState [ 1 ] =
asmState [ 1 ] ; simState [ 2 ] = asmState [ 2 ] ; simState [ 3 ] = asmState
[ 3 ] ; simState [ 4 ] = asmState [ 4 ] ; simState [ 5 ] = asmState [ 5 ] ;
simState [ 6 ] = asmState [ 6 ] ; simState [ 7 ] = asmState [ 7 ] ; simState
[ 8 ] = asmState [ 8 ] ; simState [ 9 ] = asmState [ 9 ] ; simState [ 10 ] =
asmState [ 10 ] ; simState [ 11 ] = asmState [ 11 ] ; simState [ 12 ] =
asmState [ 12 ] ; simState [ 13 ] = asmState [ 13 ] ; simState [ 14 ] =
asmState [ 14 ] ; simState [ 15 ] = asmState [ 15 ] ; simState [ 16 ] =
asmState [ 16 ] ; simState [ 17 ] = asmState [ 17 ] ; }
