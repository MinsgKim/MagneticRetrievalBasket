#include <math.h>
#include <string.h>
#include "pm_std.h"
#include "sm_std.h"
#include "ne_std.h"
#include "ne_dae.h"
#include "sm_ssci_run_time_errors.h"
#include "sm_RuntimeDerivedValuesBundle.h"
void MagneticBasket_Simscape_Optimizer_dda62cd9_1_resetSimStateVector ( const
void * mech , double * state ) { double xx [ 1 ] ; ( void ) mech ; xx [ 0 ] =
0.0 ; state [ 0 ] = xx [ 0 ] ; state [ 1 ] = xx [ 0 ] ; state [ 2 ] = xx [ 0
] ; state [ 3 ] = xx [ 0 ] ; state [ 4 ] = xx [ 0 ] ; state [ 5 ] = xx [ 0 ]
; state [ 6 ] = xx [ 0 ] ; state [ 7 ] = xx [ 0 ] ; state [ 8 ] = xx [ 0 ] ;
state [ 9 ] = xx [ 0 ] ; state [ 10 ] = xx [ 0 ] ; state [ 11 ] = xx [ 0 ] ;
state [ 12 ] = xx [ 0 ] ; state [ 13 ] = xx [ 0 ] ; state [ 14 ] = xx [ 0 ] ;
state [ 15 ] = xx [ 0 ] ; state [ 16 ] = xx [ 0 ] ; state [ 17 ] = xx [ 0 ] ;
} static void perturbSimJointPrimitiveState_0_0 ( double mag , double * state
) { state [ 0 ] = state [ 0 ] + mag ; } static void
perturbSimJointPrimitiveState_0_0v ( double mag , double * state ) { state [
0 ] = state [ 0 ] + mag ; state [ 1 ] = state [ 1 ] - 0.875 * mag ; } static
void perturbSimJointPrimitiveState_1_0 ( double mag , double * state ) {
state [ 2 ] = state [ 2 ] + mag ; } static void
perturbSimJointPrimitiveState_1_0v ( double mag , double * state ) { state [
2 ] = state [ 2 ] + mag ; state [ 3 ] = state [ 3 ] - 0.875 * mag ; } static
void perturbSimJointPrimitiveState_2_0 ( double mag , double * state ) {
state [ 4 ] = state [ 4 ] + mag ; } static void
perturbSimJointPrimitiveState_2_0v ( double mag , double * state ) { state [
4 ] = state [ 4 ] + mag ; state [ 5 ] = state [ 5 ] - 0.875 * mag ; } static
void perturbSimJointPrimitiveState_3_0 ( double mag , double * state ) {
state [ 6 ] = state [ 6 ] + mag ; } static void
perturbSimJointPrimitiveState_3_0v ( double mag , double * state ) { state [
6 ] = state [ 6 ] + mag ; state [ 7 ] = state [ 7 ] - 0.875 * mag ; } static
void perturbSimJointPrimitiveState_4_0 ( double mag , double * state ) {
state [ 8 ] = state [ 8 ] + mag ; } static void
perturbSimJointPrimitiveState_4_0v ( double mag , double * state ) { state [
8 ] = state [ 8 ] + mag ; state [ 9 ] = state [ 9 ] - 0.875 * mag ; } static
void perturbSimJointPrimitiveState_5_0 ( double mag , double * state ) {
state [ 10 ] = state [ 10 ] + mag ; } static void
perturbSimJointPrimitiveState_5_0v ( double mag , double * state ) { state [
10 ] = state [ 10 ] + mag ; state [ 11 ] = state [ 11 ] - 0.875 * mag ; }
static void perturbSimJointPrimitiveState_6_0 ( double mag , double * state )
{ state [ 12 ] = state [ 12 ] + mag ; } static void
perturbSimJointPrimitiveState_6_0v ( double mag , double * state ) { state [
12 ] = state [ 12 ] + mag ; state [ 13 ] = state [ 13 ] - 0.875 * mag ; }
static void perturbSimJointPrimitiveState_7_0 ( double mag , double * state )
{ state [ 14 ] = state [ 14 ] + mag ; } static void
perturbSimJointPrimitiveState_7_0v ( double mag , double * state ) { state [
14 ] = state [ 14 ] + mag ; state [ 15 ] = state [ 15 ] - 0.875 * mag ; }
void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_perturbSimJointPrimitiveState (
const void * mech , size_t stageIdx , size_t primIdx , double mag , boolean_T
doPerturbVelocity , double * state ) { ( void ) mech ; ( void ) stageIdx ; (
void ) primIdx ; ( void ) mag ; ( void ) doPerturbVelocity ; ( void ) state ;
switch ( ( stageIdx * 6 + primIdx ) * 2 + ( doPerturbVelocity ? 1 : 0 ) ) {
case 0 : perturbSimJointPrimitiveState_0_0 ( mag , state ) ; break ; case 1 :
perturbSimJointPrimitiveState_0_0v ( mag , state ) ; break ; case 12 :
perturbSimJointPrimitiveState_1_0 ( mag , state ) ; break ; case 13 :
perturbSimJointPrimitiveState_1_0v ( mag , state ) ; break ; case 24 :
perturbSimJointPrimitiveState_2_0 ( mag , state ) ; break ; case 25 :
perturbSimJointPrimitiveState_2_0v ( mag , state ) ; break ; case 36 :
perturbSimJointPrimitiveState_3_0 ( mag , state ) ; break ; case 37 :
perturbSimJointPrimitiveState_3_0v ( mag , state ) ; break ; case 48 :
perturbSimJointPrimitiveState_4_0 ( mag , state ) ; break ; case 49 :
perturbSimJointPrimitiveState_4_0v ( mag , state ) ; break ; case 60 :
perturbSimJointPrimitiveState_5_0 ( mag , state ) ; break ; case 61 :
perturbSimJointPrimitiveState_5_0v ( mag , state ) ; break ; case 72 :
perturbSimJointPrimitiveState_6_0 ( mag , state ) ; break ; case 73 :
perturbSimJointPrimitiveState_6_0v ( mag , state ) ; break ; case 84 :
perturbSimJointPrimitiveState_7_0 ( mag , state ) ; break ; case 85 :
perturbSimJointPrimitiveState_7_0v ( mag , state ) ; break ; } } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_perturbFlexibleBodyState ( const
void * mech , size_t stageIdx , double mag , boolean_T doPerturbVelocity ,
double * state ) { ( void ) mech ; ( void ) stageIdx ; ( void ) mag ; ( void
) doPerturbVelocity ; ( void ) state ; switch ( stageIdx * 2 + (
doPerturbVelocity ? 1 : 0 ) ) { } } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_constructStateVector ( const
void * mech , const double * solverState , const double * u , const double *
uDot , double * discreteState , double * fullState ) { ( void ) mech ; ( void
) u ; ( void ) uDot ; ( void ) discreteState ; fullState [ 0 ] = solverState
[ 0 ] ; fullState [ 1 ] = solverState [ 1 ] ; fullState [ 2 ] = solverState [
2 ] ; fullState [ 3 ] = solverState [ 3 ] ; fullState [ 4 ] = solverState [ 4
] ; fullState [ 5 ] = solverState [ 5 ] ; fullState [ 6 ] = solverState [ 6 ]
; fullState [ 7 ] = solverState [ 7 ] ; fullState [ 8 ] = solverState [ 8 ] ;
fullState [ 9 ] = solverState [ 9 ] ; fullState [ 10 ] = solverState [ 10 ] ;
fullState [ 11 ] = solverState [ 11 ] ; fullState [ 12 ] = solverState [ 12 ]
; fullState [ 13 ] = solverState [ 13 ] ; fullState [ 14 ] = solverState [ 14
] ; fullState [ 15 ] = solverState [ 15 ] ; fullState [ 16 ] = solverState [
16 ] ; fullState [ 17 ] = solverState [ 17 ] ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_extractSolverStateVector ( const
void * mech , const double * fullState , double * solverState ) { ( void )
mech ; solverState [ 0 ] = fullState [ 0 ] ; solverState [ 1 ] = fullState [
1 ] ; solverState [ 2 ] = fullState [ 2 ] ; solverState [ 3 ] = fullState [ 3
] ; solverState [ 4 ] = fullState [ 4 ] ; solverState [ 5 ] = fullState [ 5 ]
; solverState [ 6 ] = fullState [ 6 ] ; solverState [ 7 ] = fullState [ 7 ] ;
solverState [ 8 ] = fullState [ 8 ] ; solverState [ 9 ] = fullState [ 9 ] ;
solverState [ 10 ] = fullState [ 10 ] ; solverState [ 11 ] = fullState [ 11 ]
; solverState [ 12 ] = fullState [ 12 ] ; solverState [ 13 ] = fullState [ 13
] ; solverState [ 14 ] = fullState [ 14 ] ; solverState [ 15 ] = fullState [
15 ] ; solverState [ 16 ] = fullState [ 16 ] ; solverState [ 17 ] = fullState
[ 17 ] ; } boolean_T
MagneticBasket_Simscape_Optimizer_dda62cd9_1_isPositionViolation ( const void
* mech , const RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags
, const double * state , const int * modeVector ) { const double * rtdvd =
rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts . mValues ;
int ii [ 1 ] ; double xx [ 47 ] ; ( void ) mech ; ( void ) rtdvd ; ( void )
rtdvi ; ( void ) eqnEnableFlags ; ( void ) modeVector ; xx [ 0 ] = 2.0 ; xx [
1 ] = 1.0e-3 ; xx [ 2 ] = 0.5 ; xx [ 3 ] = xx [ 2 ] * state [ 6 ] ; xx [ 4 ]
= sin ( xx [ 3 ] ) ; xx [ 5 ] = xx [ 1 ] * xx [ 4 ] ; xx [ 6 ] = 2.0e-3 ; xx
[ 7 ] = xx [ 0 ] * xx [ 5 ] * xx [ 4 ] - xx [ 6 ] ; xx [ 8 ] = xx [ 2 ] *
state [ 4 ] ; xx [ 9 ] = sin ( xx [ 8 ] ) ; xx [ 10 ] = xx [ 1 ] * xx [ 9 ] ;
xx [ 11 ] = xx [ 0 ] * xx [ 10 ] * xx [ 9 ] - xx [ 6 ] ; xx [ 12 ] =
0.7071067811865476 ; xx [ 13 ] = xx [ 2 ] * state [ 2 ] ; xx [ 14 ] = xx [ 12
] * cos ( xx [ 13 ] ) ; xx [ 15 ] = xx [ 12 ] * sin ( xx [ 13 ] ) ; xx [ 13 ]
= xx [ 14 ] + xx [ 15 ] ; xx [ 16 ] = xx [ 12 ] * xx [ 13 ] ; xx [ 17 ] = xx
[ 14 ] - xx [ 15 ] ; xx [ 14 ] = xx [ 17 ] * xx [ 12 ] ; xx [ 15 ] = xx [ 16
] - xx [ 14 ] ; xx [ 18 ] = cos ( xx [ 8 ] ) ; xx [ 8 ] = xx [ 0 ] * xx [ 18
] * xx [ 10 ] ; xx [ 10 ] = xx [ 15 ] * xx [ 8 ] ; xx [ 19 ] = xx [ 16 ] + xx
[ 14 ] ; xx [ 14 ] = xx [ 11 ] * xx [ 15 ] ; xx [ 16 ] = xx [ 17 ] * xx [ 1 ]
; xx [ 20 ] = xx [ 0 ] * xx [ 17 ] * xx [ 16 ] - 2.500000000000001e-3 ; xx [
17 ] = xx [ 0 ] * xx [ 16 ] * xx [ 13 ] ; xx [ 13 ] = xx [ 12 ] * xx [ 12 ] *
xx [ 17 ] ; xx [ 16 ] = xx [ 12 ] * xx [ 20 ] * xx [ 12 ] ; xx [ 12 ] = xx [
9 ] * xx [ 19 ] + xx [ 15 ] * xx [ 18 ] ; xx [ 21 ] = cos ( xx [ 3 ] ) ; xx [
3 ] = xx [ 0 ] * xx [ 21 ] * xx [ 5 ] ; xx [ 5 ] = xx [ 12 ] * xx [ 3 ] ; xx
[ 22 ] = xx [ 15 ] * xx [ 9 ] - xx [ 18 ] * xx [ 19 ] ; xx [ 9 ] = xx [ 7 ] *
xx [ 12 ] ; xx [ 18 ] = xx [ 12 ] * xx [ 21 ] - xx [ 4 ] * xx [ 22 ] ; xx [
23 ] = xx [ 1 ] * xx [ 18 ] ; xx [ 24 ] = xx [ 2 ] * state [ 14 ] ; xx [ 25 ]
= sin ( xx [ 24 ] ) ; xx [ 26 ] = xx [ 1 ] * xx [ 25 ] ; xx [ 27 ] = xx [ 6 ]
- xx [ 0 ] * xx [ 26 ] * xx [ 25 ] ; xx [ 28 ] = xx [ 2 ] * state [ 12 ] ; xx
[ 29 ] = sin ( xx [ 28 ] ) ; xx [ 30 ] = xx [ 1 ] * xx [ 29 ] ; xx [ 31 ] =
xx [ 6 ] - xx [ 0 ] * xx [ 30 ] * xx [ 29 ] ; xx [ 32 ] = xx [ 2 ] * state [
8 ] ; xx [ 33 ] = cos ( xx [ 32 ] ) ; xx [ 34 ] = xx [ 2 ] * state [ 10 ] ;
xx [ 2 ] = sin ( xx [ 34 ] ) ; xx [ 35 ] = cos ( xx [ 34 ] ) ; xx [ 34 ] =
sin ( xx [ 32 ] ) ; xx [ 32 ] = xx [ 33 ] * xx [ 2 ] + xx [ 35 ] * xx [ 34 ]
; xx [ 36 ] = cos ( xx [ 28 ] ) ; xx [ 28 ] = xx [ 0 ] * xx [ 36 ] * xx [ 30
] ; xx [ 30 ] = xx [ 32 ] * xx [ 28 ] ; xx [ 37 ] = xx [ 33 ] * xx [ 35 ] -
xx [ 34 ] * xx [ 2 ] ; xx [ 38 ] = xx [ 32 ] * xx [ 31 ] ; xx [ 39 ] = xx [ 1
] * xx [ 2 ] ; xx [ 40 ] = xx [ 6 ] - xx [ 0 ] * xx [ 39 ] * xx [ 2 ] ; xx [
2 ] = xx [ 0 ] * xx [ 35 ] * xx [ 39 ] ; xx [ 6 ] = xx [ 2 ] * xx [ 34 ] ; xx
[ 35 ] = xx [ 40 ] * xx [ 34 ] ; xx [ 39 ] = xx [ 1 ] * xx [ 34 ] ; xx [ 41 ]
= xx [ 29 ] * xx [ 37 ] + xx [ 32 ] * xx [ 36 ] ; xx [ 42 ] = xx [ 41 ] * xx
[ 27 ] ; xx [ 43 ] = cos ( xx [ 24 ] ) ; xx [ 24 ] = xx [ 0 ] * xx [ 43 ] *
xx [ 26 ] ; xx [ 26 ] = xx [ 41 ] * xx [ 24 ] ; xx [ 44 ] = xx [ 32 ] * xx [
29 ] - xx [ 36 ] * xx [ 37 ] ; xx [ 29 ] = xx [ 25 ] * xx [ 44 ] - xx [ 41 ]
* xx [ 43 ] ; xx [ 36 ] = xx [ 1 ] * xx [ 29 ] ; xx [ 45 ] = fabs ( xx [ 7 ]
+ xx [ 11 ] + xx [ 0 ] * ( xx [ 10 ] * xx [ 19 ] - xx [ 15 ] * xx [ 14 ] ) +
xx [ 20 ] - xx [ 0 ] * ( xx [ 13 ] + xx [ 16 ] ) + state [ 0 ] - ( xx [ 5 ] *
xx [ 22 ] + xx [ 12 ] * xx [ 9 ] ) * xx [ 0 ] + xx [ 0 ] * xx [ 23 ] * xx [
18 ] - ( xx [ 27 ] + xx [ 31 ] - xx [ 0 ] * ( xx [ 30 ] * xx [ 37 ] + xx [ 32
] * xx [ 38 ] ) + xx [ 40 ] - xx [ 0 ] * ( xx [ 33 ] * xx [ 6 ] + xx [ 35 ] *
xx [ 34 ] ) - xx [ 0 ] * xx [ 39 ] * xx [ 34 ] - ( xx [ 41 ] * xx [ 42 ] - xx
[ 26 ] * xx [ 44 ] ) * xx [ 0 ] - xx [ 0 ] * xx [ 36 ] * xx [ 29 ] ) + 0.011
) ; xx [ 46 ] = fabs ( xx [ 3 ] + xx [ 0 ] * ( xx [ 9 ] * xx [ 22 ] - xx [ 12
] * xx [ 5 ] ) + xx [ 8 ] - ( xx [ 14 ] * xx [ 19 ] + xx [ 15 ] * xx [ 10 ] )
* xx [ 0 ] - ( xx [ 17 ] + ( xx [ 16 ] - xx [ 13 ] ) * xx [ 0 ] ) - xx [ 0 ]
* ( xx [ 21 ] * xx [ 22 ] + xx [ 12 ] * xx [ 4 ] ) * xx [ 23 ] - ( xx [ 0 ] *
( xx [ 43 ] * xx [ 44 ] + xx [ 41 ] * xx [ 25 ] ) * xx [ 36 ] + ( xx [ 38 ] *
xx [ 37 ] - xx [ 32 ] * xx [ 30 ] ) * xx [ 0 ] + ( xx [ 33 ] * xx [ 35 ] - xx
[ 6 ] * xx [ 34 ] ) * xx [ 0 ] + xx [ 2 ] + xx [ 0 ] * xx [ 33 ] * xx [ 39 ]
+ xx [ 28 ] - xx [ 0 ] * ( xx [ 41 ] * xx [ 26 ] + xx [ 42 ] * xx [ 44 ] ) +
xx [ 24 ] ) - 1.500000000000001e-3 ) ; ii [ 0 ] = 45 ; { int ll ; for ( ll =
46 ; ll < 47 ; ++ ll ) if ( xx [ ll ] > xx [ ii [ 0 ] ] ) ii [ 0 ] = ll ; }
ii [ 0 ] -= 45 ; xx [ 0 ] = xx [ 45 + ( ii [ 0 ] ) ] ; return xx [ 0 ] >
1.0e-9 ; } boolean_T
MagneticBasket_Simscape_Optimizer_dda62cd9_1_isVelocityViolation ( const void
* mech , const RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags
, const double * state , const int * modeVector ) { const double * rtdvd =
rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts . mValues ;
int ii [ 1 ] ; double xx [ 49 ] ; ( void ) mech ; ( void ) rtdvd ; ( void )
rtdvi ; ( void ) eqnEnableFlags ; ( void ) modeVector ; xx [ 0 ] = 2.0 ; xx [
1 ] = 0.7071067811865476 ; xx [ 2 ] = 1.0e-3 ; xx [ 3 ] = 0.5 ; xx [ 4 ] = xx
[ 3 ] * state [ 2 ] ; xx [ 5 ] = xx [ 1 ] * cos ( xx [ 4 ] ) ; xx [ 6 ] = xx
[ 1 ] * sin ( xx [ 4 ] ) ; xx [ 4 ] = xx [ 5 ] - xx [ 6 ] ; xx [ 7 ] = xx [ 4
] * xx [ 2 ] ; xx [ 8 ] = ( xx [ 2 ] - xx [ 0 ] * xx [ 4 ] * xx [ 7 ] ) *
state [ 3 ] ; xx [ 9 ] = xx [ 1 ] * xx [ 1 ] * xx [ 8 ] ; xx [ 10 ] = xx [ 5
] + xx [ 6 ] ; xx [ 5 ] = xx [ 0 ] * xx [ 7 ] * xx [ 10 ] * state [ 3 ] ; xx
[ 6 ] = xx [ 1 ] * xx [ 1 ] * xx [ 5 ] ; xx [ 7 ] = xx [ 3 ] * state [ 4 ] ;
xx [ 11 ] = cos ( xx [ 7 ] ) ; xx [ 12 ] = sin ( xx [ 7 ] ) ; xx [ 7 ] = xx [
2 ] * xx [ 12 ] ; xx [ 13 ] = 1.0e-3 ; xx [ 14 ] = xx [ 13 ] * xx [ 12 ] ; xx
[ 15 ] = xx [ 0 ] * xx [ 11 ] * xx [ 7 ] * state [ 5 ] + xx [ 0 ] * xx [ 11 ]
* xx [ 14 ] * state [ 3 ] ; xx [ 16 ] = xx [ 1 ] * xx [ 10 ] ; xx [ 10 ] = xx
[ 4 ] * xx [ 1 ] ; xx [ 1 ] = xx [ 16 ] - xx [ 10 ] ; xx [ 4 ] = 2.0e-3 ; xx
[ 17 ] = ( xx [ 2 ] - xx [ 0 ] * xx [ 7 ] * xx [ 12 ] ) * state [ 5 ] - ( xx
[ 0 ] * xx [ 14 ] * xx [ 12 ] - xx [ 4 ] ) * state [ 3 ] ; xx [ 7 ] = xx [ 1
] * xx [ 17 ] ; xx [ 14 ] = xx [ 16 ] + xx [ 10 ] ; xx [ 10 ] = xx [ 15 ] *
xx [ 1 ] ; xx [ 16 ] = xx [ 3 ] * state [ 6 ] ; xx [ 18 ] = cos ( xx [ 16 ] )
; xx [ 19 ] = sin ( xx [ 16 ] ) ; xx [ 16 ] = xx [ 2 ] * xx [ 19 ] ; xx [ 20
] = state [ 3 ] + state [ 5 ] ; xx [ 21 ] = xx [ 13 ] * xx [ 19 ] ; xx [ 22 ]
= xx [ 0 ] * xx [ 18 ] * xx [ 16 ] * state [ 7 ] + xx [ 20 ] * xx [ 0 ] * xx
[ 18 ] * xx [ 21 ] ; xx [ 23 ] = xx [ 12 ] * xx [ 14 ] + xx [ 1 ] * xx [ 11 ]
; xx [ 24 ] = ( xx [ 2 ] - xx [ 0 ] * xx [ 16 ] * xx [ 19 ] ) * state [ 7 ] -
( xx [ 0 ] * xx [ 21 ] * xx [ 19 ] - xx [ 4 ] ) * xx [ 20 ] ; xx [ 16 ] = xx
[ 23 ] * xx [ 24 ] ; xx [ 21 ] = xx [ 1 ] * xx [ 12 ] - xx [ 11 ] * xx [ 14 ]
; xx [ 11 ] = xx [ 22 ] * xx [ 23 ] ; xx [ 12 ] = ( xx [ 20 ] + state [ 7 ] )
* xx [ 13 ] ; xx [ 20 ] = xx [ 23 ] * xx [ 18 ] - xx [ 19 ] * xx [ 21 ] ; xx
[ 25 ] = xx [ 12 ] * xx [ 20 ] ; xx [ 26 ] = xx [ 3 ] * state [ 12 ] ; xx [
27 ] = sin ( xx [ 26 ] ) ; xx [ 28 ] = xx [ 3 ] * state [ 8 ] ; xx [ 29 ] =
cos ( xx [ 28 ] ) ; xx [ 30 ] = xx [ 3 ] * state [ 10 ] ; xx [ 31 ] = cos (
xx [ 30 ] ) ; xx [ 32 ] = sin ( xx [ 28 ] ) ; xx [ 28 ] = sin ( xx [ 30 ] ) ;
xx [ 30 ] = xx [ 29 ] * xx [ 31 ] - xx [ 32 ] * xx [ 28 ] ; xx [ 33 ] = xx [
29 ] * xx [ 28 ] + xx [ 31 ] * xx [ 32 ] ; xx [ 34 ] = cos ( xx [ 26 ] ) ; xx
[ 26 ] = xx [ 27 ] * xx [ 30 ] + xx [ 33 ] * xx [ 34 ] ; xx [ 35 ] = xx [ 3 ]
* state [ 14 ] ; xx [ 3 ] = sin ( xx [ 35 ] ) ; xx [ 36 ] = xx [ 2 ] * xx [ 3
] ; xx [ 37 ] = state [ 9 ] + state [ 11 ] ; xx [ 38 ] = xx [ 37 ] + state [
13 ] ; xx [ 39 ] = xx [ 13 ] * xx [ 3 ] ; xx [ 40 ] = ( xx [ 2 ] - xx [ 0 ] *
xx [ 36 ] * xx [ 3 ] ) * state [ 15 ] + xx [ 38 ] * ( xx [ 4 ] - xx [ 0 ] *
xx [ 39 ] * xx [ 3 ] ) ; xx [ 41 ] = xx [ 26 ] * xx [ 40 ] ; xx [ 42 ] = xx [
33 ] * xx [ 27 ] - xx [ 34 ] * xx [ 30 ] ; xx [ 43 ] = cos ( xx [ 35 ] ) ; xx
[ 35 ] = xx [ 38 ] * xx [ 0 ] * xx [ 43 ] * xx [ 39 ] + xx [ 0 ] * xx [ 43 ]
* xx [ 36 ] * state [ 15 ] ; xx [ 36 ] = xx [ 26 ] * xx [ 35 ] ; xx [ 39 ] =
xx [ 13 ] * xx [ 28 ] ; xx [ 44 ] = xx [ 2 ] * xx [ 28 ] ; xx [ 45 ] = xx [ 0
] * xx [ 31 ] * xx [ 39 ] * state [ 9 ] + xx [ 0 ] * xx [ 31 ] * xx [ 44 ] *
state [ 11 ] ; xx [ 31 ] = ( xx [ 2 ] - xx [ 0 ] * xx [ 44 ] * xx [ 28 ] ) *
state [ 11 ] + ( xx [ 4 ] - xx [ 0 ] * xx [ 39 ] * xx [ 28 ] ) * state [ 9 ]
; xx [ 28 ] = xx [ 31 ] * xx [ 32 ] ; xx [ 39 ] = xx [ 32 ] * xx [ 45 ] ; xx
[ 44 ] = xx [ 2 ] * xx [ 32 ] ; xx [ 46 ] = xx [ 13 ] * xx [ 27 ] ; xx [ 47 ]
= xx [ 2 ] * xx [ 27 ] ; xx [ 48 ] = xx [ 37 ] * xx [ 0 ] * xx [ 34 ] * xx [
46 ] + xx [ 0 ] * xx [ 34 ] * xx [ 47 ] * state [ 13 ] ; xx [ 34 ] = ( xx [ 2
] - xx [ 0 ] * xx [ 47 ] * xx [ 27 ] ) * state [ 13 ] + ( xx [ 4 ] - xx [ 0 ]
* xx [ 46 ] * xx [ 27 ] ) * xx [ 37 ] ; xx [ 4 ] = xx [ 33 ] * xx [ 34 ] ; xx
[ 27 ] = xx [ 33 ] * xx [ 48 ] ; xx [ 37 ] = ( xx [ 38 ] + state [ 15 ] ) *
xx [ 13 ] ; xx [ 13 ] = xx [ 3 ] * xx [ 42 ] - xx [ 26 ] * xx [ 43 ] ; xx [
38 ] = xx [ 37 ] * xx [ 13 ] ; xx [ 46 ] = fabs ( state [ 1 ] + xx [ 0 ] * (
xx [ 9 ] + xx [ 6 ] ) - xx [ 5 ] + xx [ 15 ] + xx [ 0 ] * ( xx [ 7 ] * xx [
14 ] - xx [ 1 ] * xx [ 10 ] ) + xx [ 22 ] - ( xx [ 16 ] * xx [ 21 ] + xx [ 23
] * xx [ 11 ] ) * xx [ 0 ] - xx [ 0 ] * ( xx [ 18 ] * xx [ 21 ] + xx [ 23 ] *
xx [ 19 ] ) * xx [ 25 ] - ( xx [ 0 ] * ( xx [ 41 ] * xx [ 42 ] + xx [ 26 ] *
xx [ 36 ] ) - xx [ 35 ] - ( xx [ 45 ] + ( xx [ 29 ] * xx [ 28 ] - xx [ 39 ] *
xx [ 32 ] ) * xx [ 0 ] + xx [ 0 ] * xx [ 29 ] * xx [ 44 ] * state [ 9 ] + xx
[ 48 ] + ( xx [ 4 ] * xx [ 30 ] - xx [ 33 ] * xx [ 27 ] ) * xx [ 0 ] ) - xx [
0 ] * ( xx [ 43 ] * xx [ 42 ] + xx [ 26 ] * xx [ 3 ] ) * xx [ 38 ] ) ) ; xx [
47 ] = fabs ( xx [ 12 ] - xx [ 0 ] * xx [ 25 ] * xx [ 20 ] + xx [ 8 ] - ( xx
[ 9 ] - xx [ 6 ] ) * xx [ 0 ] + xx [ 17 ] - ( xx [ 10 ] * xx [ 14 ] + xx [ 1
] * xx [ 7 ] ) * xx [ 0 ] + xx [ 24 ] + xx [ 0 ] * ( xx [ 11 ] * xx [ 21 ] -
xx [ 23 ] * xx [ 16 ] ) - ( xx [ 37 ] - xx [ 0 ] * xx [ 38 ] * xx [ 13 ] + (
xx [ 2 ] - xx [ 0 ] * xx [ 44 ] * xx [ 32 ] ) * state [ 9 ] + xx [ 31 ] - xx
[ 0 ] * ( xx [ 29 ] * xx [ 39 ] + xx [ 28 ] * xx [ 32 ] ) + xx [ 34 ] - xx [
0 ] * ( xx [ 27 ] * xx [ 30 ] + xx [ 33 ] * xx [ 4 ] ) + xx [ 40 ] - ( xx [
26 ] * xx [ 41 ] - xx [ 36 ] * xx [ 42 ] ) * xx [ 0 ] ) ) ; ii [ 0 ] = 46 ; {
int ll ; for ( ll = 47 ; ll < 48 ; ++ ll ) if ( xx [ ll ] > xx [ ii [ 0 ] ] )
ii [ 0 ] = ll ; } ii [ 0 ] -= 46 ; xx [ 0 ] = xx [ 46 + ( ii [ 0 ] ) ] ;
return xx [ 0 ] > 1.0e-9 ; } PmfMessageId
MagneticBasket_Simscape_Optimizer_dda62cd9_1_projectStateSim ( const void *
mech , const RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags ,
const int * modeVector , double * state , void * neDiagMgr0 ) { const double
* rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts .
mValues ; NeuDiagnosticManager * neDiagMgr = ( NeuDiagnosticManager * )
neDiagMgr0 ; int ii [ 2 ] ; double xx [ 154 ] ; ( void ) mech ; ( void )
rtdvd ; ( void ) rtdvi ; ( void ) eqnEnableFlags ; ( void ) modeVector ; (
void ) neDiagMgr ; xx [ 0 ] = 1.0 ; xx [ 1 ] = 2.0 ; xx [ 2 ] = 0.5 ; xx [ 3
] = xx [ 2 ] * state [ 4 ] ; xx [ 4 ] = cos ( xx [ 3 ] ) ; xx [ 5 ] = xx [ 2
] * state [ 6 ] ; xx [ 6 ] = sin ( xx [ 5 ] ) ; xx [ 7 ] = cos ( xx [ 5 ] ) ;
xx [ 5 ] = sin ( xx [ 3 ] ) ; xx [ 3 ] = xx [ 4 ] * xx [ 6 ] + xx [ 7 ] * xx
[ 5 ] ; xx [ 8 ] = 0.7071067811865476 ; xx [ 9 ] = xx [ 2 ] * state [ 2 ] ;
xx [ 10 ] = xx [ 8 ] * cos ( xx [ 9 ] ) ; xx [ 11 ] = xx [ 8 ] * sin ( xx [ 9
] ) ; xx [ 9 ] = xx [ 10 ] + xx [ 11 ] ; xx [ 12 ] = xx [ 8 ] * xx [ 9 ] ; xx
[ 13 ] = xx [ 10 ] - xx [ 11 ] ; xx [ 10 ] = xx [ 13 ] * xx [ 8 ] ; xx [ 11 ]
= xx [ 12 ] + xx [ 10 ] ; xx [ 14 ] = xx [ 12 ] - xx [ 10 ] ; xx [ 10 ] = xx
[ 4 ] * xx [ 7 ] - xx [ 5 ] * xx [ 6 ] ; xx [ 12 ] = xx [ 3 ] * xx [ 11 ] +
xx [ 14 ] * xx [ 10 ] ; xx [ 15 ] = 1.0e-3 ; xx [ 16 ] = xx [ 12 ] * xx [ 15
] ; xx [ 17 ] = xx [ 15 ] * xx [ 6 ] ; xx [ 18 ] = xx [ 1 ] * xx [ 7 ] * xx [
17 ] ; xx [ 19 ] = 2.0e-3 ; xx [ 20 ] = xx [ 1 ] * xx [ 17 ] * xx [ 6 ] - xx
[ 19 ] ; xx [ 17 ] = xx [ 20 ] * xx [ 5 ] ; xx [ 21 ] = xx [ 18 ] * xx [ 5 ]
; xx [ 22 ] = xx [ 18 ] - ( xx [ 4 ] * xx [ 17 ] + xx [ 21 ] * xx [ 5 ] ) *
xx [ 1 ] ; xx [ 23 ] = xx [ 15 ] * xx [ 5 ] ; xx [ 24 ] = xx [ 1 ] * xx [ 4 ]
* xx [ 23 ] ; xx [ 25 ] = xx [ 22 ] + xx [ 24 ] ; xx [ 26 ] = xx [ 13 ] * xx
[ 25 ] ; xx [ 27 ] = xx [ 4 ] * xx [ 21 ] ; xx [ 21 ] = xx [ 17 ] * xx [ 5 ]
; xx [ 17 ] = xx [ 1 ] * xx [ 23 ] * xx [ 5 ] - xx [ 19 ] ; xx [ 23 ] = xx [
20 ] + xx [ 1 ] * ( xx [ 27 ] - xx [ 21 ] ) + xx [ 17 ] ; xx [ 28 ] = xx [ 23
] * xx [ 13 ] ; xx [ 29 ] = 1.0e-3 ; xx [ 30 ] = xx [ 13 ] * xx [ 29 ] ; xx [
31 ] = xx [ 25 ] - ( xx [ 13 ] * xx [ 26 ] - xx [ 28 ] * xx [ 9 ] ) * xx [ 1
] - xx [ 1 ] * xx [ 30 ] * xx [ 9 ] ; xx [ 25 ] = xx [ 1 ] * ( xx [ 13 ] * xx
[ 28 ] + xx [ 26 ] * xx [ 9 ] ) - ( xx [ 23 ] + xx [ 1 ] * xx [ 13 ] * xx [
30 ] ) + xx [ 29 ] ; xx [ 23 ] = xx [ 8 ] * xx [ 25 ] * xx [ 8 ] ; xx [ 26 ]
= xx [ 8 ] * xx [ 31 ] * xx [ 8 ] ; xx [ 28 ] = xx [ 29 ] * xx [ 5 ] ; xx [
30 ] = xx [ 22 ] + xx [ 1 ] * xx [ 4 ] * xx [ 28 ] ; xx [ 22 ] = xx [ 1 ] * (
xx [ 21 ] - xx [ 27 ] ) - ( xx [ 20 ] + xx [ 1 ] * xx [ 28 ] * xx [ 5 ] ) +
xx [ 29 ] ; xx [ 21 ] = xx [ 14 ] * xx [ 22 ] ; xx [ 27 ] = xx [ 14 ] * xx [
30 ] ; xx [ 28 ] = xx [ 14 ] * xx [ 5 ] - xx [ 4 ] * xx [ 11 ] ; xx [ 32 ] =
xx [ 5 ] * xx [ 11 ] + xx [ 14 ] * xx [ 4 ] ; xx [ 4 ] = xx [ 32 ] * xx [ 7 ]
- xx [ 6 ] * xx [ 28 ] ; xx [ 5 ] = xx [ 15 ] * xx [ 4 ] ; xx [ 33 ] = xx [ 1
] * ( xx [ 7 ] * xx [ 28 ] + xx [ 32 ] * xx [ 6 ] ) * xx [ 5 ] ; xx [ 34 ] =
xx [ 29 ] * xx [ 6 ] ; xx [ 35 ] = xx [ 1 ] * xx [ 7 ] * xx [ 34 ] ; xx [ 7 ]
= xx [ 29 ] - xx [ 1 ] * xx [ 34 ] * xx [ 6 ] ; xx [ 6 ] = xx [ 32 ] * xx [ 7
] ; xx [ 34 ] = xx [ 32 ] * xx [ 35 ] ; xx [ 36 ] = xx [ 2 ] * state [ 8 ] ;
xx [ 37 ] = cos ( xx [ 36 ] ) ; xx [ 38 ] = xx [ 2 ] * state [ 12 ] ; xx [ 39
] = cos ( xx [ 38 ] ) ; xx [ 40 ] = xx [ 2 ] * state [ 14 ] ; xx [ 41 ] = sin
( xx [ 40 ] ) ; xx [ 42 ] = cos ( xx [ 40 ] ) ; xx [ 40 ] = sin ( xx [ 38 ] )
; xx [ 38 ] = xx [ 39 ] * xx [ 41 ] + xx [ 42 ] * xx [ 40 ] ; xx [ 43 ] = xx
[ 2 ] * state [ 10 ] ; xx [ 44 ] = sin ( xx [ 43 ] ) ; xx [ 45 ] = cos ( xx [
43 ] ) ; xx [ 43 ] = xx [ 39 ] * xx [ 42 ] - xx [ 40 ] * xx [ 41 ] ; xx [ 46
] = xx [ 38 ] * xx [ 44 ] - xx [ 45 ] * xx [ 43 ] ; xx [ 47 ] = xx [ 38 ] *
xx [ 45 ] + xx [ 44 ] * xx [ 43 ] ; xx [ 48 ] = sin ( xx [ 36 ] ) ; xx [ 36 ]
= xx [ 48 ] * xx [ 46 ] - xx [ 47 ] * xx [ 37 ] ; xx [ 49 ] = xx [ 15 ] * xx
[ 36 ] ; xx [ 50 ] = xx [ 15 ] * xx [ 41 ] ; xx [ 51 ] = xx [ 19 ] - xx [ 1 ]
* xx [ 50 ] * xx [ 41 ] ; xx [ 52 ] = xx [ 51 ] * xx [ 40 ] ; xx [ 53 ] = xx
[ 1 ] * xx [ 42 ] * xx [ 50 ] ; xx [ 50 ] = xx [ 53 ] * xx [ 40 ] ; xx [ 54 ]
= ( xx [ 39 ] * xx [ 52 ] - xx [ 50 ] * xx [ 40 ] ) * xx [ 1 ] ; xx [ 55 ] =
xx [ 15 ] * xx [ 40 ] ; xx [ 56 ] = xx [ 1 ] * xx [ 39 ] * xx [ 55 ] ; xx [
57 ] = xx [ 54 ] + xx [ 53 ] + xx [ 56 ] ; xx [ 58 ] = xx [ 51 ] - xx [ 1 ] *
( xx [ 39 ] * xx [ 50 ] + xx [ 52 ] * xx [ 40 ] ) ; xx [ 50 ] = xx [ 19 ] -
xx [ 1 ] * xx [ 55 ] * xx [ 40 ] ; xx [ 52 ] = xx [ 58 ] + xx [ 50 ] ; xx [
55 ] = xx [ 52 ] * xx [ 44 ] ; xx [ 59 ] = xx [ 45 ] * xx [ 55 ] ; xx [ 60 ]
= xx [ 44 ] * xx [ 57 ] ; xx [ 61 ] = xx [ 60 ] * xx [ 44 ] ; xx [ 62 ] = xx
[ 15 ] * xx [ 44 ] ; xx [ 63 ] = xx [ 1 ] * xx [ 45 ] * xx [ 62 ] ; xx [ 64 ]
= xx [ 57 ] + xx [ 1 ] * ( xx [ 59 ] - xx [ 61 ] ) + xx [ 63 ] ; xx [ 65 ] =
xx [ 48 ] * xx [ 64 ] ; xx [ 66 ] = ( xx [ 45 ] * xx [ 60 ] + xx [ 55 ] * xx
[ 44 ] ) * xx [ 1 ] ; xx [ 55 ] = xx [ 19 ] - xx [ 1 ] * xx [ 62 ] * xx [ 44
] ; xx [ 60 ] = xx [ 52 ] - xx [ 66 ] + xx [ 55 ] ; xx [ 62 ] = xx [ 60 ] *
xx [ 48 ] ; xx [ 67 ] = xx [ 29 ] * xx [ 48 ] ; xx [ 68 ] = xx [ 37 ] * xx [
45 ] - xx [ 48 ] * xx [ 44 ] ; xx [ 69 ] = xx [ 37 ] * xx [ 44 ] + xx [ 45 ]
* xx [ 48 ] ; xx [ 70 ] = xx [ 38 ] * xx [ 68 ] + xx [ 69 ] * xx [ 43 ] ; xx
[ 71 ] = xx [ 70 ] * xx [ 15 ] ; xx [ 72 ] = xx [ 29 ] * xx [ 44 ] ; xx [ 73
] = xx [ 1 ] * ( xx [ 61 ] - xx [ 59 ] ) - xx [ 57 ] - xx [ 1 ] * xx [ 45 ] *
xx [ 72 ] ; xx [ 45 ] = xx [ 52 ] - ( xx [ 66 ] + xx [ 1 ] * xx [ 72 ] * xx [
44 ] ) + xx [ 29 ] ; xx [ 44 ] = xx [ 45 ] * xx [ 48 ] ; xx [ 52 ] = xx [ 48
] * xx [ 73 ] ; xx [ 57 ] = xx [ 69 ] * xx [ 40 ] - xx [ 39 ] * xx [ 68 ] ;
xx [ 59 ] = xx [ 40 ] * xx [ 68 ] + xx [ 69 ] * xx [ 39 ] ; xx [ 61 ] = xx [
41 ] * xx [ 57 ] - xx [ 59 ] * xx [ 42 ] ; xx [ 66 ] = xx [ 15 ] * xx [ 61 ]
; xx [ 72 ] = xx [ 1 ] * ( xx [ 42 ] * xx [ 57 ] + xx [ 59 ] * xx [ 41 ] ) *
xx [ 66 ] ; xx [ 74 ] = xx [ 29 ] * xx [ 40 ] ; xx [ 75 ] = xx [ 53 ] + xx [
54 ] + xx [ 1 ] * xx [ 39 ] * xx [ 74 ] ; xx [ 39 ] = xx [ 58 ] - xx [ 1 ] *
xx [ 74 ] * xx [ 40 ] + xx [ 29 ] ; xx [ 40 ] = xx [ 69 ] * xx [ 39 ] ; xx [
54 ] = xx [ 69 ] * xx [ 75 ] ; xx [ 58 ] = xx [ 29 ] * xx [ 41 ] ; xx [ 74 ]
= xx [ 29 ] - xx [ 1 ] * xx [ 58 ] * xx [ 41 ] ; xx [ 41 ] = xx [ 59 ] * xx [
74 ] ; xx [ 76 ] = xx [ 1 ] * xx [ 42 ] * xx [ 58 ] ; xx [ 42 ] = xx [ 59 ] *
xx [ 76 ] ; xx [ 58 ] = 0.0 ; xx [ 77 ] = xx [ 1 ] * xx [ 5 ] * xx [ 4 ] ; xx
[ 4 ] = xx [ 1 ] * xx [ 66 ] * xx [ 61 ] ; xx [ 78 ] = xx [ 0 ] ; xx [ 79 ] =
xx [ 1 ] * xx [ 16 ] * ( xx [ 10 ] * xx [ 11 ] - xx [ 3 ] * xx [ 14 ] ) + xx
[ 31 ] + xx [ 1 ] * ( xx [ 23 ] - xx [ 26 ] ) ; xx [ 80 ] = xx [ 30 ] + xx [
1 ] * ( xx [ 21 ] * xx [ 11 ] - xx [ 14 ] * xx [ 27 ] ) - xx [ 33 ] ; xx [ 81
] = xx [ 35 ] - ( xx [ 6 ] * xx [ 28 ] + xx [ 32 ] * xx [ 34 ] ) * xx [ 1 ] -
xx [ 33 ] ; xx [ 82 ] = xx [ 1 ] * ( xx [ 37 ] * xx [ 46 ] + xx [ 47 ] * xx [
48 ] ) * xx [ 49 ] - ( xx [ 1 ] * ( xx [ 65 ] * xx [ 48 ] - xx [ 37 ] * xx [
62 ] ) - xx [ 64 ] - xx [ 1 ] * xx [ 37 ] * xx [ 67 ] ) ; xx [ 83 ] = xx [ 1
] * xx [ 71 ] * ( xx [ 68 ] * xx [ 43 ] - xx [ 69 ] * xx [ 38 ] ) - ( xx [ 73
] - ( xx [ 37 ] * xx [ 44 ] + xx [ 52 ] * xx [ 48 ] ) * xx [ 1 ] ) ; xx [ 84
] = xx [ 72 ] + xx [ 75 ] + ( xx [ 40 ] * xx [ 68 ] - xx [ 69 ] * xx [ 54 ] )
* xx [ 1 ] ; xx [ 85 ] = xx [ 72 ] - ( ( xx [ 41 ] * xx [ 57 ] + xx [ 59 ] *
xx [ 42 ] ) * xx [ 1 ] - xx [ 76 ] ) ; xx [ 86 ] = xx [ 58 ] ; xx [ 87 ] = xx
[ 25 ] - ( xx [ 1 ] * xx [ 12 ] * xx [ 16 ] + ( xx [ 26 ] + xx [ 23 ] ) * xx
[ 1 ] ) + xx [ 15 ] ; xx [ 88 ] = xx [ 22 ] - ( xx [ 77 ] + ( xx [ 27 ] * xx
[ 11 ] + xx [ 14 ] * xx [ 21 ] ) * xx [ 1 ] ) + xx [ 15 ] ; xx [ 89 ] = xx [
7 ] + xx [ 1 ] * ( xx [ 34 ] * xx [ 28 ] - xx [ 32 ] * xx [ 6 ] ) - xx [ 77 ]
+ xx [ 15 ] ; xx [ 90 ] = - ( xx [ 60 ] - ( ( xx [ 37 ] * xx [ 65 ] + xx [ 62
] * xx [ 48 ] ) * xx [ 1 ] + xx [ 1 ] * xx [ 67 ] * xx [ 48 ] ) - xx [ 1 ] *
xx [ 49 ] * xx [ 36 ] + xx [ 19 ] ) ; xx [ 91 ] = - ( xx [ 45 ] + xx [ 1 ] *
( xx [ 37 ] * xx [ 52 ] - xx [ 44 ] * xx [ 48 ] ) - xx [ 1 ] * xx [ 70 ] * xx
[ 71 ] + xx [ 15 ] ) ; xx [ 92 ] = - ( xx [ 39 ] - xx [ 1 ] * ( xx [ 54 ] *
xx [ 68 ] + xx [ 69 ] * xx [ 40 ] ) - xx [ 4 ] + xx [ 15 ] ) ; xx [ 93 ] = -
( xx [ 74 ] + xx [ 1 ] * ( xx [ 42 ] * xx [ 57 ] - xx [ 59 ] * xx [ 41 ] ) -
xx [ 4 ] + xx [ 15 ] ) ; xx [ 3 ] = xx [ 14 ] * xx [ 24 ] ; xx [ 5 ] = xx [
17 ] * xx [ 14 ] ; xx [ 6 ] = xx [ 13 ] * xx [ 15 ] ; xx [ 7 ] =
2.500000000000001e-3 ; xx [ 10 ] = xx [ 1 ] * xx [ 13 ] * xx [ 6 ] - xx [ 7 ]
; xx [ 12 ] = xx [ 1 ] * xx [ 6 ] * xx [ 9 ] ; xx [ 6 ] = xx [ 8 ] * xx [ 8 ]
* xx [ 12 ] ; xx [ 9 ] = xx [ 8 ] * xx [ 10 ] * xx [ 8 ] ; xx [ 13 ] = xx [
32 ] * xx [ 18 ] ; xx [ 16 ] = xx [ 20 ] * xx [ 32 ] ; xx [ 21 ] = xx [ 69 ]
* xx [ 56 ] ; xx [ 22 ] = xx [ 69 ] * xx [ 50 ] ; xx [ 23 ] = xx [ 63 ] * xx
[ 48 ] ; xx [ 25 ] = xx [ 55 ] * xx [ 48 ] ; xx [ 26 ] = xx [ 15 ] * xx [ 48
] ; xx [ 27 ] = xx [ 59 ] * xx [ 51 ] ; xx [ 30 ] = xx [ 59 ] * xx [ 53 ] ;
xx [ 31 ] = 0.011 ; xx [ 34 ] = 1.500000000000001e-3 ; xx [ 35 ] = - ( xx [
20 ] + xx [ 17 ] + xx [ 1 ] * ( xx [ 3 ] * xx [ 11 ] - xx [ 14 ] * xx [ 5 ] )
+ xx [ 10 ] - xx [ 1 ] * ( xx [ 6 ] + xx [ 9 ] ) + state [ 0 ] - ( xx [ 13 ]
* xx [ 28 ] + xx [ 32 ] * xx [ 16 ] ) * xx [ 1 ] + xx [ 77 ] - ( xx [ 51 ] +
xx [ 50 ] - xx [ 1 ] * ( xx [ 21 ] * xx [ 68 ] + xx [ 69 ] * xx [ 22 ] ) + xx
[ 55 ] - xx [ 1 ] * ( xx [ 37 ] * xx [ 23 ] + xx [ 25 ] * xx [ 48 ] ) - xx [
1 ] * xx [ 26 ] * xx [ 48 ] - ( xx [ 59 ] * xx [ 27 ] - xx [ 30 ] * xx [ 57 ]
) * xx [ 1 ] - xx [ 4 ] ) + xx [ 31 ] ) ; xx [ 36 ] = - ( xx [ 18 ] + xx [ 1
] * ( xx [ 16 ] * xx [ 28 ] - xx [ 32 ] * xx [ 13 ] ) + xx [ 24 ] - ( xx [ 5
] * xx [ 11 ] + xx [ 14 ] * xx [ 3 ] ) * xx [ 1 ] - ( xx [ 12 ] + ( xx [ 9 ]
- xx [ 6 ] ) * xx [ 1 ] ) - xx [ 33 ] - ( xx [ 72 ] + ( xx [ 22 ] * xx [ 68 ]
- xx [ 69 ] * xx [ 21 ] ) * xx [ 1 ] + ( xx [ 37 ] * xx [ 25 ] - xx [ 23 ] *
xx [ 48 ] ) * xx [ 1 ] + xx [ 63 ] + xx [ 1 ] * xx [ 37 ] * xx [ 26 ] + xx [
56 ] - xx [ 1 ] * ( xx [ 59 ] * xx [ 30 ] + xx [ 27 ] * xx [ 57 ] ) + xx [ 53
] ) - xx [ 34 ] ) ; xx [ 3 ] = 1.0e-8 ; memcpy ( xx + 37 , xx + 78 , 16 *
sizeof ( double ) ) ; factorAndSolveWide ( 2 , 8 , xx + 37 , xx + 4 , xx + 9
, ii + 0 , xx + 35 , xx [ 3 ] , xx + 20 ) ; xx [ 4 ] = state [ 0 ] + xx [ 20
] ; xx [ 5 ] = state [ 4 ] + xx [ 22 ] ; xx [ 6 ] = xx [ 5 ] * xx [ 2 ] ; xx
[ 9 ] = cos ( xx [ 6 ] ) ; xx [ 10 ] = state [ 6 ] + xx [ 23 ] ; xx [ 11 ] =
xx [ 10 ] * xx [ 2 ] ; xx [ 12 ] = sin ( xx [ 11 ] ) ; xx [ 13 ] = cos ( xx [
11 ] ) ; xx [ 11 ] = sin ( xx [ 6 ] ) ; xx [ 6 ] = xx [ 9 ] * xx [ 12 ] + xx
[ 13 ] * xx [ 11 ] ; xx [ 14 ] = state [ 2 ] + xx [ 21 ] ; xx [ 16 ] = xx [
14 ] * xx [ 2 ] ; xx [ 17 ] = xx [ 8 ] * cos ( xx [ 16 ] ) ; xx [ 18 ] = xx [
8 ] * sin ( xx [ 16 ] ) ; xx [ 16 ] = xx [ 17 ] + xx [ 18 ] ; xx [ 28 ] = xx
[ 8 ] * xx [ 16 ] ; xx [ 30 ] = xx [ 17 ] - xx [ 18 ] ; xx [ 17 ] = xx [ 30 ]
* xx [ 8 ] ; xx [ 18 ] = xx [ 28 ] + xx [ 17 ] ; xx [ 32 ] = xx [ 28 ] - xx [
17 ] ; xx [ 17 ] = xx [ 9 ] * xx [ 13 ] - xx [ 11 ] * xx [ 12 ] ; xx [ 28 ] =
xx [ 6 ] * xx [ 18 ] + xx [ 32 ] * xx [ 17 ] ; xx [ 33 ] = xx [ 28 ] * xx [
15 ] ; xx [ 35 ] = xx [ 15 ] * xx [ 12 ] ; xx [ 36 ] = xx [ 1 ] * xx [ 13 ] *
xx [ 35 ] ; xx [ 37 ] = xx [ 1 ] * xx [ 35 ] * xx [ 12 ] - xx [ 19 ] ; xx [
35 ] = xx [ 37 ] * xx [ 11 ] ; xx [ 38 ] = xx [ 36 ] * xx [ 11 ] ; xx [ 39 ]
= xx [ 36 ] - ( xx [ 9 ] * xx [ 35 ] + xx [ 38 ] * xx [ 11 ] ) * xx [ 1 ] ;
xx [ 40 ] = xx [ 15 ] * xx [ 11 ] ; xx [ 41 ] = xx [ 1 ] * xx [ 9 ] * xx [ 40
] ; xx [ 42 ] = xx [ 39 ] + xx [ 41 ] ; xx [ 43 ] = xx [ 30 ] * xx [ 42 ] ;
xx [ 44 ] = xx [ 9 ] * xx [ 38 ] ; xx [ 38 ] = xx [ 35 ] * xx [ 11 ] ; xx [
35 ] = xx [ 1 ] * xx [ 40 ] * xx [ 11 ] - xx [ 19 ] ; xx [ 40 ] = xx [ 37 ] +
xx [ 1 ] * ( xx [ 44 ] - xx [ 38 ] ) + xx [ 35 ] ; xx [ 45 ] = xx [ 40 ] * xx
[ 30 ] ; xx [ 46 ] = xx [ 30 ] * xx [ 29 ] ; xx [ 47 ] = xx [ 42 ] - ( xx [
30 ] * xx [ 43 ] - xx [ 45 ] * xx [ 16 ] ) * xx [ 1 ] - xx [ 1 ] * xx [ 46 ]
* xx [ 16 ] ; xx [ 42 ] = xx [ 1 ] * ( xx [ 30 ] * xx [ 45 ] + xx [ 43 ] * xx
[ 16 ] ) - ( xx [ 40 ] + xx [ 1 ] * xx [ 30 ] * xx [ 46 ] ) + xx [ 29 ] ; xx
[ 40 ] = xx [ 8 ] * xx [ 42 ] * xx [ 8 ] ; xx [ 43 ] = xx [ 8 ] * xx [ 47 ] *
xx [ 8 ] ; xx [ 45 ] = xx [ 29 ] * xx [ 11 ] ; xx [ 46 ] = xx [ 39 ] + xx [ 1
] * xx [ 9 ] * xx [ 45 ] ; xx [ 39 ] = xx [ 1 ] * ( xx [ 38 ] - xx [ 44 ] ) -
( xx [ 37 ] + xx [ 1 ] * xx [ 45 ] * xx [ 11 ] ) + xx [ 29 ] ; xx [ 38 ] = xx
[ 32 ] * xx [ 39 ] ; xx [ 44 ] = xx [ 32 ] * xx [ 46 ] ; xx [ 45 ] = xx [ 32
] * xx [ 11 ] - xx [ 9 ] * xx [ 18 ] ; xx [ 48 ] = xx [ 11 ] * xx [ 18 ] + xx
[ 32 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 48 ] * xx [ 13 ] - xx [ 12 ] * xx [ 45 ]
; xx [ 11 ] = xx [ 15 ] * xx [ 9 ] ; xx [ 49 ] = xx [ 1 ] * ( xx [ 13 ] * xx
[ 45 ] + xx [ 48 ] * xx [ 12 ] ) * xx [ 11 ] ; xx [ 50 ] = xx [ 29 ] * xx [
12 ] ; xx [ 51 ] = xx [ 1 ] * xx [ 13 ] * xx [ 50 ] ; xx [ 13 ] = xx [ 29 ] -
xx [ 1 ] * xx [ 50 ] * xx [ 12 ] ; xx [ 12 ] = xx [ 48 ] * xx [ 13 ] ; xx [
50 ] = xx [ 48 ] * xx [ 51 ] ; xx [ 52 ] = state [ 8 ] + xx [ 24 ] ; xx [ 53
] = xx [ 52 ] * xx [ 2 ] ; xx [ 54 ] = cos ( xx [ 53 ] ) ; xx [ 55 ] = state
[ 12 ] + xx [ 26 ] ; xx [ 56 ] = xx [ 55 ] * xx [ 2 ] ; xx [ 57 ] = cos ( xx
[ 56 ] ) ; xx [ 59 ] = state [ 14 ] + xx [ 27 ] ; xx [ 60 ] = xx [ 59 ] * xx
[ 2 ] ; xx [ 61 ] = sin ( xx [ 60 ] ) ; xx [ 62 ] = cos ( xx [ 60 ] ) ; xx [
60 ] = sin ( xx [ 56 ] ) ; xx [ 56 ] = xx [ 57 ] * xx [ 61 ] + xx [ 62 ] * xx
[ 60 ] ; xx [ 20 ] = state [ 10 ] + xx [ 25 ] ; xx [ 21 ] = xx [ 20 ] * xx [
2 ] ; xx [ 22 ] = sin ( xx [ 21 ] ) ; xx [ 23 ] = cos ( xx [ 21 ] ) ; xx [ 21
] = xx [ 57 ] * xx [ 62 ] - xx [ 60 ] * xx [ 61 ] ; xx [ 24 ] = xx [ 56 ] *
xx [ 22 ] - xx [ 23 ] * xx [ 21 ] ; xx [ 25 ] = xx [ 56 ] * xx [ 23 ] + xx [
22 ] * xx [ 21 ] ; xx [ 26 ] = sin ( xx [ 53 ] ) ; xx [ 27 ] = xx [ 26 ] * xx
[ 24 ] - xx [ 25 ] * xx [ 54 ] ; xx [ 53 ] = xx [ 15 ] * xx [ 27 ] ; xx [ 63
] = xx [ 15 ] * xx [ 61 ] ; xx [ 64 ] = xx [ 19 ] - xx [ 1 ] * xx [ 63 ] * xx
[ 61 ] ; xx [ 65 ] = xx [ 64 ] * xx [ 60 ] ; xx [ 66 ] = xx [ 1 ] * xx [ 62 ]
* xx [ 63 ] ; xx [ 63 ] = xx [ 66 ] * xx [ 60 ] ; xx [ 67 ] = ( xx [ 57 ] *
xx [ 65 ] - xx [ 63 ] * xx [ 60 ] ) * xx [ 1 ] ; xx [ 68 ] = xx [ 15 ] * xx [
60 ] ; xx [ 69 ] = xx [ 1 ] * xx [ 57 ] * xx [ 68 ] ; xx [ 70 ] = xx [ 67 ] +
xx [ 66 ] + xx [ 69 ] ; xx [ 71 ] = xx [ 64 ] - xx [ 1 ] * ( xx [ 57 ] * xx [
63 ] + xx [ 65 ] * xx [ 60 ] ) ; xx [ 63 ] = xx [ 19 ] - xx [ 1 ] * xx [ 68 ]
* xx [ 60 ] ; xx [ 65 ] = xx [ 71 ] + xx [ 63 ] ; xx [ 68 ] = xx [ 65 ] * xx
[ 22 ] ; xx [ 72 ] = xx [ 23 ] * xx [ 68 ] ; xx [ 73 ] = xx [ 22 ] * xx [ 70
] ; xx [ 74 ] = xx [ 73 ] * xx [ 22 ] ; xx [ 75 ] = xx [ 15 ] * xx [ 22 ] ;
xx [ 76 ] = xx [ 1 ] * xx [ 23 ] * xx [ 75 ] ; xx [ 77 ] = xx [ 70 ] + xx [ 1
] * ( xx [ 72 ] - xx [ 74 ] ) + xx [ 76 ] ; xx [ 78 ] = xx [ 26 ] * xx [ 77 ]
; xx [ 79 ] = ( xx [ 23 ] * xx [ 73 ] + xx [ 68 ] * xx [ 22 ] ) * xx [ 1 ] ;
xx [ 68 ] = xx [ 19 ] - xx [ 1 ] * xx [ 75 ] * xx [ 22 ] ; xx [ 73 ] = xx [
65 ] - xx [ 79 ] + xx [ 68 ] ; xx [ 75 ] = xx [ 73 ] * xx [ 26 ] ; xx [ 80 ]
= xx [ 29 ] * xx [ 26 ] ; xx [ 81 ] = xx [ 54 ] * xx [ 23 ] - xx [ 26 ] * xx
[ 22 ] ; xx [ 82 ] = xx [ 54 ] * xx [ 22 ] + xx [ 23 ] * xx [ 26 ] ; xx [ 83
] = xx [ 56 ] * xx [ 81 ] + xx [ 82 ] * xx [ 21 ] ; xx [ 84 ] = xx [ 83 ] *
xx [ 15 ] ; xx [ 85 ] = xx [ 29 ] * xx [ 22 ] ; xx [ 86 ] = xx [ 1 ] * ( xx [
74 ] - xx [ 72 ] ) - xx [ 70 ] - xx [ 1 ] * xx [ 23 ] * xx [ 85 ] ; xx [ 23 ]
= xx [ 65 ] - ( xx [ 79 ] + xx [ 1 ] * xx [ 85 ] * xx [ 22 ] ) + xx [ 29 ] ;
xx [ 22 ] = xx [ 23 ] * xx [ 26 ] ; xx [ 65 ] = xx [ 26 ] * xx [ 86 ] ; xx [
70 ] = xx [ 82 ] * xx [ 60 ] - xx [ 57 ] * xx [ 81 ] ; xx [ 72 ] = xx [ 60 ]
* xx [ 81 ] + xx [ 82 ] * xx [ 57 ] ; xx [ 74 ] = xx [ 61 ] * xx [ 70 ] - xx
[ 72 ] * xx [ 62 ] ; xx [ 79 ] = xx [ 15 ] * xx [ 74 ] ; xx [ 85 ] = xx [ 1 ]
* ( xx [ 62 ] * xx [ 70 ] + xx [ 72 ] * xx [ 61 ] ) * xx [ 79 ] ; xx [ 87 ] =
xx [ 29 ] * xx [ 60 ] ; xx [ 88 ] = xx [ 66 ] + xx [ 67 ] + xx [ 1 ] * xx [
57 ] * xx [ 87 ] ; xx [ 57 ] = xx [ 71 ] - xx [ 1 ] * xx [ 87 ] * xx [ 60 ] +
xx [ 29 ] ; xx [ 60 ] = xx [ 82 ] * xx [ 57 ] ; xx [ 67 ] = xx [ 82 ] * xx [
88 ] ; xx [ 71 ] = xx [ 29 ] * xx [ 61 ] ; xx [ 87 ] = xx [ 29 ] - xx [ 1 ] *
xx [ 71 ] * xx [ 61 ] ; xx [ 61 ] = xx [ 72 ] * xx [ 87 ] ; xx [ 89 ] = xx [
1 ] * xx [ 62 ] * xx [ 71 ] ; xx [ 62 ] = xx [ 72 ] * xx [ 89 ] ; xx [ 71 ] =
xx [ 1 ] * xx [ 11 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 1 ] * xx [ 79 ] * xx [ 74 ]
; xx [ 90 ] = xx [ 0 ] ; xx [ 91 ] = xx [ 1 ] * xx [ 33 ] * ( xx [ 17 ] * xx
[ 18 ] - xx [ 6 ] * xx [ 32 ] ) + xx [ 47 ] + xx [ 1 ] * ( xx [ 40 ] - xx [
43 ] ) ; xx [ 92 ] = xx [ 46 ] + xx [ 1 ] * ( xx [ 38 ] * xx [ 18 ] - xx [ 32
] * xx [ 44 ] ) - xx [ 49 ] ; xx [ 93 ] = xx [ 51 ] - ( xx [ 12 ] * xx [ 45 ]
+ xx [ 48 ] * xx [ 50 ] ) * xx [ 1 ] - xx [ 49 ] ; xx [ 94 ] = xx [ 1 ] * (
xx [ 54 ] * xx [ 24 ] + xx [ 25 ] * xx [ 26 ] ) * xx [ 53 ] - ( xx [ 1 ] * (
xx [ 78 ] * xx [ 26 ] - xx [ 54 ] * xx [ 75 ] ) - xx [ 77 ] - xx [ 1 ] * xx [
54 ] * xx [ 80 ] ) ; xx [ 95 ] = xx [ 1 ] * xx [ 84 ] * ( xx [ 81 ] * xx [ 21
] - xx [ 82 ] * xx [ 56 ] ) - ( xx [ 86 ] - ( xx [ 54 ] * xx [ 22 ] + xx [ 65
] * xx [ 26 ] ) * xx [ 1 ] ) ; xx [ 96 ] = xx [ 85 ] + xx [ 88 ] + ( xx [ 60
] * xx [ 81 ] - xx [ 82 ] * xx [ 67 ] ) * xx [ 1 ] ; xx [ 97 ] = xx [ 85 ] -
( ( xx [ 61 ] * xx [ 70 ] + xx [ 72 ] * xx [ 62 ] ) * xx [ 1 ] - xx [ 89 ] )
; xx [ 98 ] = xx [ 58 ] ; xx [ 99 ] = xx [ 42 ] - ( xx [ 1 ] * xx [ 28 ] * xx
[ 33 ] + ( xx [ 43 ] + xx [ 40 ] ) * xx [ 1 ] ) + xx [ 15 ] ; xx [ 100 ] = xx
[ 39 ] - ( xx [ 71 ] + ( xx [ 44 ] * xx [ 18 ] + xx [ 32 ] * xx [ 38 ] ) * xx
[ 1 ] ) + xx [ 15 ] ; xx [ 101 ] = xx [ 13 ] + xx [ 1 ] * ( xx [ 50 ] * xx [
45 ] - xx [ 48 ] * xx [ 12 ] ) - xx [ 71 ] + xx [ 15 ] ; xx [ 102 ] = - ( xx
[ 73 ] - ( ( xx [ 54 ] * xx [ 78 ] + xx [ 75 ] * xx [ 26 ] ) * xx [ 1 ] + xx
[ 1 ] * xx [ 80 ] * xx [ 26 ] ) - xx [ 1 ] * xx [ 53 ] * xx [ 27 ] + xx [ 19
] ) ; xx [ 103 ] = - ( xx [ 23 ] + xx [ 1 ] * ( xx [ 54 ] * xx [ 65 ] - xx [
22 ] * xx [ 26 ] ) - xx [ 1 ] * xx [ 83 ] * xx [ 84 ] + xx [ 15 ] ) ; xx [
104 ] = - ( xx [ 57 ] - xx [ 1 ] * ( xx [ 67 ] * xx [ 81 ] + xx [ 82 ] * xx [
60 ] ) - xx [ 9 ] + xx [ 15 ] ) ; xx [ 105 ] = - ( xx [ 87 ] + xx [ 1 ] * (
xx [ 62 ] * xx [ 70 ] - xx [ 72 ] * xx [ 61 ] ) - xx [ 9 ] + xx [ 15 ] ) ; xx
[ 6 ] = xx [ 32 ] * xx [ 41 ] ; xx [ 11 ] = xx [ 35 ] * xx [ 32 ] ; xx [ 12 ]
= xx [ 30 ] * xx [ 15 ] ; xx [ 13 ] = xx [ 1 ] * xx [ 30 ] * xx [ 12 ] - xx [
7 ] ; xx [ 17 ] = xx [ 1 ] * xx [ 12 ] * xx [ 16 ] ; xx [ 12 ] = xx [ 8 ] *
xx [ 8 ] * xx [ 17 ] ; xx [ 16 ] = xx [ 8 ] * xx [ 13 ] * xx [ 8 ] ; xx [ 21
] = xx [ 48 ] * xx [ 36 ] ; xx [ 22 ] = xx [ 37 ] * xx [ 48 ] ; xx [ 23 ] =
xx [ 82 ] * xx [ 69 ] ; xx [ 24 ] = xx [ 82 ] * xx [ 63 ] ; xx [ 25 ] = xx [
76 ] * xx [ 26 ] ; xx [ 27 ] = xx [ 68 ] * xx [ 26 ] ; xx [ 28 ] = xx [ 15 ]
* xx [ 26 ] ; xx [ 30 ] = xx [ 72 ] * xx [ 64 ] ; xx [ 33 ] = xx [ 72 ] * xx
[ 66 ] ; xx [ 38 ] = - ( xx [ 37 ] + xx [ 35 ] + xx [ 1 ] * ( xx [ 6 ] * xx [
18 ] - xx [ 32 ] * xx [ 11 ] ) + xx [ 13 ] - xx [ 1 ] * ( xx [ 12 ] + xx [ 16
] ) + xx [ 4 ] - ( xx [ 21 ] * xx [ 45 ] + xx [ 48 ] * xx [ 22 ] ) * xx [ 1 ]
+ xx [ 71 ] - ( xx [ 64 ] + xx [ 63 ] - xx [ 1 ] * ( xx [ 23 ] * xx [ 81 ] +
xx [ 82 ] * xx [ 24 ] ) + xx [ 68 ] - xx [ 1 ] * ( xx [ 54 ] * xx [ 25 ] + xx
[ 27 ] * xx [ 26 ] ) - xx [ 1 ] * xx [ 28 ] * xx [ 26 ] - ( xx [ 72 ] * xx [
30 ] - xx [ 33 ] * xx [ 70 ] ) * xx [ 1 ] - xx [ 9 ] ) + xx [ 31 ] ) ; xx [
39 ] = - ( xx [ 36 ] + xx [ 1 ] * ( xx [ 22 ] * xx [ 45 ] - xx [ 48 ] * xx [
21 ] ) + xx [ 41 ] - ( xx [ 11 ] * xx [ 18 ] + xx [ 32 ] * xx [ 6 ] ) * xx [
1 ] - ( xx [ 17 ] + ( xx [ 16 ] - xx [ 12 ] ) * xx [ 1 ] ) - xx [ 49 ] - ( xx
[ 85 ] + ( xx [ 24 ] * xx [ 81 ] - xx [ 82 ] * xx [ 23 ] ) * xx [ 1 ] + ( xx
[ 54 ] * xx [ 27 ] - xx [ 25 ] * xx [ 26 ] ) * xx [ 1 ] + xx [ 76 ] + xx [ 1
] * xx [ 54 ] * xx [ 28 ] + xx [ 69 ] - xx [ 1 ] * ( xx [ 72 ] * xx [ 33 ] +
xx [ 30 ] * xx [ 70 ] ) + xx [ 66 ] ) - xx [ 34 ] ) ; memcpy ( xx + 60 , xx +
90 , 16 * sizeof ( double ) ) ; factorAndSolveWide ( 2 , 8 , xx + 60 , xx +
11 , xx + 16 , ii + 0 , xx + 38 , xx [ 3 ] , xx + 21 ) ; xx [ 6 ] = xx [ 4 ]
+ xx [ 21 ] ; xx [ 4 ] = xx [ 14 ] + xx [ 22 ] ; xx [ 9 ] = xx [ 5 ] + xx [
23 ] ; xx [ 5 ] = xx [ 10 ] + xx [ 24 ] ; xx [ 10 ] = xx [ 52 ] + xx [ 25 ] ;
xx [ 11 ] = xx [ 20 ] + xx [ 26 ] ; xx [ 12 ] = xx [ 55 ] + xx [ 27 ] ; xx [
13 ] = xx [ 59 ] + xx [ 28 ] ; xx [ 35 ] = xx [ 6 ] ; xx [ 36 ] = state [ 1 ]
; xx [ 37 ] = xx [ 4 ] ; xx [ 38 ] = state [ 3 ] ; xx [ 39 ] = xx [ 9 ] ; xx
[ 40 ] = state [ 5 ] ; xx [ 41 ] = xx [ 5 ] ; xx [ 42 ] = state [ 7 ] ; xx [
43 ] = xx [ 10 ] ; xx [ 44 ] = state [ 9 ] ; xx [ 45 ] = xx [ 11 ] ; xx [ 46
] = state [ 11 ] ; xx [ 47 ] = xx [ 12 ] ; xx [ 48 ] = state [ 13 ] ; xx [ 49
] = xx [ 13 ] ; xx [ 50 ] = state [ 15 ] ; xx [ 51 ] = state [ 16 ] ; xx [ 52
] = state [ 17 ] ; xx [ 14 ] = xx [ 5 ] * xx [ 2 ] ; xx [ 5 ] = sin ( xx [ 14
] ) ; xx [ 16 ] = xx [ 15 ] * xx [ 5 ] ; xx [ 17 ] = xx [ 1 ] * xx [ 16 ] *
xx [ 5 ] - xx [ 19 ] ; xx [ 18 ] = xx [ 9 ] * xx [ 2 ] ; xx [ 9 ] = sin ( xx
[ 18 ] ) ; xx [ 20 ] = xx [ 15 ] * xx [ 9 ] ; xx [ 21 ] = xx [ 1 ] * xx [ 20
] * xx [ 9 ] - xx [ 19 ] ; xx [ 22 ] = xx [ 4 ] * xx [ 2 ] ; xx [ 4 ] = xx [
8 ] * cos ( xx [ 22 ] ) ; xx [ 23 ] = xx [ 8 ] * sin ( xx [ 22 ] ) ; xx [ 22
] = xx [ 4 ] + xx [ 23 ] ; xx [ 24 ] = xx [ 8 ] * xx [ 22 ] ; xx [ 25 ] = xx
[ 4 ] - xx [ 23 ] ; xx [ 4 ] = xx [ 25 ] * xx [ 8 ] ; xx [ 23 ] = xx [ 24 ] -
xx [ 4 ] ; xx [ 26 ] = cos ( xx [ 18 ] ) ; xx [ 18 ] = xx [ 1 ] * xx [ 26 ] *
xx [ 20 ] ; xx [ 20 ] = xx [ 23 ] * xx [ 18 ] ; xx [ 27 ] = xx [ 24 ] + xx [
4 ] ; xx [ 4 ] = xx [ 21 ] * xx [ 23 ] ; xx [ 24 ] = xx [ 25 ] * xx [ 15 ] ;
xx [ 28 ] = xx [ 1 ] * xx [ 25 ] * xx [ 24 ] - xx [ 7 ] ; xx [ 7 ] = xx [ 1 ]
* xx [ 24 ] * xx [ 22 ] ; xx [ 22 ] = xx [ 8 ] * xx [ 8 ] * xx [ 7 ] ; xx [
24 ] = xx [ 8 ] * xx [ 28 ] * xx [ 8 ] ; xx [ 25 ] = xx [ 9 ] * xx [ 27 ] +
xx [ 23 ] * xx [ 26 ] ; xx [ 30 ] = cos ( xx [ 14 ] ) ; xx [ 14 ] = xx [ 1 ]
* xx [ 30 ] * xx [ 16 ] ; xx [ 16 ] = xx [ 25 ] * xx [ 14 ] ; xx [ 32 ] = xx
[ 23 ] * xx [ 9 ] - xx [ 26 ] * xx [ 27 ] ; xx [ 9 ] = xx [ 17 ] * xx [ 25 ]
; xx [ 26 ] = xx [ 25 ] * xx [ 30 ] - xx [ 5 ] * xx [ 32 ] ; xx [ 33 ] = xx [
15 ] * xx [ 26 ] ; xx [ 53 ] = xx [ 13 ] * xx [ 2 ] ; xx [ 13 ] = sin ( xx [
53 ] ) ; xx [ 54 ] = xx [ 15 ] * xx [ 13 ] ; xx [ 55 ] = xx [ 19 ] - xx [ 1 ]
* xx [ 54 ] * xx [ 13 ] ; xx [ 56 ] = xx [ 12 ] * xx [ 2 ] ; xx [ 12 ] = sin
( xx [ 56 ] ) ; xx [ 57 ] = xx [ 15 ] * xx [ 12 ] ; xx [ 59 ] = xx [ 19 ] -
xx [ 1 ] * xx [ 57 ] * xx [ 12 ] ; xx [ 60 ] = xx [ 10 ] * xx [ 2 ] ; xx [ 10
] = cos ( xx [ 60 ] ) ; xx [ 61 ] = xx [ 11 ] * xx [ 2 ] ; xx [ 11 ] = sin (
xx [ 61 ] ) ; xx [ 62 ] = cos ( xx [ 61 ] ) ; xx [ 61 ] = sin ( xx [ 60 ] ) ;
xx [ 60 ] = xx [ 10 ] * xx [ 11 ] + xx [ 62 ] * xx [ 61 ] ; xx [ 63 ] = cos (
xx [ 56 ] ) ; xx [ 56 ] = xx [ 1 ] * xx [ 63 ] * xx [ 57 ] ; xx [ 57 ] = xx [
60 ] * xx [ 56 ] ; xx [ 64 ] = xx [ 10 ] * xx [ 62 ] - xx [ 61 ] * xx [ 11 ]
; xx [ 65 ] = xx [ 60 ] * xx [ 59 ] ; xx [ 66 ] = xx [ 15 ] * xx [ 11 ] ; xx
[ 67 ] = xx [ 19 ] - xx [ 1 ] * xx [ 66 ] * xx [ 11 ] ; xx [ 11 ] = xx [ 1 ]
* xx [ 62 ] * xx [ 66 ] ; xx [ 62 ] = xx [ 11 ] * xx [ 61 ] ; xx [ 66 ] = xx
[ 67 ] * xx [ 61 ] ; xx [ 68 ] = xx [ 15 ] * xx [ 61 ] ; xx [ 69 ] = xx [ 12
] * xx [ 64 ] + xx [ 60 ] * xx [ 63 ] ; xx [ 70 ] = xx [ 69 ] * xx [ 55 ] ;
xx [ 71 ] = cos ( xx [ 53 ] ) ; xx [ 53 ] = xx [ 1 ] * xx [ 71 ] * xx [ 54 ]
; xx [ 54 ] = xx [ 69 ] * xx [ 53 ] ; xx [ 72 ] = xx [ 60 ] * xx [ 12 ] - xx
[ 63 ] * xx [ 64 ] ; xx [ 12 ] = xx [ 13 ] * xx [ 72 ] - xx [ 69 ] * xx [ 71
] ; xx [ 63 ] = xx [ 15 ] * xx [ 12 ] ; xx [ 73 ] = fabs ( xx [ 17 ] + xx [
21 ] + xx [ 1 ] * ( xx [ 20 ] * xx [ 27 ] - xx [ 23 ] * xx [ 4 ] ) + xx [ 28
] - xx [ 1 ] * ( xx [ 22 ] + xx [ 24 ] ) + xx [ 6 ] - ( xx [ 16 ] * xx [ 32 ]
+ xx [ 25 ] * xx [ 9 ] ) * xx [ 1 ] + xx [ 1 ] * xx [ 33 ] * xx [ 26 ] - ( xx
[ 55 ] + xx [ 59 ] - xx [ 1 ] * ( xx [ 57 ] * xx [ 64 ] + xx [ 60 ] * xx [ 65
] ) + xx [ 67 ] - xx [ 1 ] * ( xx [ 10 ] * xx [ 62 ] + xx [ 66 ] * xx [ 61 ]
) - xx [ 1 ] * xx [ 68 ] * xx [ 61 ] - ( xx [ 69 ] * xx [ 70 ] - xx [ 54 ] *
xx [ 72 ] ) * xx [ 1 ] - xx [ 1 ] * xx [ 63 ] * xx [ 12 ] ) + xx [ 31 ] ) ;
xx [ 74 ] = fabs ( xx [ 14 ] + xx [ 1 ] * ( xx [ 9 ] * xx [ 32 ] - xx [ 25 ]
* xx [ 16 ] ) + xx [ 18 ] - ( xx [ 4 ] * xx [ 27 ] + xx [ 23 ] * xx [ 20 ] )
* xx [ 1 ] - ( xx [ 7 ] + ( xx [ 24 ] - xx [ 22 ] ) * xx [ 1 ] ) - xx [ 1 ] *
( xx [ 30 ] * xx [ 32 ] + xx [ 25 ] * xx [ 5 ] ) * xx [ 33 ] - ( xx [ 1 ] * (
xx [ 71 ] * xx [ 72 ] + xx [ 69 ] * xx [ 13 ] ) * xx [ 63 ] + ( xx [ 65 ] *
xx [ 64 ] - xx [ 60 ] * xx [ 57 ] ) * xx [ 1 ] + ( xx [ 10 ] * xx [ 66 ] - xx
[ 62 ] * xx [ 61 ] ) * xx [ 1 ] + xx [ 11 ] + xx [ 1 ] * xx [ 10 ] * xx [ 68
] + xx [ 56 ] - xx [ 1 ] * ( xx [ 69 ] * xx [ 54 ] + xx [ 70 ] * xx [ 72 ] )
+ xx [ 53 ] ) - xx [ 34 ] ) ; ii [ 0 ] = 73 ; { int ll ; for ( ll = 74 ; ll <
75 ; ++ ll ) if ( xx [ ll ] > xx [ ii [ 0 ] ] ) ii [ 0 ] = ll ; } ii [ 0 ] -=
73 ; xx [ 4 ] = xx [ 73 + ( ii [ 0 ] ) ] ; xx [ 5 ] = 1.0e-9 ; if ( xx [ 4 ]
> xx [ 5 ] ) { switch ( ii [ 0 ] ) { case 0 : case 1 : { return
sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:constraintViolation" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint5' kinematic constraints cannot be maintained. Check solver type and consistency tolerance in the Simscape Solver Configuration block. Check Simulink solver type and tolerances in Model Configuration Parameters. A kinematic singularity might be the source of this problem."
, neDiagMgr ) ; } } } xx [ 4 ] = xx [ 2 ] * xx [ 39 ] ; xx [ 6 ] = cos ( xx [
4 ] ) ; xx [ 7 ] = xx [ 2 ] * xx [ 41 ] ; xx [ 9 ] = sin ( xx [ 7 ] ) ; xx [
10 ] = cos ( xx [ 7 ] ) ; xx [ 7 ] = sin ( xx [ 4 ] ) ; xx [ 4 ] = xx [ 6 ] *
xx [ 9 ] + xx [ 10 ] * xx [ 7 ] ; xx [ 11 ] = xx [ 2 ] * xx [ 37 ] ; xx [ 12
] = xx [ 8 ] * cos ( xx [ 11 ] ) ; xx [ 13 ] = xx [ 8 ] * sin ( xx [ 11 ] ) ;
xx [ 11 ] = xx [ 12 ] + xx [ 13 ] ; xx [ 14 ] = xx [ 8 ] * xx [ 11 ] ; xx [
16 ] = xx [ 12 ] - xx [ 13 ] ; xx [ 12 ] = xx [ 16 ] * xx [ 8 ] ; xx [ 13 ] =
xx [ 14 ] + xx [ 12 ] ; xx [ 17 ] = xx [ 14 ] - xx [ 12 ] ; xx [ 12 ] = xx [
6 ] * xx [ 10 ] - xx [ 7 ] * xx [ 9 ] ; xx [ 14 ] = xx [ 4 ] * xx [ 13 ] + xx
[ 17 ] * xx [ 12 ] ; xx [ 18 ] = xx [ 14 ] * xx [ 15 ] ; xx [ 20 ] = xx [ 15
] * xx [ 9 ] ; xx [ 21 ] = xx [ 1 ] * xx [ 10 ] * xx [ 20 ] ; xx [ 22 ] = xx
[ 1 ] * xx [ 20 ] * xx [ 9 ] - xx [ 19 ] ; xx [ 20 ] = xx [ 22 ] * xx [ 7 ] ;
xx [ 23 ] = xx [ 21 ] * xx [ 7 ] ; xx [ 24 ] = xx [ 21 ] - ( xx [ 6 ] * xx [
20 ] + xx [ 23 ] * xx [ 7 ] ) * xx [ 1 ] ; xx [ 25 ] = xx [ 15 ] * xx [ 7 ] ;
xx [ 26 ] = xx [ 1 ] * xx [ 6 ] * xx [ 25 ] ; xx [ 27 ] = xx [ 24 ] + xx [ 26
] ; xx [ 28 ] = xx [ 16 ] * xx [ 27 ] ; xx [ 30 ] = xx [ 6 ] * xx [ 23 ] ; xx
[ 23 ] = xx [ 20 ] * xx [ 7 ] ; xx [ 20 ] = xx [ 1 ] * xx [ 25 ] * xx [ 7 ] -
xx [ 19 ] ; xx [ 25 ] = xx [ 22 ] + xx [ 1 ] * ( xx [ 30 ] - xx [ 23 ] ) + xx
[ 20 ] ; xx [ 31 ] = xx [ 25 ] * xx [ 16 ] ; xx [ 32 ] = xx [ 16 ] * xx [ 29
] ; xx [ 33 ] = xx [ 32 ] * xx [ 11 ] ; xx [ 34 ] = xx [ 1 ] * xx [ 33 ] ; xx
[ 53 ] = xx [ 27 ] - ( xx [ 16 ] * xx [ 28 ] - xx [ 31 ] * xx [ 11 ] ) * xx [
1 ] - xx [ 34 ] ; xx [ 27 ] = xx [ 1 ] * xx [ 16 ] * xx [ 32 ] ; xx [ 32 ] =
xx [ 1 ] * ( xx [ 16 ] * xx [ 31 ] + xx [ 28 ] * xx [ 11 ] ) - ( xx [ 25 ] +
xx [ 27 ] ) + xx [ 29 ] ; xx [ 11 ] = xx [ 8 ] * xx [ 32 ] * xx [ 8 ] ; xx [
16 ] = xx [ 8 ] * xx [ 53 ] * xx [ 8 ] ; xx [ 25 ] = xx [ 29 ] * xx [ 7 ] ;
xx [ 28 ] = xx [ 6 ] * xx [ 25 ] ; xx [ 31 ] = xx [ 1 ] * xx [ 28 ] ; xx [ 54
] = xx [ 24 ] + xx [ 31 ] ; xx [ 24 ] = xx [ 1 ] * xx [ 25 ] * xx [ 7 ] ; xx
[ 25 ] = xx [ 1 ] * ( xx [ 23 ] - xx [ 30 ] ) - ( xx [ 22 ] + xx [ 24 ] ) +
xx [ 29 ] ; xx [ 23 ] = xx [ 17 ] * xx [ 25 ] ; xx [ 30 ] = xx [ 17 ] * xx [
54 ] ; xx [ 55 ] = xx [ 17 ] * xx [ 7 ] - xx [ 6 ] * xx [ 13 ] ; xx [ 56 ] =
xx [ 7 ] * xx [ 13 ] + xx [ 17 ] * xx [ 6 ] ; xx [ 6 ] = xx [ 10 ] * xx [ 55
] + xx [ 56 ] * xx [ 9 ] ; xx [ 7 ] = xx [ 56 ] * xx [ 10 ] - xx [ 9 ] * xx [
55 ] ; xx [ 57 ] = xx [ 15 ] * xx [ 7 ] ; xx [ 59 ] = xx [ 1 ] * xx [ 6 ] *
xx [ 57 ] ; xx [ 60 ] = xx [ 29 ] * xx [ 9 ] ; xx [ 61 ] = xx [ 10 ] * xx [
60 ] ; xx [ 10 ] = xx [ 1 ] * xx [ 61 ] ; xx [ 62 ] = xx [ 29 ] - xx [ 1 ] *
xx [ 60 ] * xx [ 9 ] ; xx [ 9 ] = xx [ 56 ] * xx [ 62 ] ; xx [ 60 ] = xx [ 56
] * xx [ 10 ] ; xx [ 63 ] = xx [ 2 ] * xx [ 43 ] ; xx [ 64 ] = cos ( xx [ 63
] ) ; xx [ 65 ] = xx [ 2 ] * xx [ 47 ] ; xx [ 66 ] = cos ( xx [ 65 ] ) ; xx [
67 ] = xx [ 2 ] * xx [ 49 ] ; xx [ 68 ] = sin ( xx [ 67 ] ) ; xx [ 69 ] = cos
( xx [ 67 ] ) ; xx [ 67 ] = sin ( xx [ 65 ] ) ; xx [ 65 ] = xx [ 66 ] * xx [
68 ] + xx [ 69 ] * xx [ 67 ] ; xx [ 70 ] = xx [ 2 ] * xx [ 45 ] ; xx [ 71 ] =
sin ( xx [ 70 ] ) ; xx [ 72 ] = cos ( xx [ 70 ] ) ; xx [ 70 ] = xx [ 66 ] *
xx [ 69 ] - xx [ 67 ] * xx [ 68 ] ; xx [ 73 ] = xx [ 65 ] * xx [ 71 ] - xx [
72 ] * xx [ 70 ] ; xx [ 74 ] = xx [ 65 ] * xx [ 72 ] + xx [ 71 ] * xx [ 70 ]
; xx [ 75 ] = sin ( xx [ 63 ] ) ; xx [ 63 ] = xx [ 75 ] * xx [ 73 ] - xx [ 74
] * xx [ 64 ] ; xx [ 76 ] = xx [ 15 ] * xx [ 63 ] ; xx [ 77 ] = xx [ 15 ] *
xx [ 68 ] ; xx [ 78 ] = xx [ 19 ] - xx [ 1 ] * xx [ 77 ] * xx [ 68 ] ; xx [
79 ] = xx [ 78 ] * xx [ 67 ] ; xx [ 80 ] = xx [ 1 ] * xx [ 69 ] * xx [ 77 ] ;
xx [ 77 ] = xx [ 80 ] * xx [ 67 ] ; xx [ 81 ] = ( xx [ 66 ] * xx [ 79 ] - xx
[ 77 ] * xx [ 67 ] ) * xx [ 1 ] ; xx [ 82 ] = xx [ 15 ] * xx [ 67 ] ; xx [ 83
] = xx [ 1 ] * xx [ 66 ] * xx [ 82 ] ; xx [ 84 ] = xx [ 81 ] + xx [ 80 ] + xx
[ 83 ] ; xx [ 85 ] = xx [ 78 ] - xx [ 1 ] * ( xx [ 66 ] * xx [ 77 ] + xx [ 79
] * xx [ 67 ] ) ; xx [ 77 ] = xx [ 19 ] - xx [ 1 ] * xx [ 82 ] * xx [ 67 ] ;
xx [ 79 ] = xx [ 85 ] + xx [ 77 ] ; xx [ 82 ] = xx [ 79 ] * xx [ 71 ] ; xx [
86 ] = xx [ 72 ] * xx [ 82 ] ; xx [ 87 ] = xx [ 71 ] * xx [ 84 ] ; xx [ 88 ]
= xx [ 87 ] * xx [ 71 ] ; xx [ 89 ] = xx [ 15 ] * xx [ 71 ] ; xx [ 90 ] = xx
[ 1 ] * xx [ 72 ] * xx [ 89 ] ; xx [ 91 ] = xx [ 84 ] + xx [ 1 ] * ( xx [ 86
] - xx [ 88 ] ) + xx [ 90 ] ; xx [ 92 ] = xx [ 75 ] * xx [ 91 ] ; xx [ 93 ] =
( xx [ 72 ] * xx [ 87 ] + xx [ 82 ] * xx [ 71 ] ) * xx [ 1 ] ; xx [ 82 ] = xx
[ 19 ] - xx [ 1 ] * xx [ 89 ] * xx [ 71 ] ; xx [ 87 ] = xx [ 79 ] - xx [ 93 ]
+ xx [ 82 ] ; xx [ 89 ] = xx [ 87 ] * xx [ 75 ] ; xx [ 94 ] = xx [ 29 ] * xx
[ 75 ] ; xx [ 95 ] = xx [ 64 ] * xx [ 94 ] ; xx [ 96 ] = xx [ 1 ] * xx [ 95 ]
; xx [ 97 ] = xx [ 64 ] * xx [ 72 ] - xx [ 75 ] * xx [ 71 ] ; xx [ 98 ] = xx
[ 64 ] * xx [ 71 ] + xx [ 72 ] * xx [ 75 ] ; xx [ 99 ] = xx [ 65 ] * xx [ 97
] + xx [ 98 ] * xx [ 70 ] ; xx [ 100 ] = xx [ 99 ] * xx [ 15 ] ; xx [ 101 ] =
xx [ 29 ] * xx [ 71 ] ; xx [ 102 ] = xx [ 72 ] * xx [ 101 ] ; xx [ 72 ] = xx
[ 1 ] * xx [ 102 ] ; xx [ 103 ] = xx [ 1 ] * ( xx [ 88 ] - xx [ 86 ] ) - xx [
84 ] - xx [ 72 ] ; xx [ 84 ] = xx [ 1 ] * xx [ 101 ] * xx [ 71 ] ; xx [ 71 ]
= xx [ 79 ] - ( xx [ 93 ] + xx [ 84 ] ) + xx [ 29 ] ; xx [ 79 ] = xx [ 71 ] *
xx [ 75 ] ; xx [ 86 ] = xx [ 75 ] * xx [ 103 ] ; xx [ 88 ] = xx [ 98 ] * xx [
67 ] - xx [ 66 ] * xx [ 97 ] ; xx [ 93 ] = xx [ 67 ] * xx [ 97 ] + xx [ 98 ]
* xx [ 66 ] ; xx [ 101 ] = xx [ 69 ] * xx [ 88 ] + xx [ 93 ] * xx [ 68 ] ; xx
[ 104 ] = xx [ 68 ] * xx [ 88 ] - xx [ 93 ] * xx [ 69 ] ; xx [ 105 ] = xx [
15 ] * xx [ 104 ] ; xx [ 106 ] = xx [ 1 ] * xx [ 101 ] * xx [ 105 ] ; xx [
107 ] = xx [ 29 ] * xx [ 67 ] ; xx [ 108 ] = xx [ 66 ] * xx [ 107 ] ; xx [ 66
] = xx [ 1 ] * xx [ 108 ] ; xx [ 109 ] = xx [ 80 ] + xx [ 81 ] + xx [ 66 ] ;
xx [ 81 ] = xx [ 1 ] * xx [ 107 ] * xx [ 67 ] ; xx [ 67 ] = xx [ 85 ] - xx [
81 ] + xx [ 29 ] ; xx [ 85 ] = xx [ 98 ] * xx [ 67 ] ; xx [ 107 ] = xx [ 98 ]
* xx [ 109 ] ; xx [ 110 ] = xx [ 29 ] * xx [ 68 ] ; xx [ 111 ] = xx [ 29 ] -
xx [ 1 ] * xx [ 110 ] * xx [ 68 ] ; xx [ 68 ] = xx [ 93 ] * xx [ 111 ] ; xx [
112 ] = xx [ 69 ] * xx [ 110 ] ; xx [ 69 ] = xx [ 1 ] * xx [ 112 ] ; xx [ 110
] = xx [ 93 ] * xx [ 69 ] ; xx [ 113 ] = xx [ 1 ] * xx [ 57 ] * xx [ 7 ] ; xx
[ 57 ] = xx [ 1 ] * xx [ 94 ] * xx [ 75 ] ; xx [ 94 ] = xx [ 1 ] * xx [ 105 ]
* xx [ 104 ] ; xx [ 114 ] = xx [ 0 ] ; xx [ 115 ] = xx [ 1 ] * xx [ 18 ] * (
xx [ 12 ] * xx [ 13 ] - xx [ 4 ] * xx [ 17 ] ) + xx [ 53 ] + xx [ 1 ] * ( xx
[ 11 ] - xx [ 16 ] ) ; xx [ 116 ] = xx [ 54 ] + xx [ 1 ] * ( xx [ 23 ] * xx [
13 ] - xx [ 17 ] * xx [ 30 ] ) - xx [ 59 ] ; xx [ 117 ] = xx [ 10 ] - ( xx [
9 ] * xx [ 55 ] + xx [ 56 ] * xx [ 60 ] ) * xx [ 1 ] - xx [ 59 ] ; xx [ 118 ]
= xx [ 1 ] * ( xx [ 64 ] * xx [ 73 ] + xx [ 74 ] * xx [ 75 ] ) * xx [ 76 ] -
( xx [ 1 ] * ( xx [ 92 ] * xx [ 75 ] - xx [ 64 ] * xx [ 89 ] ) - xx [ 91 ] -
xx [ 96 ] ) ; xx [ 119 ] = xx [ 1 ] * xx [ 100 ] * ( xx [ 97 ] * xx [ 70 ] -
xx [ 98 ] * xx [ 65 ] ) - ( xx [ 103 ] - ( xx [ 64 ] * xx [ 79 ] + xx [ 86 ]
* xx [ 75 ] ) * xx [ 1 ] ) ; xx [ 120 ] = xx [ 106 ] + xx [ 109 ] + ( xx [ 85
] * xx [ 97 ] - xx [ 98 ] * xx [ 107 ] ) * xx [ 1 ] ; xx [ 121 ] = xx [ 106 ]
- ( ( xx [ 68 ] * xx [ 88 ] + xx [ 93 ] * xx [ 110 ] ) * xx [ 1 ] - xx [ 69 ]
) ; xx [ 122 ] = xx [ 58 ] ; xx [ 123 ] = xx [ 32 ] - ( xx [ 1 ] * xx [ 14 ]
* xx [ 18 ] + ( xx [ 16 ] + xx [ 11 ] ) * xx [ 1 ] ) + xx [ 15 ] ; xx [ 124 ]
= xx [ 25 ] - ( xx [ 113 ] + ( xx [ 30 ] * xx [ 13 ] + xx [ 17 ] * xx [ 23 ]
) * xx [ 1 ] ) + xx [ 15 ] ; xx [ 125 ] = xx [ 62 ] + xx [ 1 ] * ( xx [ 60 ]
* xx [ 55 ] - xx [ 56 ] * xx [ 9 ] ) - xx [ 113 ] + xx [ 15 ] ; xx [ 126 ] =
- ( xx [ 87 ] - ( ( xx [ 64 ] * xx [ 92 ] + xx [ 89 ] * xx [ 75 ] ) * xx [ 1
] + xx [ 57 ] ) - xx [ 1 ] * xx [ 76 ] * xx [ 63 ] + xx [ 19 ] ) ; xx [ 127 ]
= - ( xx [ 71 ] + xx [ 1 ] * ( xx [ 64 ] * xx [ 86 ] - xx [ 79 ] * xx [ 75 ]
) - xx [ 1 ] * xx [ 99 ] * xx [ 100 ] + xx [ 15 ] ) ; xx [ 128 ] = - ( xx [
67 ] - xx [ 1 ] * ( xx [ 107 ] * xx [ 97 ] + xx [ 98 ] * xx [ 85 ] ) - xx [
94 ] + xx [ 15 ] ) ; xx [ 129 ] = - ( xx [ 111 ] + xx [ 1 ] * ( xx [ 110 ] *
xx [ 88 ] - xx [ 93 ] * xx [ 68 ] ) - xx [ 94 ] + xx [ 15 ] ) ; xx [ 0 ] = xx
[ 44 ] + xx [ 46 ] ; xx [ 4 ] = xx [ 0 ] + xx [ 48 ] ; xx [ 9 ] = xx [ 50 ] *
xx [ 111 ] + xx [ 4 ] * xx [ 78 ] ; xx [ 11 ] = xx [ 9 ] * xx [ 93 ] ; xx [
12 ] = xx [ 4 ] * xx [ 80 ] + xx [ 1 ] * xx [ 50 ] * xx [ 112 ] ; xx [ 14 ] =
xx [ 93 ] * xx [ 12 ] ; xx [ 16 ] = xx [ 44 ] * xx [ 90 ] + xx [ 1 ] * xx [
46 ] * xx [ 102 ] ; xx [ 18 ] = xx [ 29 ] - xx [ 84 ] ; xx [ 19 ] = xx [ 46 ]
* xx [ 18 ] + xx [ 82 ] * xx [ 44 ] ; xx [ 23 ] = xx [ 19 ] * xx [ 75 ] ; xx
[ 25 ] = xx [ 75 ] * xx [ 16 ] ; xx [ 30 ] = xx [ 0 ] * xx [ 83 ] + xx [ 1 ]
* xx [ 48 ] * xx [ 108 ] ; xx [ 32 ] = xx [ 29 ] - xx [ 81 ] ; xx [ 53 ] = xx
[ 48 ] * xx [ 32 ] + xx [ 0 ] * xx [ 77 ] ; xx [ 0 ] = xx [ 98 ] * xx [ 53 ]
; xx [ 54 ] = xx [ 98 ] * xx [ 30 ] ; xx [ 58 ] = ( xx [ 4 ] + xx [ 50 ] ) *
xx [ 15 ] ; xx [ 4 ] = xx [ 58 ] * xx [ 104 ] ; xx [ 59 ] = xx [ 29 ] - xx [
27 ] ; xx [ 27 ] = xx [ 38 ] * xx [ 59 ] ; xx [ 60 ] = xx [ 8 ] * xx [ 8 ] *
xx [ 27 ] ; xx [ 63 ] = xx [ 1 ] * xx [ 38 ] * xx [ 33 ] ; xx [ 33 ] = xx [ 8
] * xx [ 8 ] * xx [ 63 ] ; xx [ 65 ] = xx [ 1 ] * xx [ 40 ] * xx [ 28 ] + xx
[ 38 ] * xx [ 26 ] ; xx [ 28 ] = xx [ 29 ] - xx [ 24 ] ; xx [ 24 ] = xx [ 40
] * xx [ 28 ] - xx [ 20 ] * xx [ 38 ] ; xx [ 67 ] = xx [ 17 ] * xx [ 24 ] ;
xx [ 68 ] = xx [ 65 ] * xx [ 17 ] ; xx [ 70 ] = xx [ 38 ] + xx [ 40 ] ; xx [
71 ] = xx [ 1 ] * xx [ 42 ] * xx [ 61 ] + xx [ 70 ] * xx [ 21 ] ; xx [ 61 ] =
xx [ 42 ] * xx [ 62 ] - xx [ 70 ] * xx [ 22 ] ; xx [ 73 ] = xx [ 56 ] * xx [
61 ] ; xx [ 74 ] = xx [ 71 ] * xx [ 56 ] ; xx [ 76 ] = ( xx [ 70 ] + xx [ 42
] ) * xx [ 15 ] ; xx [ 70 ] = xx [ 76 ] * xx [ 7 ] ; xx [ 79 ] = xx [ 29 ] -
xx [ 57 ] ; xx [ 84 ] = xx [ 1 ] * ( xx [ 11 ] * xx [ 88 ] + xx [ 93 ] * xx [
14 ] ) - xx [ 12 ] - ( xx [ 16 ] + ( xx [ 64 ] * xx [ 23 ] - xx [ 25 ] * xx [
75 ] ) * xx [ 1 ] + xx [ 1 ] * xx [ 44 ] * xx [ 95 ] + xx [ 30 ] + ( xx [ 0 ]
* xx [ 97 ] - xx [ 98 ] * xx [ 54 ] ) * xx [ 1 ] ) - xx [ 1 ] * xx [ 101 ] *
xx [ 4 ] - ( xx [ 36 ] + xx [ 1 ] * ( xx [ 60 ] + xx [ 33 ] ) - xx [ 63 ] +
xx [ 65 ] + xx [ 1 ] * ( xx [ 67 ] * xx [ 13 ] - xx [ 17 ] * xx [ 68 ] ) + xx
[ 71 ] - ( xx [ 73 ] * xx [ 55 ] + xx [ 56 ] * xx [ 74 ] ) * xx [ 1 ] - xx [
1 ] * xx [ 6 ] * xx [ 70 ] ) ; xx [ 85 ] = xx [ 58 ] - xx [ 1 ] * xx [ 4 ] *
xx [ 104 ] + xx [ 44 ] * xx [ 79 ] + xx [ 19 ] - xx [ 1 ] * ( xx [ 64 ] * xx
[ 25 ] + xx [ 23 ] * xx [ 75 ] ) + xx [ 53 ] - xx [ 1 ] * ( xx [ 54 ] * xx [
97 ] + xx [ 98 ] * xx [ 0 ] ) + xx [ 9 ] - ( xx [ 93 ] * xx [ 11 ] - xx [ 14
] * xx [ 88 ] ) * xx [ 1 ] - ( xx [ 76 ] - xx [ 1 ] * xx [ 70 ] * xx [ 7 ] +
xx [ 27 ] - ( xx [ 60 ] - xx [ 33 ] ) * xx [ 1 ] + xx [ 24 ] - ( xx [ 68 ] *
xx [ 13 ] + xx [ 17 ] * xx [ 67 ] ) * xx [ 1 ] + xx [ 61 ] + xx [ 1 ] * ( xx
[ 74 ] * xx [ 55 ] - xx [ 56 ] * xx [ 73 ] ) ) ; memcpy ( xx + 138 , xx + 114
, 16 * sizeof ( double ) ) ; factorAndSolveWide ( 2 , 8 , xx + 138 , xx + 11
, xx + 23 , ii + 0 , xx + 84 , xx [ 3 ] , xx + 130 ) ; xx [ 0 ] = xx [ 36 ] +
xx [ 130 ] ; xx [ 3 ] = xx [ 38 ] + xx [ 131 ] ; xx [ 4 ] = xx [ 40 ] + xx [
132 ] ; xx [ 9 ] = xx [ 42 ] + xx [ 133 ] ; xx [ 11 ] = xx [ 44 ] + xx [ 134
] ; xx [ 12 ] = xx [ 46 ] + xx [ 135 ] ; xx [ 14 ] = xx [ 48 ] + xx [ 136 ] ;
xx [ 16 ] = xx [ 50 ] + xx [ 137 ] ; xx [ 112 ] = xx [ 35 ] ; xx [ 113 ] = xx
[ 0 ] ; xx [ 114 ] = xx [ 37 ] ; xx [ 115 ] = xx [ 3 ] ; xx [ 116 ] = xx [ 39
] ; xx [ 117 ] = xx [ 4 ] ; xx [ 118 ] = xx [ 41 ] ; xx [ 119 ] = xx [ 9 ] ;
xx [ 120 ] = xx [ 43 ] ; xx [ 121 ] = xx [ 11 ] ; xx [ 122 ] = xx [ 45 ] ; xx
[ 123 ] = xx [ 12 ] ; xx [ 124 ] = xx [ 47 ] ; xx [ 125 ] = xx [ 14 ] ; xx [
126 ] = xx [ 49 ] ; xx [ 127 ] = xx [ 16 ] ; xx [ 128 ] = xx [ 51 ] ; xx [
129 ] = xx [ 52 ] ; xx [ 19 ] = xx [ 3 ] * xx [ 59 ] ; xx [ 23 ] = xx [ 8 ] *
xx [ 8 ] * xx [ 19 ] ; xx [ 24 ] = xx [ 3 ] * xx [ 34 ] ; xx [ 25 ] = xx [ 8
] * xx [ 8 ] * xx [ 24 ] ; xx [ 27 ] = xx [ 4 ] * xx [ 31 ] + xx [ 3 ] * xx [
26 ] ; xx [ 26 ] = xx [ 4 ] * xx [ 28 ] - xx [ 3 ] * xx [ 20 ] ; xx [ 20 ] =
xx [ 17 ] * xx [ 26 ] ; xx [ 28 ] = xx [ 27 ] * xx [ 17 ] ; xx [ 29 ] = xx [
3 ] + xx [ 4 ] ; xx [ 3 ] = xx [ 9 ] * xx [ 10 ] + xx [ 29 ] * xx [ 21 ] ; xx
[ 4 ] = xx [ 9 ] * xx [ 62 ] - xx [ 29 ] * xx [ 22 ] ; xx [ 10 ] = xx [ 56 ]
* xx [ 4 ] ; xx [ 21 ] = xx [ 3 ] * xx [ 56 ] ; xx [ 22 ] = ( xx [ 29 ] + xx
[ 9 ] ) * xx [ 15 ] ; xx [ 9 ] = xx [ 22 ] * xx [ 7 ] ; xx [ 29 ] = xx [ 11 ]
+ xx [ 12 ] ; xx [ 30 ] = xx [ 29 ] + xx [ 14 ] ; xx [ 31 ] = xx [ 16 ] * xx
[ 111 ] + xx [ 30 ] * xx [ 78 ] ; xx [ 33 ] = xx [ 31 ] * xx [ 93 ] ; xx [ 34
] = xx [ 30 ] * xx [ 80 ] + xx [ 16 ] * xx [ 69 ] ; xx [ 35 ] = xx [ 93 ] *
xx [ 34 ] ; xx [ 36 ] = xx [ 11 ] * xx [ 90 ] + xx [ 12 ] * xx [ 72 ] ; xx [
37 ] = xx [ 12 ] * xx [ 18 ] + xx [ 11 ] * xx [ 82 ] ; xx [ 12 ] = xx [ 37 ]
* xx [ 75 ] ; xx [ 18 ] = xx [ 75 ] * xx [ 36 ] ; xx [ 38 ] = xx [ 29 ] * xx
[ 83 ] + xx [ 14 ] * xx [ 66 ] ; xx [ 39 ] = xx [ 14 ] * xx [ 32 ] + xx [ 29
] * xx [ 77 ] ; xx [ 14 ] = xx [ 39 ] * xx [ 98 ] ; xx [ 29 ] = xx [ 98 ] *
xx [ 38 ] ; xx [ 32 ] = ( xx [ 30 ] + xx [ 16 ] ) * xx [ 15 ] ; xx [ 15 ] =
xx [ 32 ] * xx [ 104 ] ; xx [ 40 ] = fabs ( xx [ 0 ] + xx [ 1 ] * ( xx [ 23 ]
+ xx [ 25 ] ) - xx [ 24 ] + xx [ 27 ] + xx [ 1 ] * ( xx [ 20 ] * xx [ 13 ] -
xx [ 17 ] * xx [ 28 ] ) + xx [ 3 ] - ( xx [ 10 ] * xx [ 55 ] + xx [ 56 ] * xx
[ 21 ] ) * xx [ 1 ] - xx [ 1 ] * xx [ 6 ] * xx [ 9 ] - ( xx [ 1 ] * ( xx [ 33
] * xx [ 88 ] + xx [ 93 ] * xx [ 35 ] ) - xx [ 34 ] - ( xx [ 36 ] + ( xx [ 64
] * xx [ 12 ] - xx [ 18 ] * xx [ 75 ] ) * xx [ 1 ] + xx [ 11 ] * xx [ 96 ] +
xx [ 38 ] + ( xx [ 14 ] * xx [ 97 ] - xx [ 98 ] * xx [ 29 ] ) * xx [ 1 ] ) -
xx [ 1 ] * xx [ 101 ] * xx [ 15 ] ) ) ; xx [ 41 ] = fabs ( xx [ 22 ] - xx [ 1
] * xx [ 9 ] * xx [ 7 ] + xx [ 19 ] - ( xx [ 23 ] - xx [ 25 ] ) * xx [ 1 ] +
xx [ 26 ] - ( xx [ 28 ] * xx [ 13 ] + xx [ 17 ] * xx [ 20 ] ) * xx [ 1 ] + xx
[ 4 ] + xx [ 1 ] * ( xx [ 21 ] * xx [ 55 ] - xx [ 56 ] * xx [ 10 ] ) - ( xx [
32 ] - xx [ 1 ] * xx [ 15 ] * xx [ 104 ] + xx [ 11 ] * xx [ 79 ] + xx [ 37 ]
- xx [ 1 ] * ( xx [ 64 ] * xx [ 18 ] + xx [ 12 ] * xx [ 75 ] ) + xx [ 39 ] -
xx [ 1 ] * ( xx [ 29 ] * xx [ 97 ] + xx [ 98 ] * xx [ 14 ] ) + xx [ 31 ] - (
xx [ 93 ] * xx [ 33 ] - xx [ 35 ] * xx [ 88 ] ) * xx [ 1 ] ) ) ; ii [ 0 ] =
40 ; { int ll ; for ( ll = 41 ; ll < 42 ; ++ ll ) if ( xx [ ll ] > xx [ ii [
0 ] ] ) ii [ 0 ] = ll ; } ii [ 0 ] -= 40 ; xx [ 0 ] = xx [ 40 + ( ii [ 0 ] )
] ; if ( xx [ 0 ] > xx [ 5 ] ) { switch ( ii [ 0 ] ) { case 0 : case 1 : {
return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:constraintViolation" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint5' kinematic constraints cannot be maintained. Check solver type and consistency tolerance in the Simscape Solver Configuration block. Check Simulink solver type and tolerances in Model Configuration Parameters. A kinematic singularity might be the source of this problem."
, neDiagMgr ) ; } } } xx [ 0 ] = xx [ 2 ] * xx [ 126 ] ; xx [ 3 ] = cos ( xx
[ 0 ] ) ; xx [ 4 ] = xx [ 2 ] * xx [ 120 ] ; xx [ 5 ] = cos ( xx [ 4 ] ) ; xx
[ 6 ] = xx [ 2 ] * xx [ 122 ] ; xx [ 7 ] = sin ( xx [ 6 ] ) ; xx [ 9 ] = cos
( xx [ 6 ] ) ; xx [ 6 ] = sin ( xx [ 4 ] ) ; xx [ 4 ] = xx [ 5 ] * xx [ 7 ] +
xx [ 9 ] * xx [ 6 ] ; xx [ 10 ] = xx [ 2 ] * xx [ 124 ] ; xx [ 11 ] = sin (
xx [ 10 ] ) ; xx [ 12 ] = cos ( xx [ 10 ] ) ; xx [ 10 ] = xx [ 5 ] * xx [ 9 ]
- xx [ 6 ] * xx [ 7 ] ; xx [ 5 ] = xx [ 4 ] * xx [ 11 ] - xx [ 12 ] * xx [ 10
] ; xx [ 6 ] = xx [ 11 ] * xx [ 10 ] + xx [ 4 ] * xx [ 12 ] ; xx [ 4 ] = sin
( xx [ 0 ] ) ; xx [ 0 ] = xx [ 3 ] * xx [ 5 ] + xx [ 6 ] * xx [ 4 ] ; xx [ 7
] = xx [ 2 ] * xx [ 116 ] ; xx [ 9 ] = sin ( xx [ 7 ] ) ; xx [ 10 ] = xx [ 2
] * xx [ 114 ] ; xx [ 11 ] = xx [ 8 ] * cos ( xx [ 10 ] ) ; xx [ 12 ] = xx [
8 ] * sin ( xx [ 10 ] ) ; xx [ 10 ] = xx [ 8 ] * ( xx [ 11 ] + xx [ 12 ] ) ;
xx [ 13 ] = ( xx [ 11 ] - xx [ 12 ] ) * xx [ 8 ] ; xx [ 8 ] = xx [ 10 ] + xx
[ 13 ] ; xx [ 11 ] = xx [ 10 ] - xx [ 13 ] ; xx [ 10 ] = cos ( xx [ 7 ] ) ;
xx [ 7 ] = xx [ 9 ] * xx [ 8 ] + xx [ 11 ] * xx [ 10 ] ; xx [ 12 ] = xx [ 2 ]
* xx [ 118 ] ; xx [ 2 ] = cos ( xx [ 12 ] ) ; xx [ 13 ] = sin ( xx [ 12 ] ) ;
xx [ 12 ] = xx [ 11 ] * xx [ 9 ] - xx [ 10 ] * xx [ 8 ] ; xx [ 8 ] = xx [ 7 ]
* xx [ 2 ] - xx [ 13 ] * xx [ 12 ] ; xx [ 9 ] = xx [ 2 ] * xx [ 12 ] + xx [ 7
] * xx [ 13 ] ; xx [ 2 ] = xx [ 4 ] * xx [ 5 ] - xx [ 6 ] * xx [ 3 ] ; xx [ 3
] = xx [ 0 ] * xx [ 8 ] - xx [ 9 ] * xx [ 2 ] ; xx [ 4 ] = xx [ 9 ] * xx [ 0
] + xx [ 8 ] * xx [ 2 ] ; state [ 0 ] = xx [ 112 ] ; state [ 1 ] = xx [ 113 ]
; state [ 2 ] = xx [ 114 ] ; state [ 3 ] = xx [ 115 ] ; state [ 4 ] = xx [
116 ] ; state [ 5 ] = xx [ 117 ] ; state [ 6 ] = xx [ 118 ] ; state [ 7 ] =
xx [ 119 ] ; state [ 8 ] = xx [ 120 ] ; state [ 9 ] = xx [ 121 ] ; state [ 10
] = xx [ 122 ] ; state [ 11 ] = xx [ 123 ] ; state [ 12 ] = xx [ 124 ] ;
state [ 13 ] = xx [ 125 ] ; state [ 14 ] = xx [ 126 ] ; state [ 15 ] = xx [
127 ] ; state [ 16 ] = xx [ 128 ] + pm_math_canonicalAngle ( xx [ 1 ] * atan2
( sqrt ( xx [ 3 ] * xx [ 3 ] ) , fabs ( - xx [ 4 ] ) ) * ( ( xx [ 4 ] * xx [
3 ] ) < 0.0 ? - 1.0 : + 1.0 ) - xx [ 128 ] ) ; state [ 17 ] = - ( xx [ 115 ]
+ xx [ 117 ] + xx [ 119 ] + xx [ 121 ] + xx [ 123 ] + xx [ 125 ] + xx [ 127 ]
) ; return NULL ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_computeConstraintError ( const
void * mech , const RuntimeDerivedValuesBundle * rtdv , const double * state
, const int * modeVector , double * error ) { const double * rtdvd = rtdv ->
mDoubles . mValues ; const int * rtdvi = rtdv -> mInts . mValues ; double xx
[ 45 ] ; ( void ) mech ; ( void ) rtdvd ; ( void ) rtdvi ; ( void )
modeVector ; xx [ 0 ] = 2.0 ; xx [ 1 ] = 1.0e-3 ; xx [ 2 ] = 0.5 ; xx [ 3 ] =
xx [ 2 ] * state [ 6 ] ; xx [ 4 ] = sin ( xx [ 3 ] ) ; xx [ 5 ] = xx [ 1 ] *
xx [ 4 ] ; xx [ 6 ] = 2.0e-3 ; xx [ 7 ] = xx [ 0 ] * xx [ 5 ] * xx [ 4 ] - xx
[ 6 ] ; xx [ 8 ] = xx [ 2 ] * state [ 4 ] ; xx [ 9 ] = sin ( xx [ 8 ] ) ; xx
[ 10 ] = xx [ 1 ] * xx [ 9 ] ; xx [ 11 ] = xx [ 0 ] * xx [ 10 ] * xx [ 9 ] -
xx [ 6 ] ; xx [ 12 ] = 0.7071067811865476 ; xx [ 13 ] = xx [ 2 ] * state [ 2
] ; xx [ 14 ] = xx [ 12 ] * cos ( xx [ 13 ] ) ; xx [ 15 ] = xx [ 12 ] * sin (
xx [ 13 ] ) ; xx [ 13 ] = xx [ 14 ] + xx [ 15 ] ; xx [ 16 ] = xx [ 12 ] * xx
[ 13 ] ; xx [ 17 ] = xx [ 14 ] - xx [ 15 ] ; xx [ 14 ] = xx [ 17 ] * xx [ 12
] ; xx [ 15 ] = xx [ 16 ] - xx [ 14 ] ; xx [ 18 ] = cos ( xx [ 8 ] ) ; xx [ 8
] = xx [ 0 ] * xx [ 18 ] * xx [ 10 ] ; xx [ 10 ] = xx [ 15 ] * xx [ 8 ] ; xx
[ 19 ] = xx [ 16 ] + xx [ 14 ] ; xx [ 14 ] = xx [ 11 ] * xx [ 15 ] ; xx [ 16
] = xx [ 17 ] * xx [ 1 ] ; xx [ 20 ] = xx [ 0 ] * xx [ 17 ] * xx [ 16 ] -
2.500000000000001e-3 ; xx [ 17 ] = xx [ 0 ] * xx [ 16 ] * xx [ 13 ] ; xx [ 13
] = xx [ 12 ] * xx [ 12 ] * xx [ 17 ] ; xx [ 16 ] = xx [ 12 ] * xx [ 20 ] *
xx [ 12 ] ; xx [ 12 ] = xx [ 9 ] * xx [ 19 ] + xx [ 15 ] * xx [ 18 ] ; xx [
21 ] = cos ( xx [ 3 ] ) ; xx [ 3 ] = xx [ 0 ] * xx [ 21 ] * xx [ 5 ] ; xx [ 5
] = xx [ 12 ] * xx [ 3 ] ; xx [ 22 ] = xx [ 15 ] * xx [ 9 ] - xx [ 18 ] * xx
[ 19 ] ; xx [ 9 ] = xx [ 7 ] * xx [ 12 ] ; xx [ 18 ] = xx [ 12 ] * xx [ 21 ]
- xx [ 4 ] * xx [ 22 ] ; xx [ 23 ] = xx [ 1 ] * xx [ 18 ] ; xx [ 24 ] = xx [
2 ] * state [ 14 ] ; xx [ 25 ] = sin ( xx [ 24 ] ) ; xx [ 26 ] = xx [ 1 ] *
xx [ 25 ] ; xx [ 27 ] = xx [ 6 ] - xx [ 0 ] * xx [ 26 ] * xx [ 25 ] ; xx [ 28
] = xx [ 2 ] * state [ 12 ] ; xx [ 29 ] = sin ( xx [ 28 ] ) ; xx [ 30 ] = xx
[ 1 ] * xx [ 29 ] ; xx [ 31 ] = xx [ 6 ] - xx [ 0 ] * xx [ 30 ] * xx [ 29 ] ;
xx [ 32 ] = xx [ 2 ] * state [ 8 ] ; xx [ 33 ] = cos ( xx [ 32 ] ) ; xx [ 34
] = xx [ 2 ] * state [ 10 ] ; xx [ 2 ] = sin ( xx [ 34 ] ) ; xx [ 35 ] = cos
( xx [ 34 ] ) ; xx [ 34 ] = sin ( xx [ 32 ] ) ; xx [ 32 ] = xx [ 33 ] * xx [
2 ] + xx [ 35 ] * xx [ 34 ] ; xx [ 36 ] = cos ( xx [ 28 ] ) ; xx [ 28 ] = xx
[ 0 ] * xx [ 36 ] * xx [ 30 ] ; xx [ 30 ] = xx [ 32 ] * xx [ 28 ] ; xx [ 37 ]
= xx [ 33 ] * xx [ 35 ] - xx [ 34 ] * xx [ 2 ] ; xx [ 38 ] = xx [ 32 ] * xx [
31 ] ; xx [ 39 ] = xx [ 1 ] * xx [ 2 ] ; xx [ 40 ] = xx [ 6 ] - xx [ 0 ] * xx
[ 39 ] * xx [ 2 ] ; xx [ 2 ] = xx [ 0 ] * xx [ 35 ] * xx [ 39 ] ; xx [ 6 ] =
xx [ 2 ] * xx [ 34 ] ; xx [ 35 ] = xx [ 40 ] * xx [ 34 ] ; xx [ 39 ] = xx [ 1
] * xx [ 34 ] ; xx [ 41 ] = xx [ 29 ] * xx [ 37 ] + xx [ 32 ] * xx [ 36 ] ;
xx [ 42 ] = xx [ 41 ] * xx [ 27 ] ; xx [ 43 ] = cos ( xx [ 24 ] ) ; xx [ 24 ]
= xx [ 0 ] * xx [ 43 ] * xx [ 26 ] ; xx [ 26 ] = xx [ 41 ] * xx [ 24 ] ; xx [
44 ] = xx [ 32 ] * xx [ 29 ] - xx [ 36 ] * xx [ 37 ] ; xx [ 29 ] = xx [ 25 ]
* xx [ 44 ] - xx [ 41 ] * xx [ 43 ] ; xx [ 36 ] = xx [ 1 ] * xx [ 29 ] ;
error [ 0 ] = xx [ 7 ] + xx [ 11 ] + xx [ 0 ] * ( xx [ 10 ] * xx [ 19 ] - xx
[ 15 ] * xx [ 14 ] ) + xx [ 20 ] - xx [ 0 ] * ( xx [ 13 ] + xx [ 16 ] ) +
state [ 0 ] - ( xx [ 5 ] * xx [ 22 ] + xx [ 12 ] * xx [ 9 ] ) * xx [ 0 ] + xx
[ 0 ] * xx [ 23 ] * xx [ 18 ] - ( xx [ 27 ] + xx [ 31 ] - xx [ 0 ] * ( xx [
30 ] * xx [ 37 ] + xx [ 32 ] * xx [ 38 ] ) + xx [ 40 ] - xx [ 0 ] * ( xx [ 33
] * xx [ 6 ] + xx [ 35 ] * xx [ 34 ] ) - xx [ 0 ] * xx [ 39 ] * xx [ 34 ] - (
xx [ 41 ] * xx [ 42 ] - xx [ 26 ] * xx [ 44 ] ) * xx [ 0 ] - xx [ 0 ] * xx [
36 ] * xx [ 29 ] ) + 0.011 ; error [ 1 ] = xx [ 3 ] + xx [ 0 ] * ( xx [ 9 ] *
xx [ 22 ] - xx [ 12 ] * xx [ 5 ] ) + xx [ 8 ] - ( xx [ 14 ] * xx [ 19 ] + xx
[ 15 ] * xx [ 10 ] ) * xx [ 0 ] - ( xx [ 17 ] + ( xx [ 16 ] - xx [ 13 ] ) *
xx [ 0 ] ) - xx [ 0 ] * ( xx [ 21 ] * xx [ 22 ] + xx [ 12 ] * xx [ 4 ] ) * xx
[ 23 ] - ( xx [ 0 ] * ( xx [ 43 ] * xx [ 44 ] + xx [ 41 ] * xx [ 25 ] ) * xx
[ 36 ] + ( xx [ 38 ] * xx [ 37 ] - xx [ 32 ] * xx [ 30 ] ) * xx [ 0 ] + ( xx
[ 33 ] * xx [ 35 ] - xx [ 6 ] * xx [ 34 ] ) * xx [ 0 ] + xx [ 2 ] + xx [ 0 ]
* xx [ 33 ] * xx [ 39 ] + xx [ 28 ] - xx [ 0 ] * ( xx [ 41 ] * xx [ 26 ] + xx
[ 42 ] * xx [ 44 ] ) + xx [ 24 ] ) - 1.500000000000001e-3 ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_resetModeVector ( const void *
mech , int * modeVector ) { ( void ) mech ; ( void ) modeVector ; } boolean_T
MagneticBasket_Simscape_Optimizer_dda62cd9_1_hasJointDisToNormModeChange (
const void * mech , const int * prevModeVector , const int * modeVector ) { (
void ) mech ; ( void ) prevModeVector ; ( void ) modeVector ; return 0 ; }
PmfMessageId
MagneticBasket_Simscape_Optimizer_dda62cd9_1_performJointDisToNormModeChange
( const void * mech , const RuntimeDerivedValuesBundle * rtdv , const int *
eqnEnableFlags , const int * prevModeVector , const int * modeVector , const
double * input , double * state , void * neDiagMgr0 ) { const double * rtdvd
= rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts . mValues ;
NeuDiagnosticManager * neDiagMgr = ( NeuDiagnosticManager * ) neDiagMgr0 ; (
void ) mech ; ( void ) rtdvd ; ( void ) rtdvi ; ( void ) eqnEnableFlags ; (
void ) prevModeVector ; ( void ) modeVector ; ( void ) input ; ( void ) state
; ( void ) neDiagMgr ; return NULL ; } void
MagneticBasket_Simscape_Optimizer_dda62cd9_1_onModeChangedCutJoints ( const
void * mech , const int * prevModeVector , const int * modeVector , double *
state ) { ( void ) mech ; ( void ) prevModeVector ; ( void ) modeVector ; (
void ) state ; }
