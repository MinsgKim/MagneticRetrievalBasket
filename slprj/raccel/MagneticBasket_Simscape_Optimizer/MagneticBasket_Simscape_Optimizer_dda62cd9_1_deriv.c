#include <math.h>
#include <string.h>
#include "pm_std.h"
#include "sm_std.h"
#include "ne_std.h"
#include "ne_dae.h"
#include "sm_ssci_run_time_errors.h"
#include "sm_RuntimeDerivedValuesBundle.h"
#include "MagneticBasket_Simscape_Optimizer_dda62cd9_1_geometries.h"
PmfMessageId MagneticBasket_Simscape_Optimizer_dda62cd9_1_compDerivs ( const
RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags , const double
* state , const int * modeVector , const double * input , const double *
inputDot , const double * inputDdot , const double * discreteState , double *
deriv , double * errorResult , NeuDiagnosticManager * neDiagMgr ) { const
double * rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv ->
mInts . mValues ; int ii [ 2 ] ; double xx [ 245 ] ; ( void ) rtdvd ; ( void
) rtdvi ; ( void ) eqnEnableFlags ; ( void ) modeVector ; ( void ) inputDot ;
( void ) inputDdot ; ( void ) discreteState ; ( void ) neDiagMgr ; xx [ 0 ] =
1.0 ; xx [ 1 ] = 3.000000000000001e-9 ; xx [ 2 ] = 2.0 ; xx [ 3 ] = 0.5 ; xx
[ 4 ] = xx [ 3 ] * state [ 6 ] ; xx [ 5 ] = cos ( xx [ 4 ] ) ; xx [ 6 ] =
1.0e-3 ; xx [ 7 ] = sin ( xx [ 4 ] ) ; xx [ 4 ] = xx [ 6 ] * xx [ 7 ] ; xx [
8 ] = xx [ 5 ] * xx [ 4 ] ; xx [ 9 ] = xx [ 2 ] * xx [ 8 ] ; xx [ 10 ] = xx [
3 ] * state [ 4 ] ; xx [ 11 ] = sin ( xx [ 10 ] ) ; xx [ 12 ] =
0.7071067811865476 ; xx [ 13 ] = xx [ 3 ] * state [ 2 ] ; xx [ 14 ] = xx [ 12
] * cos ( xx [ 13 ] ) ; xx [ 15 ] = xx [ 12 ] * sin ( xx [ 13 ] ) ; xx [ 13 ]
= xx [ 14 ] + xx [ 15 ] ; xx [ 16 ] = xx [ 12 ] * xx [ 13 ] ; xx [ 17 ] = xx
[ 14 ] - xx [ 15 ] ; xx [ 14 ] = xx [ 17 ] * xx [ 12 ] ; xx [ 15 ] = xx [ 16
] + xx [ 14 ] ; xx [ 18 ] = xx [ 16 ] - xx [ 14 ] ; xx [ 14 ] = cos ( xx [ 10
] ) ; xx [ 10 ] = xx [ 11 ] * xx [ 15 ] + xx [ 18 ] * xx [ 14 ] ; xx [ 16 ] =
xx [ 6 ] - xx [ 2 ] * xx [ 4 ] * xx [ 7 ] ; xx [ 4 ] = xx [ 10 ] * xx [ 16 ]
; xx [ 19 ] = xx [ 18 ] * xx [ 11 ] - xx [ 14 ] * xx [ 15 ] ; xx [ 20 ] = xx
[ 10 ] * xx [ 9 ] ; xx [ 21 ] = xx [ 5 ] * xx [ 19 ] + xx [ 10 ] * xx [ 7 ] ;
xx [ 22 ] = 1.0e-3 ; xx [ 23 ] = xx [ 10 ] * xx [ 5 ] - xx [ 7 ] * xx [ 19 ]
; xx [ 24 ] = xx [ 22 ] * xx [ 23 ] ; xx [ 25 ] = xx [ 2 ] * xx [ 21 ] * xx [
24 ] ; xx [ 26 ] = xx [ 9 ] - ( xx [ 4 ] * xx [ 19 ] + xx [ 10 ] * xx [ 20 ]
) * xx [ 2 ] - xx [ 25 ] ; xx [ 9 ] = 4.062500000000001e-12 ; xx [ 27 ] = xx
[ 26 ] / xx [ 9 ] ; xx [ 28 ] = xx [ 1 ] * xx [ 27 ] ; xx [ 29 ] = xx [ 28 ]
* xx [ 7 ] ; xx [ 30 ] = xx [ 28 ] - xx [ 2 ] * xx [ 29 ] * xx [ 7 ] ; xx [
28 ] = 3.0e-6 ; xx [ 31 ] = xx [ 5 ] * xx [ 7 ] ; xx [ 32 ] = xx [ 2 ] * xx [
31 ] ; xx [ 33 ] = xx [ 28 ] * xx [ 32 ] ; xx [ 34 ] = 7.846153846153842e-7 ;
xx [ 35 ] = xx [ 5 ] * xx [ 5 ] ; xx [ 36 ] = xx [ 2 ] * xx [ 35 ] - xx [ 0 ]
; xx [ 37 ] = xx [ 34 ] * xx [ 36 ] ; xx [ 38 ] = xx [ 33 ] * xx [ 32 ] + xx
[ 37 ] * xx [ 36 ] ; xx [ 39 ] = xx [ 28 ] + xx [ 38 ] ; xx [ 40 ] =
7.846153846153845e-10 ; xx [ 41 ] = ( xx [ 35 ] + xx [ 7 ] * xx [ 7 ] ) * xx
[ 2 ] - xx [ 0 ] ; xx [ 35 ] = xx [ 40 ] * xx [ 36 ] * xx [ 41 ] ; xx [ 42 ]
= xx [ 22 ] * xx [ 7 ] ; xx [ 43 ] = 2.0e-3 ; xx [ 44 ] = xx [ 2 ] * xx [ 42
] * xx [ 7 ] - xx [ 43 ] ; xx [ 45 ] = xx [ 2 ] * xx [ 5 ] * xx [ 42 ] ; xx [
42 ] = xx [ 34 ] * xx [ 32 ] ; xx [ 46 ] = xx [ 28 ] * xx [ 36 ] ; xx [ 47 ]
= xx [ 42 ] * xx [ 36 ] - xx [ 46 ] * xx [ 32 ] ; xx [ 48 ] = xx [ 44 ] * xx
[ 38 ] - xx [ 45 ] * xx [ 47 ] ; xx [ 38 ] = xx [ 35 ] + xx [ 48 ] ; xx [ 49
] = xx [ 39 ] * xx [ 6 ] - xx [ 38 ] ; xx [ 50 ] = xx [ 44 ] * xx [ 11 ] ; xx
[ 51 ] = xx [ 45 ] * xx [ 11 ] ; xx [ 52 ] = xx [ 45 ] - ( xx [ 14 ] * xx [
50 ] + xx [ 51 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 53 ] = xx [ 6 ] * xx [ 11 ]
; xx [ 54 ] = xx [ 14 ] * xx [ 53 ] ; xx [ 55 ] = xx [ 52 ] + xx [ 2 ] * xx [
54 ] ; xx [ 56 ] = xx [ 50 ] * xx [ 11 ] ; xx [ 50 ] = xx [ 14 ] * xx [ 51 ]
; xx [ 51 ] = xx [ 2 ] * xx [ 53 ] * xx [ 11 ] ; xx [ 53 ] = xx [ 2 ] * ( xx
[ 56 ] - xx [ 50 ] ) - ( xx [ 44 ] + xx [ 51 ] ) + xx [ 6 ] ; xx [ 57 ] = xx
[ 18 ] * xx [ 53 ] ; xx [ 58 ] = xx [ 18 ] * xx [ 55 ] ; xx [ 59 ] = xx [ 55
] + xx [ 2 ] * ( xx [ 57 ] * xx [ 15 ] - xx [ 18 ] * xx [ 58 ] ) - xx [ 25 ]
; xx [ 25 ] = 1.0625e-12 ; xx [ 55 ] = xx [ 2 ] * xx [ 5 ] * xx [ 29 ] ; xx [
29 ] = xx [ 25 ] * xx [ 27 ] + xx [ 45 ] * xx [ 55 ] - xx [ 44 ] * xx [ 30 ]
; xx [ 60 ] = 7.846153846153846e-13 ; xx [ 61 ] = xx [ 40 ] * xx [ 32 ] * xx
[ 41 ] ; xx [ 62 ] = xx [ 45 ] * xx [ 61 ] - xx [ 44 ] * xx [ 35 ] ; xx [ 35
] = xx [ 37 ] * xx [ 32 ] - xx [ 33 ] * xx [ 36 ] ; xx [ 33 ] = xx [ 46 ] *
xx [ 36 ] + xx [ 42 ] * xx [ 32 ] ; xx [ 37 ] = xx [ 44 ] * xx [ 35 ] - xx [
33 ] * xx [ 45 ] ; xx [ 42 ] = xx [ 25 ] + xx [ 60 ] * xx [ 41 ] * xx [ 41 ]
- xx [ 62 ] - xx [ 62 ] - ( xx [ 45 ] * xx [ 37 ] - xx [ 44 ] * xx [ 48 ] ) ;
xx [ 41 ] = xx [ 42 ] - xx [ 6 ] * xx [ 38 ] ; xx [ 46 ] = xx [ 41 ] + xx [
49 ] * xx [ 6 ] ; ii [ 0 ] = factorSymmetricPosDef ( xx + 46 , 1 , xx + 48 )
; if ( ii [ 0 ] != 0 ) { return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassBase" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint7' has a degenerate mass distribution on its base side."
, neDiagMgr ) ; } xx [ 48 ] = ( xx [ 59 ] - ( xx [ 29 ] + xx [ 6 ] * xx [ 30
] ) ) / xx [ 46 ] ; xx [ 62 ] = xx [ 30 ] + xx [ 49 ] * xx [ 48 ] ; xx [ 30 ]
= xx [ 61 ] + xx [ 37 ] ; xx [ 37 ] = xx [ 6 ] * xx [ 35 ] - xx [ 30 ] ; xx [
61 ] = xx [ 55 ] + xx [ 37 ] * xx [ 48 ] ; xx [ 55 ] = xx [ 61 ] * xx [ 11 ]
; xx [ 63 ] = xx [ 62 ] * xx [ 11 ] ; xx [ 64 ] = xx [ 62 ] - ( xx [ 14 ] *
xx [ 55 ] + xx [ 63 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 62 ] = xx [ 14 ] * xx [
14 ] ; xx [ 65 ] = ( xx [ 62 ] + xx [ 11 ] * xx [ 11 ] ) * xx [ 2 ] - xx [ 0
] ; xx [ 66 ] = xx [ 14 ] * xx [ 11 ] ; xx [ 67 ] = xx [ 2 ] * xx [ 66 ] ; xx
[ 68 ] = xx [ 37 ] / xx [ 46 ] ; xx [ 69 ] = xx [ 30 ] + xx [ 41 ] * xx [ 68
] ; xx [ 70 ] = xx [ 2 ] * xx [ 62 ] - xx [ 0 ] ; xx [ 62 ] = xx [ 49 ] / xx
[ 46 ] ; xx [ 71 ] = xx [ 38 ] + xx [ 41 ] * xx [ 62 ] ; xx [ 72 ] = xx [ 65
] * ( xx [ 67 ] * xx [ 69 ] - xx [ 70 ] * xx [ 71 ] ) ; xx [ 73 ] = xx [ 22 ]
* xx [ 11 ] ; xx [ 74 ] = xx [ 2 ] * xx [ 73 ] * xx [ 11 ] - xx [ 43 ] ; xx [
75 ] = xx [ 39 ] - xx [ 49 ] * xx [ 62 ] ; xx [ 76 ] = xx [ 49 ] * xx [ 68 ]
; xx [ 77 ] = xx [ 47 ] - xx [ 76 ] ; xx [ 78 ] = xx [ 75 ] * xx [ 70 ] - xx
[ 67 ] * xx [ 77 ] ; xx [ 79 ] = xx [ 35 ] - xx [ 76 ] ; xx [ 76 ] = xx [ 28
] + xx [ 33 ] ; xx [ 33 ] = xx [ 76 ] - xx [ 37 ] * xx [ 68 ] ; xx [ 80 ] =
xx [ 70 ] * xx [ 79 ] - xx [ 67 ] * xx [ 33 ] ; xx [ 81 ] = xx [ 70 ] * xx [
78 ] - xx [ 67 ] * xx [ 80 ] ; xx [ 82 ] = xx [ 2 ] * xx [ 14 ] * xx [ 73 ] ;
xx [ 73 ] = xx [ 70 ] * xx [ 77 ] + xx [ 67 ] * xx [ 75 ] ; xx [ 75 ] = xx [
33 ] * xx [ 70 ] + xx [ 67 ] * xx [ 79 ] ; xx [ 33 ] = xx [ 73 ] * xx [ 70 ]
- xx [ 75 ] * xx [ 67 ] ; xx [ 77 ] = xx [ 74 ] * xx [ 81 ] - xx [ 82 ] * xx
[ 33 ] ; xx [ 79 ] = xx [ 72 ] - xx [ 77 ] ; xx [ 83 ] = xx [ 28 ] + xx [ 81
] ; xx [ 81 ] = xx [ 79 ] + xx [ 83 ] * xx [ 6 ] ; xx [ 84 ] = xx [ 14 ] * xx
[ 7 ] + xx [ 5 ] * xx [ 11 ] ; xx [ 85 ] = xx [ 14 ] * xx [ 5 ] - xx [ 11 ] *
xx [ 7 ] ; xx [ 86 ] = xx [ 84 ] * xx [ 15 ] + xx [ 18 ] * xx [ 85 ] ; xx [
87 ] = xx [ 86 ] * xx [ 22 ] ; xx [ 88 ] = xx [ 52 ] + xx [ 82 ] ; xx [ 52 ]
= xx [ 17 ] * xx [ 88 ] ; xx [ 89 ] = xx [ 44 ] + xx [ 2 ] * ( xx [ 50 ] - xx
[ 56 ] ) + xx [ 74 ] ; xx [ 50 ] = xx [ 89 ] * xx [ 17 ] ; xx [ 56 ] = xx [
17 ] * xx [ 6 ] ; xx [ 90 ] = xx [ 88 ] - ( xx [ 17 ] * xx [ 52 ] - xx [ 50 ]
* xx [ 13 ] ) * xx [ 2 ] - xx [ 2 ] * xx [ 56 ] * xx [ 13 ] ; xx [ 88 ] = xx
[ 2 ] * ( xx [ 17 ] * xx [ 50 ] + xx [ 52 ] * xx [ 13 ] ) - ( xx [ 89 ] + xx
[ 2 ] * xx [ 17 ] * xx [ 56 ] ) + xx [ 6 ] ; xx [ 50 ] = xx [ 12 ] * xx [ 88
] * xx [ 12 ] ; xx [ 52 ] = xx [ 12 ] * xx [ 90 ] * xx [ 12 ] ; xx [ 56 ] =
xx [ 2 ] * xx [ 87 ] * ( xx [ 85 ] * xx [ 15 ] - xx [ 84 ] * xx [ 18 ] ) + xx
[ 90 ] + xx [ 2 ] * ( xx [ 50 ] - xx [ 52 ] ) ; xx [ 84 ] = xx [ 61 ] + xx [
2 ] * ( xx [ 14 ] * xx [ 63 ] - xx [ 55 ] * xx [ 11 ] ) ; xx [ 55 ] = xx [ 41
] / xx [ 46 ] ; xx [ 61 ] = ( xx [ 70 ] * xx [ 69 ] + xx [ 67 ] * xx [ 71 ] )
* xx [ 65 ] ; xx [ 63 ] = xx [ 74 ] * xx [ 72 ] + xx [ 61 ] * xx [ 82 ] ; xx
[ 69 ] = xx [ 70 ] * xx [ 80 ] + xx [ 67 ] * xx [ 78 ] ; xx [ 71 ] = xx [ 75
] * xx [ 70 ] + xx [ 73 ] * xx [ 67 ] ; xx [ 72 ] = xx [ 74 ] * xx [ 69 ] -
xx [ 71 ] * xx [ 82 ] ; xx [ 73 ] = ( xx [ 42 ] - xx [ 41 ] * xx [ 55 ] ) *
xx [ 65 ] * xx [ 65 ] - xx [ 63 ] - xx [ 63 ] - ( xx [ 82 ] * xx [ 72 ] - xx
[ 74 ] * xx [ 77 ] ) + xx [ 6 ] * xx [ 79 ] + xx [ 81 ] * xx [ 6 ] + xx [ 25
] ; ii [ 0 ] = factorSymmetricPosDef ( xx + 73 , 1 , xx + 42 ) ; if ( ii [ 0
] != 0 ) { return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassBase" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint8' has a degenerate mass distribution on its base side."
, neDiagMgr ) ; } xx [ 42 ] = ( xx [ 56 ] - ( xx [ 29 ] + xx [ 41 ] * xx [ 48
] + xx [ 84 ] * xx [ 82 ] - xx [ 74 ] * xx [ 64 ] + xx [ 6 ] * xx [ 64 ] ) )
/ xx [ 73 ] ; xx [ 29 ] = xx [ 64 ] + xx [ 81 ] * xx [ 42 ] ; xx [ 63 ] = xx
[ 61 ] + xx [ 72 ] ; xx [ 61 ] = xx [ 69 ] * xx [ 6 ] - xx [ 63 ] ; xx [ 64 ]
= xx [ 2 ] * xx [ 13 ] * xx [ 13 ] - xx [ 0 ] ; xx [ 65 ] = xx [ 81 ] / xx [
73 ] ; xx [ 72 ] = xx [ 17 ] * xx [ 13 ] ; xx [ 75 ] = xx [ 2 ] * xx [ 72 ] ;
xx [ 77 ] = xx [ 61 ] / xx [ 73 ] ; xx [ 78 ] = xx [ 81 ] * xx [ 77 ] ; xx [
79 ] = xx [ 28 ] + xx [ 71 ] ; xx [ 71 ] = 4.500000000000001e-6 + xx [ 64 ] *
( ( xx [ 83 ] - xx [ 81 ] * xx [ 65 ] ) * xx [ 64 ] + xx [ 75 ] * ( xx [ 33 ]
- xx [ 78 ] ) ) + xx [ 75 ] * ( ( xx [ 69 ] - xx [ 78 ] ) * xx [ 64 ] + xx [
75 ] * ( xx [ 79 ] - xx [ 61 ] * xx [ 77 ] ) ) ; ii [ 0 ] =
factorSymmetricPosDef ( xx + 71 , 1 , xx + 69 ) ; if ( ii [ 0 ] != 0 ) {
return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Prismatic Joint1' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 69 ] = ( xx [ 0 ] - ( xx [ 29 ] - ( xx [ 17 ] * xx [
17 ] * xx [ 29 ] - ( xx [ 84 ] + xx [ 61 ] * xx [ 42 ] ) * xx [ 17 ] * xx [
13 ] ) * xx [ 2 ] ) ) / xx [ 71 ] ; xx [ 29 ] = xx [ 17 ] * xx [ 69 ] ; xx [
75 ] = xx [ 69 ] - xx [ 2 ] * xx [ 17 ] * xx [ 29 ] ; xx [ 78 ] = xx [ 2 ] *
xx [ 29 ] * xx [ 13 ] ; xx [ 29 ] = xx [ 42 ] - ( xx [ 65 ] * xx [ 75 ] + xx
[ 77 ] * xx [ 78 ] ) ; xx [ 42 ] = xx [ 82 ] * xx [ 29 ] + xx [ 78 ] ; xx [
78 ] = xx [ 75 ] + xx [ 6 ] * xx [ 29 ] - xx [ 74 ] * xx [ 29 ] ; xx [ 75 ] =
xx [ 11 ] * xx [ 78 ] ; xx [ 80 ] = xx [ 11 ] * xx [ 42 ] ; xx [ 83 ] = xx [
42 ] - ( xx [ 14 ] * xx [ 75 ] + xx [ 80 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 42
] = xx [ 78 ] + xx [ 2 ] * ( xx [ 14 ] * xx [ 80 ] - xx [ 75 ] * xx [ 11 ] )
; xx [ 75 ] = xx [ 48 ] - ( xx [ 55 ] * xx [ 29 ] + xx [ 68 ] * xx [ 83 ] +
xx [ 42 ] * xx [ 62 ] ) ; xx [ 48 ] = xx [ 29 ] + xx [ 75 ] ; xx [ 78 ] =
0.2615384615384614 ; xx [ 80 ] = xx [ 42 ] + xx [ 6 ] * xx [ 75 ] - xx [ 44 ]
* xx [ 48 ] ; xx [ 42 ] = 738.4615384615385 ; xx [ 84 ] = xx [ 3 ] * state [
14 ] ; xx [ 85 ] = sin ( xx [ 84 ] ) ; xx [ 89 ] = xx [ 22 ] * xx [ 85 ] ; xx
[ 90 ] = xx [ 43 ] - xx [ 2 ] * xx [ 89 ] * xx [ 85 ] ; xx [ 91 ] = xx [ 3 ]
* state [ 12 ] ; xx [ 92 ] = sin ( xx [ 91 ] ) ; xx [ 93 ] = xx [ 3 ] * state
[ 8 ] ; xx [ 94 ] = cos ( xx [ 93 ] ) ; xx [ 95 ] = xx [ 3 ] * state [ 10 ] ;
xx [ 3 ] = cos ( xx [ 95 ] ) ; xx [ 96 ] = sin ( xx [ 93 ] ) ; xx [ 93 ] =
sin ( xx [ 95 ] ) ; xx [ 95 ] = xx [ 94 ] * xx [ 3 ] - xx [ 96 ] * xx [ 93 ]
; xx [ 97 ] = xx [ 94 ] * xx [ 93 ] + xx [ 3 ] * xx [ 96 ] ; xx [ 98 ] = cos
( xx [ 91 ] ) ; xx [ 91 ] = xx [ 92 ] * xx [ 95 ] + xx [ 97 ] * xx [ 98 ] ;
xx [ 99 ] = xx [ 6 ] * xx [ 85 ] ; xx [ 100 ] = xx [ 6 ] - xx [ 2 ] * xx [ 99
] * xx [ 85 ] ; xx [ 101 ] = xx [ 91 ] * xx [ 100 ] ; xx [ 102 ] = xx [ 97 ]
* xx [ 92 ] - xx [ 98 ] * xx [ 95 ] ; xx [ 103 ] = cos ( xx [ 84 ] ) ; xx [
84 ] = xx [ 103 ] * xx [ 99 ] ; xx [ 99 ] = xx [ 2 ] * xx [ 84 ] ; xx [ 104 ]
= xx [ 91 ] * xx [ 99 ] ; xx [ 105 ] = xx [ 103 ] * xx [ 102 ] + xx [ 91 ] *
xx [ 85 ] ; xx [ 106 ] = xx [ 85 ] * xx [ 102 ] - xx [ 91 ] * xx [ 103 ] ; xx
[ 107 ] = xx [ 22 ] * xx [ 106 ] ; xx [ 108 ] = xx [ 2 ] * xx [ 105 ] * xx [
107 ] ; xx [ 109 ] = ( xx [ 101 ] * xx [ 102 ] + xx [ 91 ] * xx [ 104 ] ) *
xx [ 2 ] - xx [ 99 ] - xx [ 108 ] ; xx [ 99 ] = xx [ 109 ] / xx [ 9 ] ; xx [
110 ] = xx [ 1 ] * xx [ 99 ] ; xx [ 111 ] = xx [ 110 ] * xx [ 85 ] ; xx [ 112
] = xx [ 2 ] * xx [ 111 ] * xx [ 85 ] - xx [ 110 ] ; xx [ 110 ] = xx [ 2 ] *
xx [ 103 ] * xx [ 89 ] ; xx [ 89 ] = xx [ 2 ] * xx [ 103 ] * xx [ 111 ] ; xx
[ 111 ] = xx [ 90 ] * xx [ 112 ] - xx [ 110 ] * xx [ 89 ] - xx [ 25 ] * xx [
99 ] ; xx [ 113 ] = xx [ 90 ] * xx [ 92 ] ; xx [ 114 ] = xx [ 110 ] * xx [ 92
] ; xx [ 115 ] = ( xx [ 98 ] * xx [ 113 ] - xx [ 114 ] * xx [ 92 ] ) * xx [ 2
] ; xx [ 116 ] = xx [ 6 ] * xx [ 92 ] ; xx [ 117 ] = xx [ 98 ] * xx [ 116 ] ;
xx [ 118 ] = xx [ 110 ] + xx [ 115 ] + xx [ 2 ] * xx [ 117 ] ; xx [ 119 ] =
xx [ 90 ] - xx [ 2 ] * ( xx [ 98 ] * xx [ 114 ] + xx [ 113 ] * xx [ 92 ] ) ;
xx [ 113 ] = xx [ 2 ] * xx [ 116 ] * xx [ 92 ] ; xx [ 114 ] = xx [ 119 ] - xx
[ 113 ] + xx [ 6 ] ; xx [ 116 ] = xx [ 97 ] * xx [ 114 ] ; xx [ 120 ] = xx [
97 ] * xx [ 118 ] ; xx [ 121 ] = xx [ 118 ] + ( xx [ 116 ] * xx [ 95 ] - xx [
97 ] * xx [ 120 ] ) * xx [ 2 ] + xx [ 108 ] ; xx [ 108 ] = xx [ 103 ] * xx [
85 ] ; xx [ 118 ] = xx [ 2 ] * xx [ 108 ] ; xx [ 122 ] = xx [ 28 ] * xx [ 118
] ; xx [ 123 ] = xx [ 103 ] * xx [ 103 ] ; xx [ 124 ] = xx [ 2 ] * xx [ 123 ]
- xx [ 0 ] ; xx [ 125 ] = xx [ 34 ] * xx [ 124 ] ; xx [ 126 ] = xx [ 122 ] *
xx [ 118 ] + xx [ 125 ] * xx [ 124 ] ; xx [ 127 ] = xx [ 28 ] + xx [ 126 ] ;
xx [ 128 ] = ( xx [ 123 ] + xx [ 85 ] * xx [ 85 ] ) * xx [ 2 ] - xx [ 0 ] ;
xx [ 123 ] = xx [ 40 ] * xx [ 124 ] * xx [ 128 ] ; xx [ 129 ] = xx [ 28 ] *
xx [ 124 ] ; xx [ 130 ] = xx [ 34 ] * xx [ 118 ] ; xx [ 34 ] = xx [ 129 ] *
xx [ 118 ] - xx [ 130 ] * xx [ 124 ] ; xx [ 131 ] = xx [ 126 ] * xx [ 90 ] -
xx [ 110 ] * xx [ 34 ] ; xx [ 126 ] = xx [ 123 ] - xx [ 131 ] ; xx [ 132 ] =
xx [ 127 ] * xx [ 6 ] - xx [ 126 ] ; xx [ 133 ] = xx [ 40 ] * xx [ 118 ] * xx
[ 128 ] ; xx [ 40 ] = xx [ 90 ] * xx [ 123 ] + xx [ 110 ] * xx [ 133 ] ; xx [
123 ] = xx [ 122 ] * xx [ 124 ] - xx [ 125 ] * xx [ 118 ] ; xx [ 122 ] = xx [
129 ] * xx [ 124 ] + xx [ 130 ] * xx [ 118 ] ; xx [ 125 ] = xx [ 90 ] * xx [
123 ] - xx [ 122 ] * xx [ 110 ] ; xx [ 129 ] = xx [ 25 ] + xx [ 60 ] * xx [
128 ] * xx [ 128 ] - xx [ 40 ] - xx [ 40 ] + xx [ 131 ] * xx [ 90 ] - xx [
125 ] * xx [ 110 ] ; xx [ 40 ] = xx [ 6 ] * xx [ 126 ] - xx [ 129 ] ; xx [ 60
] = xx [ 6 ] * xx [ 132 ] - xx [ 40 ] ; ii [ 0 ] = factorSymmetricPosDef ( xx
+ 60 , 1 , xx + 128 ) ; if ( ii [ 0 ] != 0 ) { return
sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint3' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 128 ] = ( xx [ 111 ] + xx [ 6 ] * xx [ 112 ] - xx [
121 ] ) / xx [ 60 ] ; xx [ 130 ] = xx [ 22 ] * xx [ 92 ] ; xx [ 131 ] = xx [
43 ] - xx [ 2 ] * xx [ 130 ] * xx [ 92 ] ; xx [ 134 ] = xx [ 112 ] - xx [ 128
] * xx [ 132 ] ; xx [ 112 ] = xx [ 133 ] + xx [ 125 ] ; xx [ 125 ] = xx [ 112
] + xx [ 6 ] * xx [ 123 ] ; xx [ 133 ] = xx [ 89 ] - xx [ 125 ] * xx [ 128 ]
; xx [ 89 ] = xx [ 92 ] * xx [ 133 ] ; xx [ 135 ] = xx [ 92 ] * xx [ 134 ] ;
xx [ 136 ] = xx [ 134 ] + xx [ 2 ] * ( xx [ 98 ] * xx [ 89 ] - xx [ 135 ] *
xx [ 92 ] ) ; xx [ 134 ] = xx [ 2 ] * xx [ 98 ] * xx [ 130 ] ; xx [ 130 ] =
xx [ 133 ] - ( xx [ 98 ] * xx [ 135 ] + xx [ 89 ] * xx [ 92 ] ) * xx [ 2 ] ;
xx [ 89 ] = xx [ 111 ] + xx [ 128 ] * xx [ 40 ] + xx [ 131 ] * xx [ 136 ] -
xx [ 134 ] * xx [ 130 ] ; xx [ 111 ] = xx [ 115 ] + xx [ 110 ] + xx [ 134 ] ;
xx [ 115 ] = xx [ 93 ] * xx [ 111 ] ; xx [ 133 ] = xx [ 115 ] * xx [ 93 ] ;
xx [ 135 ] = xx [ 119 ] + xx [ 131 ] ; xx [ 119 ] = xx [ 135 ] * xx [ 93 ] ;
xx [ 137 ] = xx [ 3 ] * xx [ 119 ] ; xx [ 138 ] = xx [ 6 ] * xx [ 93 ] ; xx [
139 ] = xx [ 3 ] * xx [ 138 ] ; xx [ 140 ] = xx [ 2 ] * ( xx [ 133 ] - xx [
137 ] ) - xx [ 111 ] - xx [ 2 ] * xx [ 139 ] ; xx [ 141 ] = ( xx [ 3 ] * xx [
115 ] + xx [ 119 ] * xx [ 93 ] ) * xx [ 2 ] ; xx [ 115 ] = xx [ 2 ] * xx [
138 ] * xx [ 93 ] ; xx [ 119 ] = xx [ 135 ] - ( xx [ 141 ] + xx [ 115 ] ) +
xx [ 6 ] ; xx [ 138 ] = xx [ 119 ] * xx [ 96 ] ; xx [ 142 ] = xx [ 96 ] * xx
[ 140 ] ; xx [ 143 ] = xx [ 98 ] * xx [ 85 ] + xx [ 103 ] * xx [ 92 ] ; xx [
144 ] = xx [ 98 ] * xx [ 103 ] - xx [ 92 ] * xx [ 85 ] ; xx [ 145 ] = xx [
143 ] * xx [ 95 ] + xx [ 97 ] * xx [ 144 ] ; xx [ 146 ] = xx [ 145 ] * xx [
22 ] ; xx [ 147 ] = xx [ 140 ] - ( xx [ 94 ] * xx [ 138 ] + xx [ 142 ] * xx [
96 ] ) * xx [ 2 ] - xx [ 2 ] * xx [ 146 ] * ( xx [ 95 ] * xx [ 144 ] - xx [
97 ] * xx [ 143 ] ) ; xx [ 140 ] = xx [ 98 ] * xx [ 92 ] ; xx [ 148 ] = xx [
2 ] * xx [ 140 ] ; xx [ 149 ] = xx [ 28 ] + xx [ 122 ] ; xx [ 122 ] = xx [
125 ] / xx [ 60 ] ; xx [ 150 ] = xx [ 149 ] - xx [ 125 ] * xx [ 122 ] ; xx [
151 ] = xx [ 98 ] * xx [ 98 ] ; xx [ 152 ] = xx [ 2 ] * xx [ 151 ] - xx [ 0 ]
; xx [ 153 ] = xx [ 122 ] * xx [ 132 ] ; xx [ 154 ] = xx [ 123 ] - xx [ 153 ]
; xx [ 155 ] = xx [ 148 ] * xx [ 150 ] + xx [ 152 ] * xx [ 154 ] ; xx [ 156 ]
= xx [ 34 ] - xx [ 153 ] ; xx [ 153 ] = xx [ 132 ] / xx [ 60 ] ; xx [ 157 ] =
xx [ 127 ] - xx [ 153 ] * xx [ 132 ] ; xx [ 158 ] = xx [ 148 ] * xx [ 156 ] +
xx [ 157 ] * xx [ 152 ] ; xx [ 159 ] = xx [ 155 ] * xx [ 148 ] + xx [ 158 ] *
xx [ 152 ] ; xx [ 160 ] = xx [ 28 ] + xx [ 159 ] ; xx [ 161 ] = ( xx [ 151 ]
+ xx [ 92 ] * xx [ 92 ] ) * xx [ 2 ] - xx [ 0 ] ; xx [ 151 ] = xx [ 126 ] -
xx [ 153 ] * xx [ 40 ] ; xx [ 162 ] = xx [ 112 ] + xx [ 122 ] * xx [ 40 ] ;
xx [ 163 ] = xx [ 161 ] * ( xx [ 152 ] * xx [ 151 ] - xx [ 162 ] * xx [ 148 ]
) ; xx [ 164 ] = xx [ 150 ] * xx [ 152 ] - xx [ 148 ] * xx [ 154 ] ; xx [ 150
] = xx [ 152 ] * xx [ 156 ] - xx [ 148 ] * xx [ 157 ] ; xx [ 154 ] = xx [ 148
] * xx [ 164 ] + xx [ 152 ] * xx [ 150 ] ; xx [ 156 ] = xx [ 159 ] * xx [ 131
] - xx [ 154 ] * xx [ 134 ] ; xx [ 157 ] = xx [ 163 ] - xx [ 156 ] ; xx [ 159
] = xx [ 160 ] * xx [ 6 ] - xx [ 157 ] ; xx [ 165 ] = xx [ 40 ] / xx [ 60 ] ;
xx [ 166 ] = ( xx [ 162 ] * xx [ 152 ] + xx [ 148 ] * xx [ 151 ] ) * xx [ 161
] ; xx [ 151 ] = xx [ 131 ] * xx [ 163 ] + xx [ 166 ] * xx [ 134 ] ; xx [ 162
] = xx [ 155 ] * xx [ 152 ] - xx [ 158 ] * xx [ 148 ] ; xx [ 155 ] = xx [ 152
] * xx [ 164 ] - xx [ 148 ] * xx [ 150 ] ; xx [ 150 ] = xx [ 131 ] * xx [ 162
] - xx [ 134 ] * xx [ 155 ] ; xx [ 158 ] = xx [ 25 ] + ( xx [ 129 ] - xx [
165 ] * xx [ 40 ] ) * xx [ 161 ] * xx [ 161 ] - xx [ 151 ] - xx [ 151 ] + xx
[ 156 ] * xx [ 131 ] - xx [ 150 ] * xx [ 134 ] ; xx [ 129 ] = xx [ 6 ] * xx [
157 ] - xx [ 158 ] ; xx [ 151 ] = xx [ 6 ] * xx [ 159 ] - xx [ 129 ] ; ii [ 0
] = factorSymmetricPosDef ( xx + 151 , 1 , xx + 156 ) ; if ( ii [ 0 ] != 0 )
{ return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint2' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 156 ] = ( xx [ 89 ] + xx [ 136 ] * xx [ 6 ] + xx [ 147
] ) / xx [ 151 ] ; xx [ 161 ] = xx [ 22 ] * xx [ 93 ] ; xx [ 163 ] = xx [ 43
] - xx [ 2 ] * xx [ 161 ] * xx [ 93 ] ; xx [ 164 ] = xx [ 136 ] - xx [ 156 ]
* xx [ 159 ] ; xx [ 136 ] = xx [ 166 ] + xx [ 150 ] ; xx [ 150 ] = xx [ 136 ]
+ xx [ 6 ] * xx [ 162 ] ; xx [ 166 ] = xx [ 130 ] - xx [ 150 ] * xx [ 156 ] ;
xx [ 130 ] = xx [ 93 ] * xx [ 166 ] ; xx [ 167 ] = xx [ 93 ] * xx [ 164 ] ;
xx [ 168 ] = xx [ 164 ] + xx [ 2 ] * ( xx [ 3 ] * xx [ 130 ] - xx [ 167 ] *
xx [ 93 ] ) ; xx [ 164 ] = xx [ 2 ] * xx [ 3 ] * xx [ 161 ] ; xx [ 161 ] = xx
[ 111 ] + xx [ 2 ] * ( xx [ 137 ] - xx [ 133 ] ) + xx [ 164 ] ; xx [ 111 ] =
xx [ 96 ] * xx [ 161 ] ; xx [ 133 ] = xx [ 135 ] - xx [ 141 ] + xx [ 163 ] ;
xx [ 135 ] = xx [ 133 ] * xx [ 96 ] ; xx [ 137 ] = xx [ 6 ] * xx [ 96 ] ; xx
[ 141 ] = xx [ 143 ] * xx [ 93 ] - xx [ 3 ] * xx [ 144 ] ; xx [ 169 ] = xx [
143 ] * xx [ 3 ] + xx [ 93 ] * xx [ 144 ] ; xx [ 143 ] = xx [ 96 ] * xx [ 141
] - xx [ 169 ] * xx [ 94 ] ; xx [ 144 ] = xx [ 22 ] * xx [ 143 ] ; xx [ 170 ]
= xx [ 2 ] * ( xx [ 111 ] * xx [ 96 ] - xx [ 94 ] * xx [ 135 ] ) - xx [ 161 ]
- xx [ 2 ] * xx [ 94 ] * xx [ 137 ] - xx [ 2 ] * ( xx [ 94 ] * xx [ 141 ] +
xx [ 169 ] * xx [ 96 ] ) * xx [ 144 ] ; xx [ 141 ] = xx [ 3 ] * xx [ 93 ] ;
xx [ 161 ] = xx [ 2 ] * xx [ 141 ] ; xx [ 169 ] = xx [ 28 ] + xx [ 155 ] ; xx
[ 155 ] = xx [ 150 ] / xx [ 151 ] ; xx [ 171 ] = xx [ 169 ] - xx [ 150 ] * xx
[ 155 ] ; xx [ 172 ] = xx [ 3 ] * xx [ 3 ] ; xx [ 173 ] = xx [ 2 ] * xx [ 172
] - xx [ 0 ] ; xx [ 174 ] = xx [ 155 ] * xx [ 159 ] ; xx [ 175 ] = xx [ 162 ]
- xx [ 174 ] ; xx [ 176 ] = xx [ 161 ] * xx [ 171 ] + xx [ 173 ] * xx [ 175 ]
; xx [ 177 ] = xx [ 154 ] - xx [ 174 ] ; xx [ 174 ] = xx [ 159 ] / xx [ 151 ]
; xx [ 178 ] = xx [ 160 ] - xx [ 174 ] * xx [ 159 ] ; xx [ 179 ] = xx [ 161 ]
* xx [ 177 ] + xx [ 178 ] * xx [ 173 ] ; xx [ 180 ] = xx [ 176 ] * xx [ 161 ]
+ xx [ 179 ] * xx [ 173 ] ; xx [ 181 ] = ( xx [ 172 ] + xx [ 93 ] * xx [ 93 ]
) * xx [ 2 ] - xx [ 0 ] ; xx [ 172 ] = xx [ 157 ] - xx [ 174 ] * xx [ 129 ] ;
xx [ 182 ] = xx [ 136 ] + xx [ 155 ] * xx [ 129 ] ; xx [ 183 ] = xx [ 181 ] *
( xx [ 173 ] * xx [ 172 ] - xx [ 182 ] * xx [ 161 ] ) ; xx [ 184 ] = xx [ 171
] * xx [ 173 ] - xx [ 161 ] * xx [ 175 ] ; xx [ 171 ] = xx [ 177 ] * xx [ 173
] - xx [ 161 ] * xx [ 178 ] ; xx [ 175 ] = xx [ 161 ] * xx [ 184 ] + xx [ 173
] * xx [ 171 ] ; xx [ 177 ] = xx [ 180 ] * xx [ 163 ] - xx [ 175 ] * xx [ 164
] ; xx [ 178 ] = xx [ 183 ] - xx [ 177 ] ; xx [ 185 ] = xx [ 129 ] / xx [ 151
] ; xx [ 186 ] = ( xx [ 182 ] * xx [ 173 ] + xx [ 161 ] * xx [ 172 ] ) * xx [
181 ] ; xx [ 172 ] = xx [ 163 ] * xx [ 183 ] + xx [ 186 ] * xx [ 164 ] ; xx [
182 ] = xx [ 163 ] * ( xx [ 176 ] * xx [ 173 ] - xx [ 179 ] * xx [ 161 ] ) -
xx [ 164 ] * ( xx [ 173 ] * xx [ 184 ] - xx [ 161 ] * xx [ 171 ] ) ; xx [ 171
] = xx [ 6 ] * ( ( xx [ 28 ] + xx [ 180 ] ) * xx [ 6 ] - xx [ 178 ] ) - ( xx
[ 6 ] * xx [ 178 ] - ( ( xx [ 158 ] - xx [ 185 ] * xx [ 129 ] ) * xx [ 181 ]
* xx [ 181 ] - xx [ 172 ] - xx [ 172 ] + xx [ 177 ] * xx [ 163 ] - xx [ 182 ]
* xx [ 164 ] ) ) + xx [ 25 ] ; ii [ 0 ] = factorSymmetricPosDef ( xx + 171 ,
1 , xx + 158 ) ; if ( ii [ 0 ] != 0 ) { return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint1' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 158 ] = ( xx [ 89 ] + xx [ 156 ] * xx [ 129 ] + xx [
163 ] * xx [ 168 ] - xx [ 164 ] * ( xx [ 166 ] - ( xx [ 3 ] * xx [ 167 ] + xx
[ 130 ] * xx [ 93 ] ) * xx [ 2 ] ) + xx [ 168 ] * xx [ 6 ] + xx [ 170 ] ) /
xx [ 171 ] ; xx [ 89 ] = xx [ 158 ] * xx [ 164 ] ; xx [ 130 ] = xx [ 89 ] *
xx [ 93 ] ; xx [ 166 ] = xx [ 6 ] * xx [ 158 ] + xx [ 163 ] * xx [ 158 ] ; xx
[ 167 ] = xx [ 166 ] * xx [ 93 ] ; xx [ 168 ] = xx [ 89 ] - xx [ 2 ] * ( xx [
130 ] * xx [ 93 ] + xx [ 3 ] * xx [ 167 ] ) ; xx [ 89 ] = ( xx [ 167 ] * xx [
93 ] - xx [ 3 ] * xx [ 130 ] ) * xx [ 2 ] - xx [ 166 ] ; xx [ 130 ] = xx [
156 ] + xx [ 158 ] * xx [ 185 ] + xx [ 155 ] * xx [ 168 ] + xx [ 174 ] * xx [
89 ] ; xx [ 156 ] = xx [ 158 ] + xx [ 130 ] ; xx [ 166 ] = xx [ 168 ] + xx [
156 ] * xx [ 134 ] ; xx [ 167 ] = xx [ 89 ] - xx [ 130 ] * xx [ 6 ] - xx [
156 ] * xx [ 131 ] ; xx [ 89 ] = xx [ 92 ] * xx [ 167 ] ; xx [ 168 ] = xx [
92 ] * xx [ 166 ] ; xx [ 172 ] = xx [ 166 ] + xx [ 2 ] * ( xx [ 98 ] * xx [
89 ] - xx [ 168 ] * xx [ 92 ] ) ; xx [ 166 ] = xx [ 167 ] - ( xx [ 98 ] * xx
[ 168 ] + xx [ 89 ] * xx [ 92 ] ) * xx [ 2 ] ; xx [ 89 ] = xx [ 128 ] + xx [
156 ] * xx [ 165 ] + xx [ 172 ] * xx [ 122 ] + xx [ 153 ] * xx [ 166 ] ; xx [
128 ] = xx [ 156 ] + xx [ 89 ] ; xx [ 156 ] = xx [ 166 ] - xx [ 89 ] * xx [ 6
] - xx [ 128 ] * xx [ 90 ] ; xx [ 166 ] = xx [ 88 ] - ( xx [ 2 ] * xx [ 86 ]
* xx [ 87 ] + ( xx [ 52 ] + xx [ 50 ] ) * xx [ 2 ] ) + xx [ 22 ] ; xx [ 50 ]
= xx [ 2 ] * xx [ 24 ] * xx [ 23 ] ; xx [ 24 ] = xx [ 16 ] + xx [ 2 ] * ( xx
[ 20 ] * xx [ 19 ] - xx [ 10 ] * xx [ 4 ] ) - xx [ 50 ] + xx [ 22 ] ; xx [ 4
] = xx [ 24 ] / xx [ 9 ] ; xx [ 20 ] = xx [ 1 ] * xx [ 4 ] ; xx [ 52 ] = xx [
20 ] * xx [ 7 ] ; xx [ 86 ] = xx [ 2 ] * xx [ 5 ] * xx [ 52 ] ; xx [ 87 ] =
xx [ 20 ] - xx [ 2 ] * xx [ 52 ] * xx [ 7 ] ; xx [ 20 ] = xx [ 25 ] * xx [ 4
] + xx [ 45 ] * xx [ 86 ] - xx [ 44 ] * xx [ 87 ] ; xx [ 52 ] = xx [ 53 ] - (
xx [ 50 ] + ( xx [ 58 ] * xx [ 15 ] + xx [ 18 ] * xx [ 57 ] ) * xx [ 2 ] ) +
xx [ 22 ] ; xx [ 50 ] = ( xx [ 52 ] - ( xx [ 20 ] + xx [ 6 ] * xx [ 87 ] ) )
/ xx [ 46 ] ; xx [ 53 ] = xx [ 86 ] + xx [ 37 ] * xx [ 50 ] ; xx [ 57 ] = xx
[ 87 ] + xx [ 49 ] * xx [ 50 ] ; xx [ 58 ] = xx [ 57 ] * xx [ 11 ] ; xx [ 86
] = xx [ 53 ] * xx [ 11 ] ; xx [ 87 ] = xx [ 53 ] + xx [ 2 ] * ( xx [ 14 ] *
xx [ 58 ] - xx [ 86 ] * xx [ 11 ] ) ; xx [ 53 ] = xx [ 57 ] - ( xx [ 14 ] *
xx [ 86 ] + xx [ 58 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 57 ] = ( xx [ 166 ] - (
xx [ 20 ] + xx [ 41 ] * xx [ 50 ] + xx [ 87 ] * xx [ 82 ] - xx [ 74 ] * xx [
53 ] + xx [ 6 ] * xx [ 53 ] ) ) / xx [ 73 ] ; xx [ 20 ] = xx [ 53 ] + xx [ 81
] * xx [ 57 ] ; xx [ 53 ] = ( xx [ 20 ] - ( xx [ 17 ] * xx [ 17 ] * xx [ 20 ]
- ( xx [ 87 ] + xx [ 61 ] * xx [ 57 ] ) * xx [ 17 ] * xx [ 13 ] ) * xx [ 2 ]
) / xx [ 71 ] ; xx [ 20 ] = xx [ 17 ] * xx [ 53 ] ; xx [ 58 ] = xx [ 2 ] * xx
[ 17 ] * xx [ 20 ] - xx [ 53 ] ; xx [ 86 ] = xx [ 2 ] * xx [ 20 ] * xx [ 13 ]
; xx [ 20 ] = xx [ 57 ] - ( xx [ 65 ] * xx [ 58 ] - xx [ 77 ] * xx [ 86 ] ) ;
xx [ 57 ] = xx [ 82 ] * xx [ 20 ] - xx [ 86 ] ; xx [ 86 ] = xx [ 58 ] + xx [
6 ] * xx [ 20 ] - xx [ 74 ] * xx [ 20 ] ; xx [ 58 ] = xx [ 11 ] * xx [ 86 ] ;
xx [ 87 ] = xx [ 57 ] * xx [ 11 ] ; xx [ 88 ] = xx [ 57 ] - ( xx [ 14 ] * xx
[ 58 ] + xx [ 87 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 57 ] = xx [ 86 ] + xx [ 2
] * ( xx [ 14 ] * xx [ 87 ] - xx [ 58 ] * xx [ 11 ] ) ; xx [ 58 ] = xx [ 50 ]
- ( xx [ 55 ] * xx [ 20 ] + xx [ 68 ] * xx [ 88 ] + xx [ 57 ] * xx [ 62 ] ) ;
xx [ 50 ] = xx [ 20 ] + xx [ 58 ] ; xx [ 86 ] = xx [ 57 ] + xx [ 6 ] * xx [
58 ] - xx [ 44 ] * xx [ 50 ] ; xx [ 57 ] = xx [ 4 ] - ( xx [ 50 ] * xx [ 78 ]
+ ( xx [ 86 ] + xx [ 2 ] * ( xx [ 5 ] * ( xx [ 88 ] + xx [ 50 ] * xx [ 45 ] )
* xx [ 7 ] - xx [ 7 ] * xx [ 86 ] * xx [ 7 ] ) ) * xx [ 42 ] ) ; xx [ 4 ] =
xx [ 2 ] * xx [ 107 ] * xx [ 106 ] ; xx [ 50 ] = xx [ 100 ] + xx [ 2 ] * ( xx
[ 104 ] * xx [ 102 ] - xx [ 91 ] * xx [ 101 ] ) - xx [ 4 ] + xx [ 22 ] ; xx [
86 ] = xx [ 50 ] / xx [ 9 ] ; xx [ 87 ] = xx [ 1 ] * xx [ 86 ] ; xx [ 88 ] =
xx [ 87 ] * xx [ 85 ] ; xx [ 101 ] = xx [ 2 ] * xx [ 88 ] * xx [ 85 ] - xx [
87 ] ; xx [ 87 ] = xx [ 2 ] * xx [ 103 ] * xx [ 88 ] ; xx [ 88 ] = xx [ 90 ]
* xx [ 101 ] - xx [ 110 ] * xx [ 87 ] - xx [ 25 ] * xx [ 86 ] ; xx [ 104 ] =
xx [ 114 ] - xx [ 2 ] * ( xx [ 120 ] * xx [ 95 ] + xx [ 97 ] * xx [ 116 ] ) -
xx [ 4 ] + xx [ 22 ] ; xx [ 4 ] = ( xx [ 88 ] + xx [ 6 ] * xx [ 101 ] + xx [
104 ] ) / xx [ 60 ] ; xx [ 107 ] = xx [ 101 ] - xx [ 4 ] * xx [ 132 ] ; xx [
101 ] = xx [ 87 ] - xx [ 125 ] * xx [ 4 ] ; xx [ 87 ] = xx [ 92 ] * xx [ 101
] ; xx [ 114 ] = xx [ 92 ] * xx [ 107 ] ; xx [ 116 ] = xx [ 107 ] + xx [ 2 ]
* ( xx [ 98 ] * xx [ 87 ] - xx [ 114 ] * xx [ 92 ] ) ; xx [ 107 ] = xx [ 101
] - ( xx [ 98 ] * xx [ 114 ] + xx [ 87 ] * xx [ 92 ] ) * xx [ 2 ] ; xx [ 87 ]
= xx [ 88 ] + xx [ 4 ] * xx [ 40 ] + xx [ 131 ] * xx [ 116 ] - xx [ 134 ] *
xx [ 107 ] ; xx [ 88 ] = xx [ 119 ] + xx [ 2 ] * ( xx [ 94 ] * xx [ 142 ] -
xx [ 138 ] * xx [ 96 ] ) - xx [ 2 ] * xx [ 145 ] * xx [ 146 ] + xx [ 22 ] ;
xx [ 101 ] = ( xx [ 87 ] + xx [ 116 ] * xx [ 6 ] + xx [ 88 ] ) / xx [ 151 ] ;
xx [ 114 ] = xx [ 116 ] - xx [ 101 ] * xx [ 159 ] ; xx [ 116 ] = xx [ 107 ] -
xx [ 150 ] * xx [ 101 ] ; xx [ 107 ] = xx [ 93 ] * xx [ 116 ] ; xx [ 119 ] =
xx [ 93 ] * xx [ 114 ] ; xx [ 120 ] = xx [ 114 ] + xx [ 2 ] * ( xx [ 3 ] * xx
[ 107 ] - xx [ 119 ] * xx [ 93 ] ) ; xx [ 114 ] = xx [ 133 ] - ( ( xx [ 94 ]
* xx [ 111 ] + xx [ 135 ] * xx [ 96 ] ) * xx [ 2 ] + xx [ 2 ] * xx [ 137 ] *
xx [ 96 ] ) - xx [ 2 ] * xx [ 144 ] * xx [ 143 ] + xx [ 43 ] ; xx [ 43 ] = (
xx [ 87 ] + xx [ 101 ] * xx [ 129 ] + xx [ 163 ] * xx [ 120 ] - xx [ 164 ] *
( xx [ 116 ] - ( xx [ 3 ] * xx [ 119 ] + xx [ 107 ] * xx [ 93 ] ) * xx [ 2 ]
) + xx [ 120 ] * xx [ 6 ] + xx [ 114 ] ) / xx [ 171 ] ; xx [ 87 ] = xx [ 43 ]
* xx [ 164 ] ; xx [ 107 ] = xx [ 87 ] * xx [ 93 ] ; xx [ 111 ] = xx [ 6 ] *
xx [ 43 ] + xx [ 163 ] * xx [ 43 ] ; xx [ 116 ] = xx [ 111 ] * xx [ 93 ] ; xx
[ 119 ] = xx [ 87 ] - xx [ 2 ] * ( xx [ 107 ] * xx [ 93 ] + xx [ 3 ] * xx [
116 ] ) ; xx [ 87 ] = ( xx [ 116 ] * xx [ 93 ] - xx [ 3 ] * xx [ 107 ] ) * xx
[ 2 ] - xx [ 111 ] ; xx [ 107 ] = xx [ 101 ] + xx [ 43 ] * xx [ 185 ] + xx [
155 ] * xx [ 119 ] + xx [ 174 ] * xx [ 87 ] ; xx [ 101 ] = xx [ 43 ] + xx [
107 ] ; xx [ 111 ] = xx [ 119 ] + xx [ 101 ] * xx [ 134 ] ; xx [ 116 ] = xx [
87 ] - xx [ 107 ] * xx [ 6 ] - xx [ 101 ] * xx [ 131 ] ; xx [ 87 ] = xx [ 92
] * xx [ 116 ] ; xx [ 119 ] = xx [ 92 ] * xx [ 111 ] ; xx [ 120 ] = xx [ 111
] + xx [ 2 ] * ( xx [ 98 ] * xx [ 87 ] - xx [ 119 ] * xx [ 92 ] ) ; xx [ 111
] = xx [ 116 ] - ( xx [ 98 ] * xx [ 119 ] + xx [ 87 ] * xx [ 92 ] ) * xx [ 2
] ; xx [ 87 ] = xx [ 4 ] + xx [ 101 ] * xx [ 165 ] + xx [ 120 ] * xx [ 122 ]
+ xx [ 153 ] * xx [ 111 ] ; xx [ 4 ] = xx [ 101 ] + xx [ 87 ] ; xx [ 101 ] =
xx [ 111 ] - xx [ 87 ] * xx [ 6 ] - xx [ 4 ] * xx [ 90 ] ; xx [ 111 ] = xx [
86 ] + xx [ 42 ] * ( xx [ 101 ] - ( xx [ 103 ] * xx [ 85 ] * ( xx [ 120 ] +
xx [ 4 ] * xx [ 110 ] ) + xx [ 85 ] * xx [ 101 ] * xx [ 85 ] ) * xx [ 2 ] ) -
xx [ 4 ] * xx [ 78 ] ; xx [ 4 ] = xx [ 56 ] * xx [ 20 ] - xx [ 53 ] + xx [ 59
] * xx [ 58 ] + xx [ 57 ] * xx [ 26 ] + xx [ 43 ] * xx [ 170 ] + xx [ 107 ] *
xx [ 147 ] - xx [ 87 ] * xx [ 121 ] + xx [ 111 ] * xx [ 109 ] ; xx [ 142 ] =
xx [ 69 ] + xx [ 56 ] * xx [ 29 ] + xx [ 59 ] * xx [ 75 ] + ( xx [ 27 ] - (
xx [ 48 ] * xx [ 78 ] + ( xx [ 80 ] + xx [ 2 ] * ( xx [ 5 ] * ( xx [ 83 ] +
xx [ 48 ] * xx [ 45 ] ) * xx [ 7 ] - xx [ 7 ] * xx [ 80 ] * xx [ 7 ] ) ) * xx
[ 42 ] ) ) * xx [ 26 ] + xx [ 158 ] * xx [ 170 ] + xx [ 130 ] * xx [ 147 ] -
xx [ 89 ] * xx [ 121 ] + ( xx [ 99 ] + xx [ 42 ] * ( xx [ 156 ] - ( xx [ 103
] * xx [ 85 ] * ( xx [ 172 ] + xx [ 128 ] * xx [ 110 ] ) + xx [ 85 ] * xx [
156 ] * xx [ 85 ] ) * xx [ 2 ] ) - xx [ 128 ] * xx [ 78 ] ) * xx [ 109 ] ; xx
[ 143 ] = xx [ 4 ] ; xx [ 144 ] = xx [ 4 ] ; xx [ 145 ] = xx [ 166 ] * xx [
20 ] + xx [ 52 ] * xx [ 58 ] + xx [ 24 ] * xx [ 57 ] + xx [ 114 ] * xx [ 43 ]
+ xx [ 107 ] * xx [ 88 ] + xx [ 87 ] * xx [ 104 ] + xx [ 111 ] * xx [ 50 ] ;
xx [ 4 ] = state [ 3 ] + state [ 5 ] ; xx [ 20 ] = xx [ 4 ] * xx [ 4 ] * xx [
45 ] ; xx [ 27 ] = xx [ 20 ] * xx [ 7 ] ; xx [ 29 ] = xx [ 4 ] * xx [ 44 ] *
xx [ 4 ] ; xx [ 43 ] = xx [ 29 ] * xx [ 7 ] ; xx [ 48 ] = xx [ 28 ] * ( xx [
2 ] * ( xx [ 27 ] * xx [ 7 ] - xx [ 5 ] * xx [ 43 ] ) - xx [ 20 ] ) - input [
26 ] ; xx [ 20 ] = 4.010704565915762e-6 ; xx [ 53 ] = 4.010704565915763e-7 ;
xx [ 57 ] = xx [ 4 ] + state [ 7 ] ; xx [ 58 ] = state [ 9 ] + state [ 11 ] ;
xx [ 69 ] = xx [ 58 ] + state [ 13 ] ; xx [ 75 ] = xx [ 69 ] + state [ 15 ] ;
xx [ 80 ] = xx [ 105 ] * xx [ 23 ] - xx [ 21 ] * xx [ 106 ] ; xx [ 83 ] = xx
[ 21 ] * xx [ 105 ] + xx [ 23 ] * xx [ 106 ] ; xx [ 86 ] = ( xx [ 57 ] + xx [
75 ] ) * xx [ 53 ] - ( state [ 16 ] + pm_math_canonicalAngle ( xx [ 2 ] *
atan2 ( sqrt ( xx [ 80 ] * xx [ 80 ] ) , fabs ( - xx [ 83 ] ) ) * ( ( xx [ 83
] * xx [ 80 ] ) < 0.0 ? - 1.0 : + 1.0 ) - state [ 16 ] ) ) * xx [ 20 ] ; xx [
80 ] = xx [ 86 ] - input [ 28 ] ; xx [ 83 ] = xx [ 20 ] * state [ 6 ] + xx [
53 ] * state [ 7 ] + xx [ 80 ] + xx [ 6 ] * xx [ 48 ] ; xx [ 87 ] = xx [ 83 ]
/ xx [ 9 ] ; xx [ 89 ] = xx [ 48 ] - xx [ 1 ] * xx [ 87 ] ; xx [ 99 ] = ( xx
[ 6 ] * ( xx [ 4 ] + xx [ 57 ] ) * state [ 7 ] + ( xx [ 5 ] * xx [ 27 ] + xx
[ 43 ] * xx [ 7 ] ) * xx [ 2 ] - xx [ 29 ] ) * xx [ 28 ] - input [ 24 ] ; xx
[ 27 ] = xx [ 7 ] * xx [ 99 ] ; xx [ 29 ] = xx [ 5 ] * xx [ 27 ] ; xx [ 43 ]
= xx [ 7 ] * xx [ 89 ] ; xx [ 57 ] = xx [ 89 ] - ( xx [ 29 ] + xx [ 43 ] * xx
[ 7 ] ) * xx [ 2 ] ; xx [ 89 ] = xx [ 82 ] * state [ 3 ] * state [ 3 ] ; xx [
101 ] = xx [ 89 ] * xx [ 11 ] ; xx [ 107 ] = xx [ 74 ] * state [ 3 ] * state
[ 3 ] ; xx [ 111 ] = xx [ 107 ] * xx [ 11 ] ; xx [ 116 ] = xx [ 6 ] * ( state
[ 3 ] + xx [ 4 ] ) * state [ 5 ] + ( xx [ 14 ] * xx [ 101 ] + xx [ 111 ] * xx
[ 11 ] ) * xx [ 2 ] - xx [ 107 ] ; xx [ 107 ] = xx [ 2 ] * ( xx [ 101 ] * xx
[ 11 ] - xx [ 14 ] * xx [ 111 ] ) - xx [ 89 ] ; xx [ 89 ] = xx [ 116 ] * xx [
47 ] + xx [ 39 ] * xx [ 107 ] ; xx [ 39 ] = xx [ 57 ] - input [ 32 ] + xx [
89 ] ; xx [ 47 ] = xx [ 20 ] * state [ 4 ] + xx [ 53 ] * state [ 5 ] ; xx [
101 ] = xx [ 27 ] * xx [ 7 ] ; xx [ 27 ] = xx [ 99 ] + xx [ 2 ] * ( xx [ 5 ]
* xx [ 43 ] - xx [ 101 ] ) ; xx [ 43 ] = xx [ 116 ] * xx [ 30 ] + xx [ 107 ]
* xx [ 38 ] ; xx [ 30 ] = xx [ 80 ] - xx [ 25 ] * xx [ 87 ] + xx [ 27 ] * xx
[ 45 ] - xx [ 44 ] * xx [ 57 ] - input [ 34 ] - xx [ 43 ] ; xx [ 38 ] = ( xx
[ 47 ] + xx [ 30 ] + xx [ 39 ] * xx [ 6 ] ) / xx [ 46 ] ; xx [ 57 ] = xx [ 39
] - xx [ 49 ] * xx [ 38 ] ; xx [ 39 ] = xx [ 76 ] * xx [ 116 ] + xx [ 107 ] *
xx [ 35 ] ; xx [ 35 ] = xx [ 27 ] - input [ 30 ] + xx [ 39 ] - xx [ 37 ] * xx
[ 38 ] ; xx [ 27 ] = xx [ 11 ] * xx [ 35 ] ; xx [ 76 ] = xx [ 11 ] * xx [ 57
] ; xx [ 111 ] = xx [ 57 ] - ( xx [ 14 ] * xx [ 27 ] + xx [ 76 ] * xx [ 11 ]
) * xx [ 2 ] ; xx [ 57 ] = state [ 3 ] * state [ 3 ] ; xx [ 119 ] = xx [ 6 ]
* xx [ 57 ] ; xx [ 120 ] = xx [ 119 ] * xx [ 33 ] ; xx [ 33 ] = xx [ 111 ] -
input [ 38 ] + xx [ 120 ] ; xx [ 128 ] = xx [ 20 ] * state [ 2 ] + xx [ 53 ]
* state [ 3 ] ; xx [ 130 ] = xx [ 35 ] + xx [ 2 ] * ( xx [ 14 ] * xx [ 76 ] -
xx [ 27 ] * xx [ 11 ] ) ; xx [ 27 ] = xx [ 119 ] * xx [ 63 ] ; xx [ 35 ] = (
xx [ 128 ] + xx [ 30 ] - xx [ 41 ] * xx [ 38 ] + xx [ 130 ] * xx [ 82 ] - xx
[ 74 ] * xx [ 111 ] - input [ 40 ] - xx [ 27 ] + xx [ 33 ] * xx [ 6 ] ) / xx
[ 73 ] ; xx [ 30 ] = xx [ 33 ] - xx [ 81 ] * xx [ 35 ] ; xx [ 33 ] = xx [ 79
] * xx [ 119 ] ; xx [ 63 ] = ( xx [ 30 ] - ( xx [ 17 ] * xx [ 17 ] * xx [ 30
] - xx [ 17 ] * ( xx [ 130 ] - input [ 36 ] + xx [ 33 ] - xx [ 61 ] * xx [ 35
] ) * xx [ 13 ] ) * xx [ 2 ] ) / xx [ 71 ] ; xx [ 30 ] = 1.0 ; xx [ 76 ] = xx
[ 64 ] * state [ 3 ] * state [ 3 ] ; xx [ 64 ] = xx [ 30 ] * xx [ 76 ] ; xx [
79 ] = 2.220446049250313e-16 ; xx [ 111 ] = xx [ 2 ] * xx [ 72 ] * state [ 3
] * state [ 3 ] ; xx [ 72 ] = xx [ 79 ] * xx [ 111 ] ; xx [ 130 ] = xx [ 64 ]
- xx [ 72 ] ; xx [ 133 ] = xx [ 30 ] * xx [ 111 ] ; xx [ 30 ] = xx [ 79 ] *
xx [ 76 ] ; xx [ 76 ] = xx [ 133 ] + xx [ 30 ] ; xx [ 79 ] = xx [ 2 ] * xx [
18 ] * xx [ 15 ] ; xx [ 111 ] = xx [ 2 ] * xx [ 66 ] * state [ 5 ] * state [
5 ] ; xx [ 66 ] = xx [ 70 ] * state [ 5 ] * state [ 5 ] ; xx [ 135 ] = xx [ 2
] * xx [ 15 ] * xx [ 15 ] - xx [ 0 ] ; xx [ 137 ] = xx [ 79 ] * xx [ 111 ] -
xx [ 66 ] * xx [ 135 ] ; xx [ 138 ] = xx [ 2 ] * xx [ 19 ] * xx [ 19 ] - xx [
0 ] ; xx [ 146 ] = xx [ 2 ] * xx [ 138 ] * state [ 3 ] * state [ 5 ] ; xx [
156 ] = xx [ 130 ] * xx [ 67 ] - xx [ 76 ] * xx [ 70 ] + xx [ 137 ] - xx [
146 ] ; xx [ 158 ] = 4.0 ; xx [ 167 ] = xx [ 10 ] * xx [ 19 ] ; xx [ 168 ] =
xx [ 158 ] * xx [ 167 ] * state [ 3 ] * state [ 5 ] ; xx [ 172 ] = xx [ 111 ]
* xx [ 135 ] ; xx [ 111 ] = xx [ 79 ] * xx [ 66 ] ; xx [ 66 ] = xx [ 168 ] -
( xx [ 67 ] * xx [ 76 ] + xx [ 130 ] * xx [ 70 ] + xx [ 172 ] + xx [ 111 ] )
; xx [ 79 ] = 1.0e-3 ; xx [ 135 ] = xx [ 17 ] * xx [ 79 ] ; xx [ 176 ] = xx [
57 ] * ( xx [ 79 ] - xx [ 2 ] * xx [ 17 ] * xx [ 135 ] ) ; xx [ 177 ] = xx [
12 ] * xx [ 12 ] * xx [ 176 ] ; xx [ 178 ] = xx [ 2 ] * xx [ 135 ] * xx [ 13
] * xx [ 57 ] ; xx [ 57 ] = xx [ 12 ] * xx [ 12 ] * xx [ 178 ] ; xx [ 12 ] =
state [ 5 ] * state [ 5 ] ; xx [ 135 ] = xx [ 79 ] * xx [ 11 ] ; xx [ 179 ] =
xx [ 12 ] * ( xx [ 79 ] - xx [ 2 ] * xx [ 135 ] * xx [ 11 ] ) ; xx [ 180 ] =
xx [ 2 ] * xx [ 14 ] * xx [ 135 ] * xx [ 12 ] ; xx [ 12 ] = xx [ 18 ] * xx [
180 ] ; xx [ 135 ] = xx [ 18 ] * xx [ 179 ] ; xx [ 181 ] = ( xx [ 6 ] - xx [
51 ] ) * state [ 5 ] * state [ 3 ] ; xx [ 51 ] = xx [ 2 ] * xx [ 54 ] * state
[ 5 ] * state [ 3 ] ; xx [ 54 ] = xx [ 18 ] * xx [ 51 ] ; xx [ 183 ] = xx [
18 ] * xx [ 181 ] ; xx [ 184 ] = state [ 7 ] * state [ 7 ] ; xx [ 187 ] = xx
[ 79 ] * xx [ 7 ] ; xx [ 188 ] = xx [ 184 ] * ( xx [ 79 ] - xx [ 2 ] * xx [
187 ] * xx [ 7 ] ) ; xx [ 189 ] = xx [ 2 ] * xx [ 5 ] * xx [ 187 ] * xx [ 184
] ; xx [ 184 ] = xx [ 10 ] * xx [ 189 ] ; xx [ 187 ] = xx [ 10 ] * xx [ 188 ]
; xx [ 190 ] = xx [ 4 ] * xx [ 16 ] * state [ 7 ] ; xx [ 16 ] = xx [ 4 ] * xx
[ 2 ] * xx [ 8 ] * state [ 7 ] ; xx [ 8 ] = xx [ 10 ] * xx [ 16 ] ; xx [ 191
] = xx [ 10 ] * xx [ 190 ] ; xx [ 192 ] = xx [ 36 ] * state [ 7 ] * state [ 7
] ; xx [ 193 ] = xx [ 2 ] * xx [ 167 ] ; xx [ 167 ] = xx [ 2 ] * xx [ 31 ] *
state [ 7 ] * state [ 7 ] ; xx [ 31 ] = xx [ 94 ] * xx [ 96 ] ; xx [ 194 ] =
xx [ 2 ] * xx [ 31 ] * state [ 9 ] * state [ 9 ] ; xx [ 195 ] = xx [ 2 ] * xx
[ 94 ] * xx [ 94 ] - xx [ 0 ] ; xx [ 196 ] = xx [ 195 ] * state [ 9 ] * state
[ 9 ] ; xx [ 197 ] = xx [ 2 ] * xx [ 31 ] ; xx [ 31 ] = xx [ 2 ] * xx [ 141 ]
* state [ 11 ] * state [ 11 ] ; xx [ 141 ] = xx [ 173 ] * state [ 11 ] *
state [ 11 ] ; xx [ 198 ] = xx [ 2 ] * xx [ 95 ] * xx [ 95 ] - xx [ 0 ] ; xx
[ 199 ] = xx [ 161 ] * xx [ 194 ] - xx [ 196 ] * xx [ 173 ] + xx [ 197 ] * xx
[ 31 ] - xx [ 141 ] * xx [ 195 ] - xx [ 2 ] * xx [ 198 ] * state [ 9 ] *
state [ 11 ] ; xx [ 200 ] = xx [ 199 ] * xx [ 152 ] ; xx [ 201 ] = xx [ 161 ]
* xx [ 196 ] ; xx [ 161 ] = xx [ 194 ] * xx [ 173 ] ; xx [ 173 ] = xx [ 31 ]
* xx [ 195 ] ; xx [ 31 ] = xx [ 197 ] * xx [ 141 ] ; xx [ 141 ] = xx [ 97 ] *
xx [ 95 ] ; xx [ 195 ] = xx [ 158 ] * xx [ 141 ] * state [ 9 ] * state [ 11 ]
; xx [ 158 ] = xx [ 201 ] + xx [ 161 ] + xx [ 173 ] + xx [ 31 ] + xx [ 195 ]
; xx [ 197 ] = xx [ 2 ] * xx [ 141 ] ; xx [ 141 ] = xx [ 2 ] * xx [ 140 ] *
state [ 13 ] * state [ 13 ] ; xx [ 140 ] = xx [ 152 ] * state [ 13 ] * state
[ 13 ] ; xx [ 202 ] = xx [ 197 ] * xx [ 141 ] - xx [ 140 ] * xx [ 198 ] ; xx
[ 203 ] = xx [ 2 ] * xx [ 102 ] * xx [ 102 ] - xx [ 0 ] ; xx [ 204 ] = xx [ 2
] * xx [ 58 ] * xx [ 203 ] * state [ 13 ] ; xx [ 205 ] = xx [ 200 ] + xx [
158 ] * xx [ 148 ] + xx [ 202 ] - xx [ 204 ] ; xx [ 206 ] = xx [ 148 ] * xx [
199 ] ; xx [ 207 ] = xx [ 141 ] * xx [ 198 ] ; xx [ 141 ] = xx [ 197 ] * xx [
140 ] ; xx [ 140 ] = xx [ 2 ] * xx [ 91 ] * xx [ 102 ] ; xx [ 197 ] = xx [ 2
] * xx [ 58 ] * xx [ 140 ] * state [ 13 ] ; xx [ 198 ] = xx [ 158 ] * xx [
152 ] - xx [ 206 ] + xx [ 207 ] + xx [ 141 ] - xx [ 197 ] ; xx [ 208 ] = xx [
124 ] * state [ 15 ] * state [ 15 ] ; xx [ 209 ] = xx [ 2 ] * xx [ 108 ] *
state [ 15 ] * state [ 15 ] ; xx [ 108 ] = state [ 9 ] * state [ 9 ] ; xx [
210 ] = xx [ 79 ] * xx [ 96 ] ; xx [ 211 ] = state [ 11 ] * state [ 11 ] ; xx
[ 212 ] = xx [ 79 ] * xx [ 93 ] ; xx [ 213 ] = xx [ 211 ] * ( xx [ 2 ] * xx [
212 ] * xx [ 93 ] - xx [ 79 ] ) ; xx [ 214 ] = xx [ 213 ] * xx [ 96 ] ; xx [
215 ] = xx [ 2 ] * xx [ 3 ] * xx [ 212 ] * xx [ 211 ] ; xx [ 211 ] = xx [ 215
] * xx [ 96 ] ; xx [ 212 ] = xx [ 2 ] * xx [ 139 ] * state [ 11 ] * state [ 9
] ; xx [ 139 ] = xx [ 212 ] * xx [ 96 ] ; xx [ 216 ] = ( xx [ 6 ] - xx [ 115
] ) * state [ 11 ] * state [ 9 ] ; xx [ 115 ] = xx [ 216 ] * xx [ 96 ] ; xx [
217 ] = state [ 13 ] * state [ 13 ] ; xx [ 218 ] = xx [ 79 ] * xx [ 92 ] ; xx
[ 219 ] = xx [ 217 ] * ( xx [ 2 ] * xx [ 218 ] * xx [ 92 ] - xx [ 79 ] ) ; xx
[ 220 ] = xx [ 97 ] * xx [ 219 ] ; xx [ 221 ] = xx [ 2 ] * xx [ 98 ] * xx [
218 ] * xx [ 217 ] ; xx [ 217 ] = xx [ 97 ] * xx [ 221 ] ; xx [ 218 ] = xx [
58 ] * xx [ 2 ] * xx [ 117 ] * state [ 13 ] ; xx [ 117 ] = xx [ 97 ] * xx [
218 ] ; xx [ 222 ] = xx [ 58 ] * ( xx [ 6 ] - xx [ 113 ] ) * state [ 13 ] ;
xx [ 113 ] = xx [ 97 ] * xx [ 222 ] ; xx [ 223 ] = state [ 15 ] * state [ 15
] ; xx [ 224 ] = xx [ 79 ] * xx [ 85 ] ; xx [ 225 ] = xx [ 223 ] * ( xx [ 2 ]
* xx [ 224 ] * xx [ 85 ] - xx [ 79 ] ) ; xx [ 226 ] = xx [ 2 ] * xx [ 103 ] *
xx [ 224 ] * xx [ 223 ] ; xx [ 223 ] = xx [ 91 ] * xx [ 226 ] ; xx [ 224 ] =
xx [ 91 ] * xx [ 225 ] ; xx [ 227 ] = xx [ 69 ] * xx [ 100 ] * state [ 15 ] ;
xx [ 100 ] = xx [ 91 ] * xx [ 227 ] ; xx [ 228 ] = xx [ 69 ] * xx [ 2 ] * xx
[ 84 ] * state [ 15 ] ; xx [ 84 ] = xx [ 91 ] * xx [ 228 ] ; xx [ 229 ] = xx
[ 17 ] * xx [ 63 ] ; xx [ 230 ] = xx [ 2 ] * xx [ 17 ] * xx [ 229 ] - xx [ 63
] ; xx [ 231 ] = xx [ 2 ] * xx [ 229 ] * xx [ 13 ] ; xx [ 229 ] = xx [ 35 ] +
xx [ 65 ] * xx [ 230 ] - xx [ 77 ] * xx [ 231 ] ; xx [ 35 ] = xx [ 119 ] - xx
[ 231 ] - xx [ 229 ] * xx [ 82 ] ; xx [ 231 ] = xx [ 230 ] - xx [ 229 ] * xx
[ 6 ] + xx [ 229 ] * xx [ 74 ] ; xx [ 230 ] = xx [ 231 ] * xx [ 11 ] ; xx [
232 ] = xx [ 11 ] * xx [ 35 ] ; xx [ 233 ] = xx [ 35 ] - ( xx [ 14 ] * xx [
230 ] + xx [ 232 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 35 ] = xx [ 231 ] + xx [ 2
] * ( xx [ 14 ] * xx [ 232 ] - xx [ 230 ] * xx [ 11 ] ) ; xx [ 230 ] = xx [
38 ] + xx [ 68 ] * xx [ 233 ] + xx [ 35 ] * xx [ 62 ] - xx [ 229 ] * xx [ 55
] ; xx [ 38 ] = xx [ 229 ] + xx [ 230 ] ; xx [ 231 ] = xx [ 35 ] - xx [ 230 ]
* xx [ 6 ] + xx [ 107 ] + xx [ 38 ] * xx [ 44 ] ; xx [ 35 ] = xx [ 87 ] + (
xx [ 231 ] + xx [ 2 ] * ( xx [ 5 ] * xx [ 7 ] * ( xx [ 233 ] + xx [ 116 ] -
xx [ 38 ] * xx [ 45 ] ) - xx [ 231 ] * xx [ 7 ] * xx [ 7 ] ) ) * xx [ 42 ] -
xx [ 38 ] * xx [ 78 ] ; xx [ 38 ] = xx [ 20 ] * state [ 8 ] + xx [ 53 ] *
state [ 9 ] ; xx [ 87 ] = xx [ 69 ] * xx [ 69 ] * xx [ 90 ] ; xx [ 231 ] = xx
[ 87 ] * xx [ 85 ] ; xx [ 232 ] = xx [ 69 ] * xx [ 69 ] * xx [ 110 ] ; xx [
233 ] = xx [ 232 ] * xx [ 85 ] ; xx [ 234 ] = ( xx [ 2 ] * ( xx [ 103 ] * xx
[ 231 ] + xx [ 233 ] * xx [ 85 ] ) - xx [ 232 ] ) * xx [ 28 ] - input [ 20 ]
; xx [ 232 ] = input [ 22 ] + xx [ 86 ] ; xx [ 86 ] = xx [ 20 ] * state [ 14
] + xx [ 53 ] * state [ 15 ] + xx [ 232 ] + xx [ 6 ] * xx [ 234 ] ; xx [ 235
] = xx [ 86 ] / xx [ 9 ] ; xx [ 236 ] = xx [ 234 ] - xx [ 1 ] * xx [ 235 ] ;
xx [ 237 ] = xx [ 28 ] * ( ( xx [ 231 ] * xx [ 85 ] - xx [ 103 ] * xx [ 233 ]
) * xx [ 2 ] - xx [ 87 ] - xx [ 6 ] * ( xx [ 69 ] + xx [ 75 ] ) * state [ 15
] ) - input [ 18 ] ; xx [ 28 ] = xx [ 85 ] * xx [ 237 ] ; xx [ 75 ] = xx [
103 ] * xx [ 28 ] ; xx [ 87 ] = xx [ 85 ] * xx [ 236 ] ; xx [ 231 ] = xx [
236 ] + xx [ 2 ] * ( xx [ 75 ] - xx [ 87 ] * xx [ 85 ] ) ; xx [ 233 ] = xx [
58 ] * xx [ 131 ] * xx [ 58 ] ; xx [ 236 ] = xx [ 233 ] * xx [ 92 ] ; xx [
238 ] = xx [ 58 ] * xx [ 58 ] * xx [ 134 ] ; xx [ 239 ] = xx [ 238 ] * xx [
92 ] ; xx [ 240 ] = ( xx [ 236 ] * xx [ 92 ] - xx [ 98 ] * xx [ 239 ] ) * xx
[ 2 ] - xx [ 233 ] - xx [ 6 ] * ( xx [ 58 ] + xx [ 69 ] ) * state [ 13 ] ; xx
[ 233 ] = xx [ 2 ] * ( xx [ 98 ] * xx [ 236 ] + xx [ 239 ] * xx [ 92 ] ) - xx
[ 238 ] ; xx [ 236 ] = xx [ 34 ] * xx [ 240 ] + xx [ 127 ] * xx [ 233 ] ; xx
[ 34 ] = xx [ 231 ] - input [ 14 ] + xx [ 236 ] ; xx [ 127 ] = xx [ 20 ] *
state [ 12 ] + xx [ 53 ] * state [ 13 ] ; xx [ 238 ] = xx [ 233 ] * xx [ 126
] - xx [ 112 ] * xx [ 240 ] ; xx [ 112 ] = xx [ 28 ] * xx [ 85 ] ; xx [ 28 ]
= xx [ 237 ] - ( xx [ 103 ] * xx [ 87 ] + xx [ 112 ] ) * xx [ 2 ] ; xx [ 87 ]
= xx [ 238 ] - ( input [ 16 ] + xx [ 232 ] - xx [ 25 ] * xx [ 235 ] + xx [ 90
] * xx [ 231 ] - xx [ 110 ] * xx [ 28 ] ) ; xx [ 126 ] = ( xx [ 127 ] + xx [
34 ] * xx [ 6 ] - xx [ 87 ] ) / xx [ 60 ] ; xx [ 231 ] = xx [ 34 ] - xx [ 126
] * xx [ 132 ] ; xx [ 34 ] = xx [ 149 ] * xx [ 240 ] + xx [ 233 ] * xx [ 123
] ; xx [ 123 ] = xx [ 28 ] - input [ 12 ] + xx [ 34 ] - xx [ 125 ] * xx [ 126
] ; xx [ 28 ] = xx [ 92 ] * xx [ 123 ] ; xx [ 149 ] = xx [ 92 ] * xx [ 231 ]
; xx [ 239 ] = xx [ 231 ] + xx [ 2 ] * ( xx [ 98 ] * xx [ 28 ] - xx [ 149 ] *
xx [ 92 ] ) ; xx [ 231 ] = xx [ 163 ] * state [ 9 ] * state [ 9 ] ; xx [ 241
] = xx [ 231 ] * xx [ 93 ] ; xx [ 242 ] = xx [ 164 ] * state [ 9 ] * state [
9 ] ; xx [ 243 ] = xx [ 242 ] * xx [ 93 ] ; xx [ 244 ] = ( xx [ 241 ] * xx [
93 ] - xx [ 3 ] * xx [ 243 ] ) * xx [ 2 ] - xx [ 231 ] - xx [ 6 ] * ( state [
9 ] + xx [ 58 ] ) * state [ 11 ] ; xx [ 58 ] = xx [ 2 ] * ( xx [ 3 ] * xx [
241 ] + xx [ 243 ] * xx [ 93 ] ) - xx [ 242 ] ; xx [ 231 ] = xx [ 154 ] * xx
[ 244 ] + xx [ 160 ] * xx [ 58 ] ; xx [ 154 ] = xx [ 239 ] - input [ 8 ] + xx
[ 231 ] ; xx [ 160 ] = xx [ 20 ] * state [ 10 ] + xx [ 53 ] * state [ 11 ] ;
xx [ 20 ] = xx [ 123 ] - ( xx [ 98 ] * xx [ 149 ] + xx [ 28 ] * xx [ 92 ] ) *
xx [ 2 ] ; xx [ 28 ] = xx [ 58 ] * xx [ 157 ] - xx [ 136 ] * xx [ 244 ] ; xx
[ 53 ] = xx [ 87 ] - xx [ 126 ] * xx [ 40 ] - ( xx [ 131 ] * xx [ 239 ] - xx
[ 134 ] * xx [ 20 ] ) - input [ 10 ] + xx [ 28 ] ; xx [ 87 ] = ( xx [ 160 ] +
xx [ 154 ] * xx [ 6 ] - xx [ 53 ] ) / xx [ 151 ] ; xx [ 123 ] = xx [ 154 ] -
xx [ 87 ] * xx [ 159 ] ; xx [ 136 ] = xx [ 169 ] * xx [ 244 ] + xx [ 58 ] *
xx [ 162 ] ; xx [ 149 ] = xx [ 20 ] - input [ 6 ] + xx [ 136 ] - xx [ 150 ] *
xx [ 87 ] ; xx [ 20 ] = xx [ 93 ] * xx [ 149 ] ; xx [ 154 ] = xx [ 93 ] * xx
[ 123 ] ; xx [ 157 ] = xx [ 123 ] + xx [ 2 ] * ( xx [ 3 ] * xx [ 20 ] - xx [
154 ] * xx [ 93 ] ) ; xx [ 123 ] = xx [ 6 ] * xx [ 108 ] ; xx [ 162 ] = xx [
175 ] * xx [ 123 ] ; xx [ 169 ] = ( xx [ 186 ] + xx [ 182 ] ) * xx [ 123 ] ;
xx [ 175 ] = ( xx [ 38 ] + xx [ 6 ] * ( xx [ 157 ] - input [ 2 ] - xx [ 162 ]
) - ( xx [ 53 ] - xx [ 87 ] * xx [ 129 ] - ( xx [ 163 ] * xx [ 157 ] - xx [
164 ] * ( xx [ 149 ] - ( xx [ 3 ] * xx [ 154 ] + xx [ 20 ] * xx [ 93 ] ) * xx
[ 2 ] ) ) - input [ 4 ] + xx [ 169 ] ) ) / xx [ 171 ] ; xx [ 20 ] = xx [ 123
] - xx [ 175 ] * xx [ 164 ] ; xx [ 53 ] = xx [ 20 ] * xx [ 93 ] ; xx [ 149 ]
= xx [ 6 ] * xx [ 175 ] + xx [ 163 ] * xx [ 175 ] ; xx [ 154 ] = xx [ 149 ] *
xx [ 93 ] ; xx [ 157 ] = xx [ 2 ] * ( xx [ 53 ] * xx [ 93 ] - xx [ 3 ] * xx [
154 ] ) - xx [ 20 ] ; xx [ 20 ] = ( xx [ 3 ] * xx [ 53 ] + xx [ 154 ] * xx [
93 ] ) * xx [ 2 ] - xx [ 149 ] ; xx [ 53 ] = xx [ 87 ] + xx [ 175 ] * xx [
185 ] + xx [ 155 ] * xx [ 157 ] + xx [ 174 ] * xx [ 20 ] ; xx [ 87 ] = xx [
175 ] + xx [ 53 ] ; xx [ 149 ] = xx [ 157 ] + xx [ 244 ] + xx [ 87 ] * xx [
134 ] ; xx [ 154 ] = xx [ 20 ] - xx [ 53 ] * xx [ 6 ] + xx [ 58 ] - xx [ 87 ]
* xx [ 131 ] ; xx [ 20 ] = xx [ 92 ] * xx [ 154 ] ; xx [ 157 ] = xx [ 92 ] *
xx [ 149 ] ; xx [ 182 ] = xx [ 149 ] + xx [ 2 ] * ( xx [ 98 ] * xx [ 20 ] -
xx [ 157 ] * xx [ 92 ] ) ; xx [ 149 ] = xx [ 154 ] - ( xx [ 98 ] * xx [ 157 ]
+ xx [ 20 ] * xx [ 92 ] ) * xx [ 2 ] ; xx [ 20 ] = xx [ 126 ] + xx [ 87 ] *
xx [ 165 ] + xx [ 182 ] * xx [ 122 ] + xx [ 153 ] * xx [ 149 ] ; xx [ 126 ] =
xx [ 87 ] + xx [ 20 ] ; xx [ 87 ] = xx [ 149 ] - xx [ 20 ] * xx [ 6 ] + xx [
233 ] - xx [ 126 ] * xx [ 90 ] ; xx [ 149 ] = xx [ 235 ] + xx [ 42 ] * ( xx [
87 ] - ( xx [ 103 ] * xx [ 85 ] * ( xx [ 182 ] + xx [ 240 ] + xx [ 126 ] * xx
[ 110 ] ) + xx [ 85 ] * xx [ 87 ] * xx [ 85 ] ) * xx [ 2 ] ) - xx [ 126 ] *
xx [ 78 ] ; xx [ 87 ] = xx [ 64 ] - xx [ 72 ] ; xx [ 64 ] = xx [ 133 ] + xx [
30 ] ; xx [ 30 ] = xx [ 87 ] * xx [ 70 ] + xx [ 64 ] * xx [ 67 ] + xx [ 111 ]
+ xx [ 172 ] - xx [ 168 ] ; xx [ 72 ] = xx [ 67 ] * xx [ 87 ] - xx [ 64 ] *
xx [ 70 ] + xx [ 137 ] - xx [ 146 ] ; xx [ 67 ] = xx [ 161 ] + xx [ 201 ] +
xx [ 31 ] + xx [ 173 ] + xx [ 195 ] ; xx [ 31 ] = xx [ 206 ] - xx [ 67 ] * xx
[ 152 ] - ( xx [ 141 ] + xx [ 207 ] ) + xx [ 197 ] ; xx [ 70 ] = xx [ 67 ] *
xx [ 148 ] + xx [ 200 ] + xx [ 202 ] - xx [ 204 ] ; xx [ 172 ] = xx [ 63 ] -
( xx [ 44 ] * xx [ 156 ] + xx [ 66 ] * xx [ 45 ] + xx [ 176 ] - ( xx [ 177 ]
- xx [ 57 ] ) * xx [ 2 ] - ( xx [ 74 ] * xx [ 76 ] + xx [ 130 ] * xx [ 82 ] )
+ xx [ 179 ] - ( xx [ 12 ] * xx [ 15 ] + xx [ 18 ] * xx [ 135 ] ) * xx [ 2 ]
+ xx [ 2 ] * ( xx [ 181 ] - ( xx [ 54 ] * xx [ 15 ] + xx [ 18 ] * xx [ 183 ]
) * xx [ 2 ] ) + xx [ 188 ] + xx [ 2 ] * ( xx [ 184 ] * xx [ 19 ] - xx [ 10 ]
* xx [ 187 ] ) + ( xx [ 190 ] + xx [ 2 ] * ( xx [ 8 ] * xx [ 19 ] - xx [ 10 ]
* xx [ 191 ] ) ) * xx [ 2 ] - xx [ 22 ] * ( xx [ 156 ] * xx [ 36 ] - xx [ 66
] * xx [ 32 ] - ( xx [ 192 ] * xx [ 138 ] + xx [ 193 ] * xx [ 167 ] ) - xx [
2 ] * xx [ 4 ] * ( xx [ 2 ] * xx [ 21 ] * xx [ 21 ] - xx [ 0 ] ) * state [ 7
] ) - ( xx [ 22 ] * ( xx [ 205 ] * xx [ 124 ] + xx [ 118 ] * xx [ 198 ] - (
xx [ 208 ] * xx [ 203 ] + xx [ 140 ] * xx [ 209 ] ) - xx [ 2 ] * xx [ 69 ] *
( xx [ 2 ] * xx [ 105 ] * xx [ 105 ] - xx [ 0 ] ) * state [ 15 ] ) + xx [ 90
] * xx [ 205 ] + xx [ 110 ] * xx [ 198 ] + xx [ 131 ] * xx [ 199 ] + xx [ 158
] * xx [ 134 ] + xx [ 108 ] * ( xx [ 2 ] * xx [ 210 ] * xx [ 96 ] - xx [ 79 ]
) - ( xx [ 163 ] * xx [ 196 ] - xx [ 164 ] * xx [ 194 ] ) + xx [ 213 ] - ( xx
[ 214 ] * xx [ 96 ] - xx [ 94 ] * xx [ 211 ] ) * xx [ 2 ] + xx [ 2 ] * ( ( xx
[ 94 ] * xx [ 139 ] + xx [ 115 ] * xx [ 96 ] ) * xx [ 2 ] - xx [ 216 ] ) + xx
[ 219 ] - ( xx [ 97 ] * xx [ 220 ] - xx [ 217 ] * xx [ 95 ] ) * xx [ 2 ] + xx
[ 2 ] * ( ( xx [ 117 ] * xx [ 95 ] + xx [ 97 ] * xx [ 113 ] ) * xx [ 2 ] - xx
[ 222 ] ) + xx [ 225 ] - xx [ 2 ] * ( xx [ 223 ] * xx [ 102 ] + xx [ 91 ] *
xx [ 224 ] ) + xx [ 2 ] * ( xx [ 2 ] * ( xx [ 91 ] * xx [ 100 ] - xx [ 84 ] *
xx [ 102 ] ) - xx [ 227 ] ) ) ) + xx [ 229 ] * xx [ 56 ] + xx [ 230 ] * xx [
59 ] + xx [ 35 ] * xx [ 26 ] - xx [ 175 ] * xx [ 170 ] - xx [ 53 ] * xx [ 147
] + xx [ 20 ] * xx [ 121 ] - xx [ 149 ] * xx [ 109 ] ; xx [ 173 ] = xx [ 229
] * xx [ 166 ] - ( xx [ 44 ] * xx [ 30 ] + xx [ 45 ] * xx [ 72 ] + xx [ 74 ]
* xx [ 87 ] - xx [ 64 ] * xx [ 82 ] + xx [ 178 ] - xx [ 2 ] * ( xx [ 57 ] +
xx [ 177 ] ) + xx [ 2 ] * ( xx [ 18 ] * xx [ 12 ] - xx [ 135 ] * xx [ 15 ] )
- xx [ 180 ] + xx [ 2 ] * ( xx [ 2 ] * ( xx [ 18 ] * xx [ 54 ] - xx [ 183 ] *
xx [ 15 ] ) - xx [ 51 ] ) + ( xx [ 187 ] * xx [ 19 ] + xx [ 10 ] * xx [ 184 ]
) * xx [ 2 ] - xx [ 189 ] + xx [ 2 ] * ( ( xx [ 191 ] * xx [ 19 ] + xx [ 10 ]
* xx [ 8 ] ) * xx [ 2 ] - xx [ 16 ] ) - xx [ 22 ] * ( xx [ 30 ] * xx [ 36 ] -
xx [ 32 ] * xx [ 72 ] + xx [ 167 ] * xx [ 138 ] - xx [ 193 ] * xx [ 192 ] -
xx [ 2 ] * xx [ 4 ] * xx [ 2 ] * xx [ 21 ] * xx [ 23 ] * state [ 7 ] ) - ( xx
[ 22 ] * ( xx [ 31 ] * xx [ 124 ] + xx [ 118 ] * xx [ 70 ] + xx [ 140 ] * xx
[ 208 ] - xx [ 209 ] * xx [ 203 ] - xx [ 2 ] * xx [ 69 ] * xx [ 2 ] * xx [
105 ] * xx [ 106 ] * state [ 15 ] ) + xx [ 90 ] * xx [ 31 ] + xx [ 110 ] * xx
[ 70 ] + xx [ 2 ] * ( xx [ 94 ] * xx [ 214 ] + xx [ 211 ] * xx [ 96 ] ) - xx
[ 215 ] - ( xx [ 164 ] * xx [ 196 ] + xx [ 163 ] * xx [ 194 ] + xx [ 2 ] * xx
[ 94 ] * xx [ 210 ] * xx [ 108 ] ) + xx [ 2 ] * ( xx [ 2 ] * ( xx [ 139 ] *
xx [ 96 ] - xx [ 94 ] * xx [ 115 ] ) - xx [ 212 ] ) - ( xx [ 67 ] * xx [ 131
] - xx [ 134 ] * xx [ 199 ] ) + xx [ 2 ] * ( xx [ 220 ] * xx [ 95 ] + xx [ 97
] * xx [ 217 ] ) - xx [ 221 ] + xx [ 2 ] * ( xx [ 2 ] * ( xx [ 97 ] * xx [
117 ] - xx [ 113 ] * xx [ 95 ] ) - xx [ 218 ] ) - ( xx [ 226 ] + ( xx [ 224 ]
* xx [ 102 ] - xx [ 91 ] * xx [ 223 ] ) * xx [ 2 ] ) + xx [ 2 ] * ( ( xx [
100 ] * xx [ 102 ] + xx [ 91 ] * xx [ 84 ] ) * xx [ 2 ] - xx [ 228 ] ) ) ) +
xx [ 230 ] * xx [ 52 ] + xx [ 35 ] * xx [ 24 ] - xx [ 114 ] * xx [ 175 ] - xx
[ 53 ] * xx [ 88 ] - xx [ 20 ] * xx [ 104 ] - xx [ 149 ] * xx [ 50 ] ; memcpy
( xx + 18 , xx + 142 , 4 * sizeof ( double ) ) ; factorAndSolveSymmetric ( xx
+ 18 , 2 , xx + 22 , ii + 0 , xx + 172 , xx + 15 , xx + 175 ) ; xx [ 0 ] = (
xx [ 15 ] * xx [ 26 ] + xx [ 24 ] * xx [ 16 ] - xx [ 83 ] ) / xx [ 9 ] ; xx [
4 ] = xx [ 48 ] + xx [ 1 ] * xx [ 0 ] ; xx [ 8 ] = xx [ 4 ] * xx [ 7 ] ; xx [
10 ] = xx [ 4 ] - ( xx [ 29 ] + xx [ 8 ] * xx [ 7 ] ) * xx [ 2 ] ; xx [ 4 ] =
xx [ 10 ] - input [ 32 ] + xx [ 89 ] ; xx [ 12 ] = xx [ 99 ] + xx [ 2 ] * (
xx [ 5 ] * xx [ 8 ] - xx [ 101 ] ) ; xx [ 8 ] = xx [ 80 ] + xx [ 25 ] * xx [
0 ] + xx [ 12 ] * xx [ 45 ] - xx [ 44 ] * xx [ 10 ] - input [ 34 ] - xx [ 43
] ; xx [ 10 ] = ( xx [ 15 ] * xx [ 59 ] + xx [ 52 ] * xx [ 16 ] - ( xx [ 47 ]
+ xx [ 8 ] + xx [ 4 ] * xx [ 6 ] ) ) / xx [ 46 ] ; xx [ 18 ] = xx [ 4 ] + xx
[ 49 ] * xx [ 10 ] ; xx [ 4 ] = xx [ 12 ] - input [ 30 ] + xx [ 39 ] + xx [
37 ] * xx [ 10 ] ; xx [ 12 ] = xx [ 4 ] * xx [ 11 ] ; xx [ 19 ] = xx [ 18 ] *
xx [ 11 ] ; xx [ 20 ] = xx [ 18 ] - ( xx [ 14 ] * xx [ 12 ] + xx [ 19 ] * xx
[ 11 ] ) * xx [ 2 ] ; xx [ 18 ] = xx [ 20 ] - input [ 38 ] + xx [ 120 ] ; xx
[ 21 ] = xx [ 4 ] + xx [ 2 ] * ( xx [ 14 ] * xx [ 19 ] - xx [ 12 ] * xx [ 11
] ) ; xx [ 4 ] = ( xx [ 56 ] * xx [ 15 ] + xx [ 166 ] * xx [ 16 ] - ( xx [
128 ] + xx [ 8 ] + xx [ 41 ] * xx [ 10 ] + xx [ 21 ] * xx [ 82 ] - xx [ 74 ]
* xx [ 20 ] - input [ 40 ] - xx [ 27 ] + xx [ 18 ] * xx [ 6 ] ) ) / xx [ 73 ]
; xx [ 8 ] = xx [ 18 ] + xx [ 81 ] * xx [ 4 ] ; xx [ 12 ] = ( xx [ 15 ] - (
xx [ 8 ] - ( xx [ 17 ] * xx [ 8 ] * xx [ 17 ] - ( xx [ 21 ] - input [ 36 ] +
xx [ 33 ] + xx [ 61 ] * xx [ 4 ] ) * xx [ 17 ] * xx [ 13 ] ) * xx [ 2 ] ) ) /
xx [ 71 ] ; xx [ 8 ] = xx [ 17 ] * xx [ 12 ] ; xx [ 18 ] = xx [ 12 ] - xx [ 2
] * xx [ 17 ] * xx [ 8 ] ; xx [ 17 ] = xx [ 2 ] * xx [ 8 ] * xx [ 13 ] ; xx [
8 ] = xx [ 4 ] - ( xx [ 65 ] * xx [ 18 ] + xx [ 77 ] * xx [ 17 ] ) ; xx [ 4 ]
= xx [ 119 ] + xx [ 17 ] + xx [ 82 ] * xx [ 8 ] ; xx [ 13 ] = xx [ 18 ] + xx
[ 6 ] * xx [ 8 ] - xx [ 74 ] * xx [ 8 ] ; xx [ 17 ] = xx [ 11 ] * xx [ 13 ] ;
xx [ 18 ] = xx [ 4 ] * xx [ 11 ] ; xx [ 19 ] = xx [ 4 ] - ( xx [ 14 ] * xx [
17 ] + xx [ 18 ] * xx [ 11 ] ) * xx [ 2 ] ; xx [ 4 ] = xx [ 13 ] + xx [ 2 ] *
( xx [ 14 ] * xx [ 18 ] - xx [ 17 ] * xx [ 11 ] ) ; xx [ 11 ] = xx [ 10 ] - (
xx [ 55 ] * xx [ 8 ] + xx [ 68 ] * xx [ 19 ] + xx [ 4 ] * xx [ 62 ] ) ; xx [
10 ] = xx [ 8 ] + xx [ 11 ] ; xx [ 13 ] = xx [ 4 ] + xx [ 6 ] * xx [ 11 ] +
xx [ 107 ] - xx [ 44 ] * xx [ 10 ] ; xx [ 4 ] = xx [ 0 ] - ( xx [ 10 ] * xx [
78 ] + ( xx [ 13 ] + xx [ 2 ] * ( xx [ 5 ] * ( xx [ 19 ] + xx [ 116 ] + xx [
10 ] * xx [ 45 ] ) * xx [ 7 ] - xx [ 7 ] * xx [ 13 ] * xx [ 7 ] ) ) * xx [ 42
] ) ; xx [ 0 ] = ( xx [ 86 ] + xx [ 15 ] * xx [ 109 ] + xx [ 50 ] * xx [ 16 ]
) / xx [ 9 ] ; xx [ 5 ] = xx [ 234 ] - xx [ 1 ] * xx [ 0 ] ; xx [ 1 ] = xx [
85 ] * xx [ 5 ] ; xx [ 7 ] = xx [ 5 ] + xx [ 2 ] * ( xx [ 75 ] - xx [ 1 ] *
xx [ 85 ] ) ; xx [ 5 ] = xx [ 7 ] - input [ 14 ] + xx [ 236 ] ; xx [ 9 ] = xx
[ 237 ] - ( xx [ 103 ] * xx [ 1 ] + xx [ 112 ] ) * xx [ 2 ] ; xx [ 1 ] = xx [
238 ] - ( input [ 16 ] + xx [ 232 ] - xx [ 25 ] * xx [ 0 ] + xx [ 90 ] * xx [
7 ] - xx [ 110 ] * xx [ 9 ] ) ; xx [ 7 ] = ( xx [ 127 ] + xx [ 5 ] * xx [ 6 ]
- xx [ 1 ] + xx [ 104 ] * xx [ 16 ] - xx [ 15 ] * xx [ 121 ] ) / xx [ 60 ] ;
xx [ 13 ] = xx [ 5 ] - xx [ 7 ] * xx [ 132 ] ; xx [ 5 ] = xx [ 9 ] - input [
12 ] + xx [ 34 ] - xx [ 125 ] * xx [ 7 ] ; xx [ 9 ] = xx [ 92 ] * xx [ 5 ] ;
xx [ 14 ] = xx [ 92 ] * xx [ 13 ] ; xx [ 17 ] = xx [ 13 ] + xx [ 2 ] * ( xx [
98 ] * xx [ 9 ] - xx [ 14 ] * xx [ 92 ] ) ; xx [ 13 ] = xx [ 17 ] - input [ 8
] + xx [ 231 ] ; xx [ 18 ] = xx [ 5 ] - ( xx [ 98 ] * xx [ 14 ] + xx [ 9 ] *
xx [ 92 ] ) * xx [ 2 ] ; xx [ 5 ] = xx [ 1 ] - xx [ 7 ] * xx [ 40 ] - ( xx [
131 ] * xx [ 17 ] - xx [ 134 ] * xx [ 18 ] ) - input [ 10 ] + xx [ 28 ] ; xx
[ 1 ] = ( xx [ 160 ] + xx [ 13 ] * xx [ 6 ] - xx [ 5 ] + xx [ 15 ] * xx [ 147
] + xx [ 88 ] * xx [ 16 ] ) / xx [ 151 ] ; xx [ 9 ] = xx [ 13 ] - xx [ 1 ] *
xx [ 159 ] ; xx [ 13 ] = xx [ 18 ] - input [ 6 ] + xx [ 136 ] - xx [ 150 ] *
xx [ 1 ] ; xx [ 14 ] = xx [ 93 ] * xx [ 13 ] ; xx [ 17 ] = xx [ 93 ] * xx [ 9
] ; xx [ 18 ] = xx [ 9 ] + xx [ 2 ] * ( xx [ 3 ] * xx [ 14 ] - xx [ 17 ] * xx
[ 93 ] ) ; xx [ 9 ] = ( xx [ 38 ] + xx [ 6 ] * ( xx [ 18 ] - input [ 2 ] - xx
[ 162 ] ) - ( xx [ 5 ] - xx [ 1 ] * xx [ 129 ] - ( xx [ 163 ] * xx [ 18 ] -
xx [ 164 ] * ( xx [ 13 ] - ( xx [ 3 ] * xx [ 17 ] + xx [ 14 ] * xx [ 93 ] ) *
xx [ 2 ] ) ) - input [ 4 ] + xx [ 169 ] ) + xx [ 15 ] * xx [ 170 ] + xx [ 114
] * xx [ 16 ] ) / xx [ 171 ] ; xx [ 5 ] = xx [ 123 ] - xx [ 9 ] * xx [ 164 ]
; xx [ 13 ] = xx [ 5 ] * xx [ 93 ] ; xx [ 14 ] = xx [ 6 ] * xx [ 9 ] + xx [
163 ] * xx [ 9 ] ; xx [ 15 ] = xx [ 14 ] * xx [ 93 ] ; xx [ 16 ] = xx [ 2 ] *
( xx [ 13 ] * xx [ 93 ] - xx [ 3 ] * xx [ 15 ] ) - xx [ 5 ] ; xx [ 5 ] = ( xx
[ 3 ] * xx [ 13 ] + xx [ 15 ] * xx [ 93 ] ) * xx [ 2 ] - xx [ 14 ] ; xx [ 3 ]
= xx [ 1 ] + xx [ 9 ] * xx [ 185 ] + xx [ 155 ] * xx [ 16 ] + xx [ 174 ] * xx
[ 5 ] ; xx [ 1 ] = xx [ 9 ] + xx [ 3 ] ; xx [ 13 ] = xx [ 16 ] + xx [ 244 ] +
xx [ 1 ] * xx [ 134 ] ; xx [ 14 ] = xx [ 5 ] - xx [ 3 ] * xx [ 6 ] + xx [ 58
] - xx [ 1 ] * xx [ 131 ] ; xx [ 5 ] = xx [ 92 ] * xx [ 14 ] ; xx [ 15 ] = xx
[ 92 ] * xx [ 13 ] ; xx [ 16 ] = xx [ 13 ] + xx [ 2 ] * ( xx [ 98 ] * xx [ 5
] - xx [ 15 ] * xx [ 92 ] ) ; xx [ 13 ] = xx [ 14 ] - ( xx [ 98 ] * xx [ 15 ]
+ xx [ 5 ] * xx [ 92 ] ) * xx [ 2 ] ; xx [ 5 ] = xx [ 7 ] + xx [ 1 ] * xx [
165 ] + xx [ 16 ] * xx [ 122 ] + xx [ 153 ] * xx [ 13 ] ; xx [ 7 ] = xx [ 1 ]
+ xx [ 5 ] ; xx [ 1 ] = xx [ 13 ] - xx [ 5 ] * xx [ 6 ] + xx [ 233 ] - xx [ 7
] * xx [ 90 ] ; xx [ 6 ] = xx [ 0 ] + xx [ 42 ] * ( xx [ 1 ] - ( xx [ 103 ] *
xx [ 85 ] * ( xx [ 16 ] + xx [ 240 ] + xx [ 7 ] * xx [ 110 ] ) + xx [ 85 ] *
xx [ 1 ] * xx [ 85 ] ) * xx [ 2 ] ) - xx [ 7 ] * xx [ 78 ] ; deriv [ 0 ] =
state [ 1 ] ; deriv [ 1 ] = xx [ 12 ] ; deriv [ 2 ] = state [ 3 ] ; deriv [ 3
] = xx [ 8 ] ; deriv [ 4 ] = state [ 5 ] ; deriv [ 5 ] = xx [ 11 ] ; deriv [
6 ] = state [ 7 ] ; deriv [ 7 ] = xx [ 4 ] ; deriv [ 8 ] = state [ 9 ] ;
deriv [ 9 ] = - xx [ 9 ] ; deriv [ 10 ] = state [ 11 ] ; deriv [ 11 ] = - xx
[ 3 ] ; deriv [ 12 ] = state [ 13 ] ; deriv [ 13 ] = - xx [ 5 ] ; deriv [ 14
] = state [ 15 ] ; deriv [ 15 ] = - xx [ 6 ] ; deriv [ 16 ] = state [ 17 ] ;
deriv [ 17 ] = xx [ 7 ] + xx [ 6 ] - ( xx [ 10 ] + xx [ 4 ] ) ; errorResult [
0 ] = 0.0 ; return NULL ; } PmfMessageId
MagneticBasket_Simscape_Optimizer_dda62cd9_1_numJacPerturbLoBounds ( const
RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags , const double
* state , const int * modeVector , const double * input , const double *
inputDot , const double * inputDdot , const double * discreteState , double *
bounds , double * errorResult , NeuDiagnosticManager * neDiagMgr ) { const
double * rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv ->
mInts . mValues ; double xx [ 2 ] ; ( void ) rtdvd ; ( void ) rtdvi ; ( void
) eqnEnableFlags ; ( void ) state ; ( void ) modeVector ; ( void ) input ; (
void ) inputDot ; ( void ) inputDdot ; ( void ) discreteState ; ( void )
neDiagMgr ; xx [ 0 ] = 1.0e-9 ; xx [ 1 ] = 1.0e-8 ; bounds [ 0 ] = xx [ 0 ] ;
bounds [ 1 ] = xx [ 0 ] ; bounds [ 2 ] = xx [ 1 ] ; bounds [ 3 ] = xx [ 1 ] ;
bounds [ 4 ] = xx [ 1 ] ; bounds [ 5 ] = xx [ 1 ] ; bounds [ 6 ] = xx [ 1 ] ;
bounds [ 7 ] = xx [ 1 ] ; bounds [ 8 ] = xx [ 1 ] ; bounds [ 9 ] = xx [ 1 ] ;
bounds [ 10 ] = xx [ 1 ] ; bounds [ 11 ] = xx [ 1 ] ; bounds [ 12 ] = xx [ 1
] ; bounds [ 13 ] = xx [ 1 ] ; bounds [ 14 ] = xx [ 1 ] ; bounds [ 15 ] = xx
[ 1 ] ; bounds [ 16 ] = xx [ 1 ] ; bounds [ 17 ] = xx [ 1 ] ; errorResult [ 0
] = 0.0 ; return NULL ; } PmfMessageId
MagneticBasket_Simscape_Optimizer_dda62cd9_1_numJacPerturbHiBounds ( const
RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags , const double
* state , const int * modeVector , const double * input , const double *
inputDot , const double * inputDdot , const double * discreteState , double *
bounds , double * errorResult , NeuDiagnosticManager * neDiagMgr ) { const
double * rtdvd = rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv ->
mInts . mValues ; double xx [ 2 ] ; ( void ) rtdvd ; ( void ) rtdvi ; ( void
) eqnEnableFlags ; ( void ) state ; ( void ) modeVector ; ( void ) input ; (
void ) inputDot ; ( void ) inputDdot ; ( void ) discreteState ; ( void )
neDiagMgr ; xx [ 0 ] = + pmf_get_inf ( ) ; xx [ 1 ] = 1.0 ; bounds [ 0 ] = xx
[ 0 ] ; bounds [ 1 ] = xx [ 0 ] ; bounds [ 2 ] = xx [ 1 ] ; bounds [ 3 ] = xx
[ 0 ] ; bounds [ 4 ] = xx [ 1 ] ; bounds [ 5 ] = xx [ 0 ] ; bounds [ 6 ] = xx
[ 1 ] ; bounds [ 7 ] = xx [ 0 ] ; bounds [ 8 ] = xx [ 1 ] ; bounds [ 9 ] = xx
[ 0 ] ; bounds [ 10 ] = xx [ 1 ] ; bounds [ 11 ] = xx [ 0 ] ; bounds [ 12 ] =
xx [ 1 ] ; bounds [ 13 ] = xx [ 0 ] ; bounds [ 14 ] = xx [ 1 ] ; bounds [ 15
] = xx [ 0 ] ; bounds [ 16 ] = xx [ 1 ] ; bounds [ 17 ] = xx [ 0 ] ;
errorResult [ 0 ] = 0.0 ; return NULL ; }
