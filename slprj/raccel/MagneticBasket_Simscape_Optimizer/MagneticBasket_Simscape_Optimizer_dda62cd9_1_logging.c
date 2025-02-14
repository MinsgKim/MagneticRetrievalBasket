#include <math.h>
#include <string.h>
#include "pm_std.h"
#include "sm_std.h"
#include "ne_std.h"
#include "ne_dae.h"
#include "sm_ssci_run_time_errors.h"
#include "sm_RuntimeDerivedValuesBundle.h"
#include "MagneticBasket_Simscape_Optimizer_dda62cd9_1_geometries.h"
PmfMessageId MagneticBasket_Simscape_Optimizer_dda62cd9_1_recordLog ( const
RuntimeDerivedValuesBundle * rtdv , const int * eqnEnableFlags , const double
* state , const int * modeVector , const double * input , const double *
inputDot , const double * inputDdot , double * logVector , double *
errorResult , NeuDiagnosticManager * neDiagMgr ) { const double * rtdvd =
rtdv -> mDoubles . mValues ; const int * rtdvi = rtdv -> mInts . mValues ;
int ii [ 2 ] ; double xx [ 248 ] ; ( void ) rtdvd ; ( void ) rtdvi ; ( void )
eqnEnableFlags ; ( void ) modeVector ; ( void ) inputDot ; ( void ) inputDdot
; ( void ) neDiagMgr ; xx [ 0 ] = 57.29577951308232 ; xx [ 1 ] = 1.0 ; xx [ 2
] = 3.000000000000001e-9 ; xx [ 3 ] = 2.0 ; xx [ 4 ] = 0.5 ; xx [ 5 ] = xx [
4 ] * state [ 6 ] ; xx [ 6 ] = cos ( xx [ 5 ] ) ; xx [ 7 ] = 1.0e-3 ; xx [ 8
] = sin ( xx [ 5 ] ) ; xx [ 5 ] = xx [ 7 ] * xx [ 8 ] ; xx [ 9 ] = xx [ 6 ] *
xx [ 5 ] ; xx [ 10 ] = xx [ 3 ] * xx [ 9 ] ; xx [ 11 ] = xx [ 4 ] * state [ 4
] ; xx [ 12 ] = sin ( xx [ 11 ] ) ; xx [ 13 ] = 0.7071067811865476 ; xx [ 14
] = xx [ 4 ] * state [ 2 ] ; xx [ 15 ] = xx [ 13 ] * cos ( xx [ 14 ] ) ; xx [
16 ] = xx [ 13 ] * sin ( xx [ 14 ] ) ; xx [ 14 ] = xx [ 15 ] + xx [ 16 ] ; xx
[ 17 ] = xx [ 13 ] * xx [ 14 ] ; xx [ 18 ] = xx [ 15 ] - xx [ 16 ] ; xx [ 15
] = xx [ 18 ] * xx [ 13 ] ; xx [ 16 ] = xx [ 17 ] + xx [ 15 ] ; xx [ 19 ] =
xx [ 17 ] - xx [ 15 ] ; xx [ 15 ] = cos ( xx [ 11 ] ) ; xx [ 11 ] = xx [ 12 ]
* xx [ 16 ] + xx [ 19 ] * xx [ 15 ] ; xx [ 17 ] = xx [ 7 ] - xx [ 3 ] * xx [
5 ] * xx [ 8 ] ; xx [ 5 ] = xx [ 11 ] * xx [ 17 ] ; xx [ 20 ] = xx [ 19 ] *
xx [ 12 ] - xx [ 15 ] * xx [ 16 ] ; xx [ 21 ] = xx [ 11 ] * xx [ 10 ] ; xx [
22 ] = xx [ 6 ] * xx [ 20 ] + xx [ 11 ] * xx [ 8 ] ; xx [ 23 ] = 1.0e-3 ; xx
[ 24 ] = xx [ 11 ] * xx [ 6 ] - xx [ 8 ] * xx [ 20 ] ; xx [ 25 ] = xx [ 23 ]
* xx [ 24 ] ; xx [ 26 ] = xx [ 3 ] * xx [ 22 ] * xx [ 25 ] ; xx [ 27 ] = xx [
10 ] - ( xx [ 5 ] * xx [ 20 ] + xx [ 11 ] * xx [ 21 ] ) * xx [ 3 ] - xx [ 26
] ; xx [ 10 ] = 4.062500000000001e-12 ; xx [ 28 ] = xx [ 27 ] / xx [ 10 ] ;
xx [ 29 ] = xx [ 2 ] * xx [ 28 ] ; xx [ 30 ] = xx [ 29 ] * xx [ 8 ] ; xx [ 31
] = xx [ 29 ] - xx [ 3 ] * xx [ 30 ] * xx [ 8 ] ; xx [ 29 ] = 3.0e-6 ; xx [
32 ] = xx [ 6 ] * xx [ 8 ] ; xx [ 33 ] = xx [ 3 ] * xx [ 32 ] ; xx [ 34 ] =
xx [ 29 ] * xx [ 33 ] ; xx [ 35 ] = 7.846153846153842e-7 ; xx [ 36 ] = xx [ 6
] * xx [ 6 ] ; xx [ 37 ] = xx [ 3 ] * xx [ 36 ] - xx [ 1 ] ; xx [ 38 ] = xx [
35 ] * xx [ 37 ] ; xx [ 39 ] = xx [ 34 ] * xx [ 33 ] + xx [ 38 ] * xx [ 37 ]
; xx [ 40 ] = xx [ 29 ] + xx [ 39 ] ; xx [ 41 ] = 7.846153846153845e-10 ; xx
[ 42 ] = ( xx [ 36 ] + xx [ 8 ] * xx [ 8 ] ) * xx [ 3 ] - xx [ 1 ] ; xx [ 36
] = xx [ 41 ] * xx [ 37 ] * xx [ 42 ] ; xx [ 43 ] = xx [ 23 ] * xx [ 8 ] ; xx
[ 44 ] = 2.0e-3 ; xx [ 45 ] = xx [ 3 ] * xx [ 43 ] * xx [ 8 ] - xx [ 44 ] ;
xx [ 46 ] = xx [ 3 ] * xx [ 6 ] * xx [ 43 ] ; xx [ 43 ] = xx [ 35 ] * xx [ 33
] ; xx [ 47 ] = xx [ 29 ] * xx [ 37 ] ; xx [ 48 ] = xx [ 43 ] * xx [ 37 ] -
xx [ 47 ] * xx [ 33 ] ; xx [ 49 ] = xx [ 45 ] * xx [ 39 ] - xx [ 46 ] * xx [
48 ] ; xx [ 39 ] = xx [ 36 ] + xx [ 49 ] ; xx [ 50 ] = xx [ 40 ] * xx [ 7 ] -
xx [ 39 ] ; xx [ 51 ] = xx [ 45 ] * xx [ 12 ] ; xx [ 52 ] = xx [ 46 ] * xx [
12 ] ; xx [ 53 ] = xx [ 46 ] - ( xx [ 15 ] * xx [ 51 ] + xx [ 52 ] * xx [ 12
] ) * xx [ 3 ] ; xx [ 54 ] = xx [ 7 ] * xx [ 12 ] ; xx [ 55 ] = xx [ 15 ] *
xx [ 54 ] ; xx [ 56 ] = xx [ 53 ] + xx [ 3 ] * xx [ 55 ] ; xx [ 57 ] = xx [
51 ] * xx [ 12 ] ; xx [ 51 ] = xx [ 15 ] * xx [ 52 ] ; xx [ 52 ] = xx [ 3 ] *
xx [ 54 ] * xx [ 12 ] ; xx [ 54 ] = xx [ 3 ] * ( xx [ 57 ] - xx [ 51 ] ) - (
xx [ 45 ] + xx [ 52 ] ) + xx [ 7 ] ; xx [ 58 ] = xx [ 19 ] * xx [ 54 ] ; xx [
59 ] = xx [ 19 ] * xx [ 56 ] ; xx [ 60 ] = xx [ 56 ] + xx [ 3 ] * ( xx [ 58 ]
* xx [ 16 ] - xx [ 19 ] * xx [ 59 ] ) - xx [ 26 ] ; xx [ 26 ] = 1.0625e-12 ;
xx [ 56 ] = xx [ 3 ] * xx [ 6 ] * xx [ 30 ] ; xx [ 30 ] = xx [ 26 ] * xx [ 28
] + xx [ 46 ] * xx [ 56 ] - xx [ 45 ] * xx [ 31 ] ; xx [ 61 ] =
7.846153846153846e-13 ; xx [ 62 ] = xx [ 41 ] * xx [ 33 ] * xx [ 42 ] ; xx [
63 ] = xx [ 46 ] * xx [ 62 ] - xx [ 45 ] * xx [ 36 ] ; xx [ 36 ] = xx [ 38 ]
* xx [ 33 ] - xx [ 34 ] * xx [ 37 ] ; xx [ 34 ] = xx [ 47 ] * xx [ 37 ] + xx
[ 43 ] * xx [ 33 ] ; xx [ 38 ] = xx [ 45 ] * xx [ 36 ] - xx [ 34 ] * xx [ 46
] ; xx [ 43 ] = xx [ 26 ] + xx [ 61 ] * xx [ 42 ] * xx [ 42 ] - xx [ 63 ] -
xx [ 63 ] - ( xx [ 46 ] * xx [ 38 ] - xx [ 45 ] * xx [ 49 ] ) ; xx [ 42 ] =
xx [ 43 ] - xx [ 7 ] * xx [ 39 ] ; xx [ 47 ] = xx [ 42 ] + xx [ 50 ] * xx [ 7
] ; ii [ 0 ] = factorSymmetricPosDef ( xx + 47 , 1 , xx + 49 ) ; if ( ii [ 0
] != 0 ) { return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassBase" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint7' has a degenerate mass distribution on its base side."
, neDiagMgr ) ; } xx [ 49 ] = ( xx [ 60 ] - ( xx [ 30 ] + xx [ 7 ] * xx [ 31
] ) ) / xx [ 47 ] ; xx [ 63 ] = xx [ 31 ] + xx [ 50 ] * xx [ 49 ] ; xx [ 31 ]
= xx [ 62 ] + xx [ 38 ] ; xx [ 38 ] = xx [ 7 ] * xx [ 36 ] - xx [ 31 ] ; xx [
62 ] = xx [ 56 ] + xx [ 38 ] * xx [ 49 ] ; xx [ 56 ] = xx [ 62 ] * xx [ 12 ]
; xx [ 64 ] = xx [ 63 ] * xx [ 12 ] ; xx [ 65 ] = xx [ 63 ] - ( xx [ 15 ] *
xx [ 56 ] + xx [ 64 ] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 63 ] = xx [ 15 ] * xx [
15 ] ; xx [ 66 ] = ( xx [ 63 ] + xx [ 12 ] * xx [ 12 ] ) * xx [ 3 ] - xx [ 1
] ; xx [ 67 ] = xx [ 15 ] * xx [ 12 ] ; xx [ 68 ] = xx [ 3 ] * xx [ 67 ] ; xx
[ 69 ] = xx [ 38 ] / xx [ 47 ] ; xx [ 70 ] = xx [ 31 ] + xx [ 42 ] * xx [ 69
] ; xx [ 71 ] = xx [ 3 ] * xx [ 63 ] - xx [ 1 ] ; xx [ 63 ] = xx [ 50 ] / xx
[ 47 ] ; xx [ 72 ] = xx [ 39 ] + xx [ 42 ] * xx [ 63 ] ; xx [ 73 ] = xx [ 66
] * ( xx [ 68 ] * xx [ 70 ] - xx [ 71 ] * xx [ 72 ] ) ; xx [ 74 ] = xx [ 23 ]
* xx [ 12 ] ; xx [ 75 ] = xx [ 3 ] * xx [ 74 ] * xx [ 12 ] - xx [ 44 ] ; xx [
76 ] = xx [ 40 ] - xx [ 50 ] * xx [ 63 ] ; xx [ 77 ] = xx [ 50 ] * xx [ 69 ]
; xx [ 78 ] = xx [ 48 ] - xx [ 77 ] ; xx [ 79 ] = xx [ 76 ] * xx [ 71 ] - xx
[ 68 ] * xx [ 78 ] ; xx [ 80 ] = xx [ 36 ] - xx [ 77 ] ; xx [ 77 ] = xx [ 29
] + xx [ 34 ] ; xx [ 34 ] = xx [ 77 ] - xx [ 38 ] * xx [ 69 ] ; xx [ 81 ] =
xx [ 71 ] * xx [ 80 ] - xx [ 68 ] * xx [ 34 ] ; xx [ 82 ] = xx [ 71 ] * xx [
79 ] - xx [ 68 ] * xx [ 81 ] ; xx [ 83 ] = xx [ 3 ] * xx [ 15 ] * xx [ 74 ] ;
xx [ 74 ] = xx [ 71 ] * xx [ 78 ] + xx [ 68 ] * xx [ 76 ] ; xx [ 76 ] = xx [
34 ] * xx [ 71 ] + xx [ 68 ] * xx [ 80 ] ; xx [ 34 ] = xx [ 74 ] * xx [ 71 ]
- xx [ 76 ] * xx [ 68 ] ; xx [ 78 ] = xx [ 75 ] * xx [ 82 ] - xx [ 83 ] * xx
[ 34 ] ; xx [ 80 ] = xx [ 73 ] - xx [ 78 ] ; xx [ 84 ] = xx [ 29 ] + xx [ 82
] ; xx [ 82 ] = xx [ 80 ] + xx [ 84 ] * xx [ 7 ] ; xx [ 85 ] = xx [ 15 ] * xx
[ 8 ] + xx [ 6 ] * xx [ 12 ] ; xx [ 86 ] = xx [ 15 ] * xx [ 6 ] - xx [ 12 ] *
xx [ 8 ] ; xx [ 87 ] = xx [ 85 ] * xx [ 16 ] + xx [ 19 ] * xx [ 86 ] ; xx [
88 ] = xx [ 87 ] * xx [ 23 ] ; xx [ 89 ] = xx [ 53 ] + xx [ 83 ] ; xx [ 53 ]
= xx [ 18 ] * xx [ 89 ] ; xx [ 90 ] = xx [ 45 ] + xx [ 3 ] * ( xx [ 51 ] - xx
[ 57 ] ) + xx [ 75 ] ; xx [ 51 ] = xx [ 90 ] * xx [ 18 ] ; xx [ 57 ] = xx [
18 ] * xx [ 7 ] ; xx [ 91 ] = xx [ 89 ] - ( xx [ 18 ] * xx [ 53 ] - xx [ 51 ]
* xx [ 14 ] ) * xx [ 3 ] - xx [ 3 ] * xx [ 57 ] * xx [ 14 ] ; xx [ 89 ] = xx
[ 3 ] * ( xx [ 18 ] * xx [ 51 ] + xx [ 53 ] * xx [ 14 ] ) - ( xx [ 90 ] + xx
[ 3 ] * xx [ 18 ] * xx [ 57 ] ) + xx [ 7 ] ; xx [ 51 ] = xx [ 13 ] * xx [ 89
] * xx [ 13 ] ; xx [ 53 ] = xx [ 13 ] * xx [ 91 ] * xx [ 13 ] ; xx [ 57 ] =
xx [ 3 ] * xx [ 88 ] * ( xx [ 86 ] * xx [ 16 ] - xx [ 85 ] * xx [ 19 ] ) + xx
[ 91 ] + xx [ 3 ] * ( xx [ 51 ] - xx [ 53 ] ) ; xx [ 85 ] = xx [ 62 ] + xx [
3 ] * ( xx [ 15 ] * xx [ 64 ] - xx [ 56 ] * xx [ 12 ] ) ; xx [ 56 ] = xx [ 42
] / xx [ 47 ] ; xx [ 62 ] = ( xx [ 71 ] * xx [ 70 ] + xx [ 68 ] * xx [ 72 ] )
* xx [ 66 ] ; xx [ 64 ] = xx [ 75 ] * xx [ 73 ] + xx [ 62 ] * xx [ 83 ] ; xx
[ 70 ] = xx [ 71 ] * xx [ 81 ] + xx [ 68 ] * xx [ 79 ] ; xx [ 72 ] = xx [ 76
] * xx [ 71 ] + xx [ 74 ] * xx [ 68 ] ; xx [ 73 ] = xx [ 75 ] * xx [ 70 ] -
xx [ 72 ] * xx [ 83 ] ; xx [ 74 ] = ( xx [ 43 ] - xx [ 42 ] * xx [ 56 ] ) *
xx [ 66 ] * xx [ 66 ] - xx [ 64 ] - xx [ 64 ] - ( xx [ 83 ] * xx [ 73 ] - xx
[ 75 ] * xx [ 78 ] ) + xx [ 7 ] * xx [ 80 ] + xx [ 82 ] * xx [ 7 ] + xx [ 26
] ; ii [ 0 ] = factorSymmetricPosDef ( xx + 74 , 1 , xx + 43 ) ; if ( ii [ 0
] != 0 ) { return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassBase" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint8' has a degenerate mass distribution on its base side."
, neDiagMgr ) ; } xx [ 43 ] = ( xx [ 57 ] - ( xx [ 30 ] + xx [ 42 ] * xx [ 49
] + xx [ 85 ] * xx [ 83 ] - xx [ 75 ] * xx [ 65 ] + xx [ 7 ] * xx [ 65 ] ) )
/ xx [ 74 ] ; xx [ 30 ] = xx [ 65 ] + xx [ 82 ] * xx [ 43 ] ; xx [ 64 ] = xx
[ 62 ] + xx [ 73 ] ; xx [ 62 ] = xx [ 70 ] * xx [ 7 ] - xx [ 64 ] ; xx [ 65 ]
= xx [ 3 ] * xx [ 14 ] * xx [ 14 ] - xx [ 1 ] ; xx [ 66 ] = xx [ 82 ] / xx [
74 ] ; xx [ 73 ] = xx [ 18 ] * xx [ 14 ] ; xx [ 76 ] = xx [ 3 ] * xx [ 73 ] ;
xx [ 78 ] = xx [ 62 ] / xx [ 74 ] ; xx [ 79 ] = xx [ 82 ] * xx [ 78 ] ; xx [
80 ] = xx [ 29 ] + xx [ 72 ] ; xx [ 72 ] = 4.500000000000001e-6 + xx [ 65 ] *
( ( xx [ 84 ] - xx [ 82 ] * xx [ 66 ] ) * xx [ 65 ] + xx [ 76 ] * ( xx [ 34 ]
- xx [ 79 ] ) ) + xx [ 76 ] * ( ( xx [ 70 ] - xx [ 79 ] ) * xx [ 65 ] + xx [
76 ] * ( xx [ 80 ] - xx [ 62 ] * xx [ 78 ] ) ) ; ii [ 0 ] =
factorSymmetricPosDef ( xx + 72 , 1 , xx + 70 ) ; if ( ii [ 0 ] != 0 ) {
return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Prismatic Joint1' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 70 ] = ( xx [ 1 ] - ( xx [ 30 ] - ( xx [ 18 ] * xx [
18 ] * xx [ 30 ] - ( xx [ 85 ] + xx [ 62 ] * xx [ 43 ] ) * xx [ 18 ] * xx [
14 ] ) * xx [ 3 ] ) ) / xx [ 72 ] ; xx [ 30 ] = xx [ 18 ] * xx [ 70 ] ; xx [
76 ] = xx [ 70 ] - xx [ 3 ] * xx [ 18 ] * xx [ 30 ] ; xx [ 79 ] = xx [ 3 ] *
xx [ 30 ] * xx [ 14 ] ; xx [ 30 ] = xx [ 43 ] - ( xx [ 66 ] * xx [ 76 ] + xx
[ 78 ] * xx [ 79 ] ) ; xx [ 43 ] = xx [ 83 ] * xx [ 30 ] + xx [ 79 ] ; xx [
79 ] = xx [ 76 ] + xx [ 7 ] * xx [ 30 ] - xx [ 75 ] * xx [ 30 ] ; xx [ 76 ] =
xx [ 12 ] * xx [ 79 ] ; xx [ 81 ] = xx [ 12 ] * xx [ 43 ] ; xx [ 84 ] = xx [
43 ] - ( xx [ 15 ] * xx [ 76 ] + xx [ 81 ] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 43
] = xx [ 79 ] + xx [ 3 ] * ( xx [ 15 ] * xx [ 81 ] - xx [ 76 ] * xx [ 12 ] )
; xx [ 76 ] = xx [ 49 ] - ( xx [ 56 ] * xx [ 30 ] + xx [ 69 ] * xx [ 84 ] +
xx [ 43 ] * xx [ 63 ] ) ; xx [ 49 ] = xx [ 30 ] + xx [ 76 ] ; xx [ 79 ] =
0.2615384615384614 ; xx [ 81 ] = xx [ 43 ] + xx [ 7 ] * xx [ 76 ] - xx [ 45 ]
* xx [ 49 ] ; xx [ 43 ] = 738.4615384615385 ; xx [ 85 ] = xx [ 4 ] * state [
14 ] ; xx [ 86 ] = sin ( xx [ 85 ] ) ; xx [ 90 ] = xx [ 23 ] * xx [ 86 ] ; xx
[ 91 ] = xx [ 44 ] - xx [ 3 ] * xx [ 90 ] * xx [ 86 ] ; xx [ 92 ] = xx [ 4 ]
* state [ 12 ] ; xx [ 93 ] = sin ( xx [ 92 ] ) ; xx [ 94 ] = xx [ 4 ] * state
[ 8 ] ; xx [ 95 ] = cos ( xx [ 94 ] ) ; xx [ 96 ] = xx [ 4 ] * state [ 10 ] ;
xx [ 4 ] = cos ( xx [ 96 ] ) ; xx [ 97 ] = sin ( xx [ 94 ] ) ; xx [ 94 ] =
sin ( xx [ 96 ] ) ; xx [ 96 ] = xx [ 95 ] * xx [ 4 ] - xx [ 97 ] * xx [ 94 ]
; xx [ 98 ] = xx [ 95 ] * xx [ 94 ] + xx [ 4 ] * xx [ 97 ] ; xx [ 99 ] = cos
( xx [ 92 ] ) ; xx [ 92 ] = xx [ 93 ] * xx [ 96 ] + xx [ 98 ] * xx [ 99 ] ;
xx [ 100 ] = xx [ 7 ] * xx [ 86 ] ; xx [ 101 ] = xx [ 7 ] - xx [ 3 ] * xx [
100 ] * xx [ 86 ] ; xx [ 102 ] = xx [ 92 ] * xx [ 101 ] ; xx [ 103 ] = xx [
98 ] * xx [ 93 ] - xx [ 99 ] * xx [ 96 ] ; xx [ 104 ] = cos ( xx [ 85 ] ) ;
xx [ 85 ] = xx [ 104 ] * xx [ 100 ] ; xx [ 100 ] = xx [ 3 ] * xx [ 85 ] ; xx
[ 105 ] = xx [ 92 ] * xx [ 100 ] ; xx [ 106 ] = xx [ 104 ] * xx [ 103 ] + xx
[ 92 ] * xx [ 86 ] ; xx [ 107 ] = xx [ 86 ] * xx [ 103 ] - xx [ 92 ] * xx [
104 ] ; xx [ 108 ] = xx [ 23 ] * xx [ 107 ] ; xx [ 109 ] = xx [ 3 ] * xx [
106 ] * xx [ 108 ] ; xx [ 110 ] = ( xx [ 102 ] * xx [ 103 ] + xx [ 92 ] * xx
[ 105 ] ) * xx [ 3 ] - xx [ 100 ] - xx [ 109 ] ; xx [ 100 ] = xx [ 110 ] / xx
[ 10 ] ; xx [ 111 ] = xx [ 2 ] * xx [ 100 ] ; xx [ 112 ] = xx [ 111 ] * xx [
86 ] ; xx [ 113 ] = xx [ 3 ] * xx [ 112 ] * xx [ 86 ] - xx [ 111 ] ; xx [ 111
] = xx [ 3 ] * xx [ 104 ] * xx [ 90 ] ; xx [ 90 ] = xx [ 3 ] * xx [ 104 ] *
xx [ 112 ] ; xx [ 112 ] = xx [ 91 ] * xx [ 113 ] - xx [ 111 ] * xx [ 90 ] -
xx [ 26 ] * xx [ 100 ] ; xx [ 114 ] = xx [ 91 ] * xx [ 93 ] ; xx [ 115 ] = xx
[ 111 ] * xx [ 93 ] ; xx [ 116 ] = ( xx [ 99 ] * xx [ 114 ] - xx [ 115 ] * xx
[ 93 ] ) * xx [ 3 ] ; xx [ 117 ] = xx [ 7 ] * xx [ 93 ] ; xx [ 118 ] = xx [
99 ] * xx [ 117 ] ; xx [ 119 ] = xx [ 111 ] + xx [ 116 ] + xx [ 3 ] * xx [
118 ] ; xx [ 120 ] = xx [ 91 ] - xx [ 3 ] * ( xx [ 99 ] * xx [ 115 ] + xx [
114 ] * xx [ 93 ] ) ; xx [ 114 ] = xx [ 3 ] * xx [ 117 ] * xx [ 93 ] ; xx [
115 ] = xx [ 120 ] - xx [ 114 ] + xx [ 7 ] ; xx [ 117 ] = xx [ 98 ] * xx [
115 ] ; xx [ 121 ] = xx [ 98 ] * xx [ 119 ] ; xx [ 122 ] = xx [ 119 ] + ( xx
[ 117 ] * xx [ 96 ] - xx [ 98 ] * xx [ 121 ] ) * xx [ 3 ] + xx [ 109 ] ; xx [
109 ] = xx [ 104 ] * xx [ 86 ] ; xx [ 119 ] = xx [ 3 ] * xx [ 109 ] ; xx [
123 ] = xx [ 29 ] * xx [ 119 ] ; xx [ 124 ] = xx [ 104 ] * xx [ 104 ] ; xx [
125 ] = xx [ 3 ] * xx [ 124 ] - xx [ 1 ] ; xx [ 126 ] = xx [ 35 ] * xx [ 125
] ; xx [ 127 ] = xx [ 123 ] * xx [ 119 ] + xx [ 126 ] * xx [ 125 ] ; xx [ 128
] = xx [ 29 ] + xx [ 127 ] ; xx [ 129 ] = ( xx [ 124 ] + xx [ 86 ] * xx [ 86
] ) * xx [ 3 ] - xx [ 1 ] ; xx [ 124 ] = xx [ 41 ] * xx [ 125 ] * xx [ 129 ]
; xx [ 130 ] = xx [ 29 ] * xx [ 125 ] ; xx [ 131 ] = xx [ 35 ] * xx [ 119 ] ;
xx [ 35 ] = xx [ 130 ] * xx [ 119 ] - xx [ 131 ] * xx [ 125 ] ; xx [ 132 ] =
xx [ 127 ] * xx [ 91 ] - xx [ 111 ] * xx [ 35 ] ; xx [ 127 ] = xx [ 124 ] -
xx [ 132 ] ; xx [ 133 ] = xx [ 128 ] * xx [ 7 ] - xx [ 127 ] ; xx [ 134 ] =
xx [ 41 ] * xx [ 119 ] * xx [ 129 ] ; xx [ 41 ] = xx [ 91 ] * xx [ 124 ] + xx
[ 111 ] * xx [ 134 ] ; xx [ 124 ] = xx [ 123 ] * xx [ 125 ] - xx [ 126 ] * xx
[ 119 ] ; xx [ 123 ] = xx [ 130 ] * xx [ 125 ] + xx [ 131 ] * xx [ 119 ] ; xx
[ 126 ] = xx [ 91 ] * xx [ 124 ] - xx [ 123 ] * xx [ 111 ] ; xx [ 130 ] = xx
[ 26 ] + xx [ 61 ] * xx [ 129 ] * xx [ 129 ] - xx [ 41 ] - xx [ 41 ] + xx [
132 ] * xx [ 91 ] - xx [ 126 ] * xx [ 111 ] ; xx [ 41 ] = xx [ 7 ] * xx [ 127
] - xx [ 130 ] ; xx [ 61 ] = xx [ 7 ] * xx [ 133 ] - xx [ 41 ] ; ii [ 0 ] =
factorSymmetricPosDef ( xx + 61 , 1 , xx + 129 ) ; if ( ii [ 0 ] != 0 ) {
return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint3' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 129 ] = ( xx [ 112 ] + xx [ 7 ] * xx [ 113 ] - xx [
122 ] ) / xx [ 61 ] ; xx [ 131 ] = xx [ 23 ] * xx [ 93 ] ; xx [ 132 ] = xx [
44 ] - xx [ 3 ] * xx [ 131 ] * xx [ 93 ] ; xx [ 135 ] = xx [ 113 ] - xx [ 129
] * xx [ 133 ] ; xx [ 113 ] = xx [ 134 ] + xx [ 126 ] ; xx [ 126 ] = xx [ 113
] + xx [ 7 ] * xx [ 124 ] ; xx [ 134 ] = xx [ 90 ] - xx [ 126 ] * xx [ 129 ]
; xx [ 90 ] = xx [ 93 ] * xx [ 134 ] ; xx [ 136 ] = xx [ 93 ] * xx [ 135 ] ;
xx [ 137 ] = xx [ 135 ] + xx [ 3 ] * ( xx [ 99 ] * xx [ 90 ] - xx [ 136 ] *
xx [ 93 ] ) ; xx [ 135 ] = xx [ 3 ] * xx [ 99 ] * xx [ 131 ] ; xx [ 131 ] =
xx [ 134 ] - ( xx [ 99 ] * xx [ 136 ] + xx [ 90 ] * xx [ 93 ] ) * xx [ 3 ] ;
xx [ 90 ] = xx [ 112 ] + xx [ 129 ] * xx [ 41 ] + xx [ 132 ] * xx [ 137 ] -
xx [ 135 ] * xx [ 131 ] ; xx [ 112 ] = xx [ 116 ] + xx [ 111 ] + xx [ 135 ] ;
xx [ 116 ] = xx [ 94 ] * xx [ 112 ] ; xx [ 134 ] = xx [ 116 ] * xx [ 94 ] ;
xx [ 136 ] = xx [ 120 ] + xx [ 132 ] ; xx [ 120 ] = xx [ 136 ] * xx [ 94 ] ;
xx [ 138 ] = xx [ 4 ] * xx [ 120 ] ; xx [ 139 ] = xx [ 7 ] * xx [ 94 ] ; xx [
140 ] = xx [ 4 ] * xx [ 139 ] ; xx [ 141 ] = xx [ 3 ] * ( xx [ 134 ] - xx [
138 ] ) - xx [ 112 ] - xx [ 3 ] * xx [ 140 ] ; xx [ 142 ] = ( xx [ 4 ] * xx [
116 ] + xx [ 120 ] * xx [ 94 ] ) * xx [ 3 ] ; xx [ 116 ] = xx [ 3 ] * xx [
139 ] * xx [ 94 ] ; xx [ 120 ] = xx [ 136 ] - ( xx [ 142 ] + xx [ 116 ] ) +
xx [ 7 ] ; xx [ 139 ] = xx [ 120 ] * xx [ 97 ] ; xx [ 143 ] = xx [ 97 ] * xx
[ 141 ] ; xx [ 144 ] = xx [ 99 ] * xx [ 86 ] + xx [ 104 ] * xx [ 93 ] ; xx [
145 ] = xx [ 99 ] * xx [ 104 ] - xx [ 93 ] * xx [ 86 ] ; xx [ 146 ] = xx [
144 ] * xx [ 96 ] + xx [ 98 ] * xx [ 145 ] ; xx [ 147 ] = xx [ 146 ] * xx [
23 ] ; xx [ 148 ] = xx [ 141 ] - ( xx [ 95 ] * xx [ 139 ] + xx [ 143 ] * xx [
97 ] ) * xx [ 3 ] - xx [ 3 ] * xx [ 147 ] * ( xx [ 96 ] * xx [ 145 ] - xx [
98 ] * xx [ 144 ] ) ; xx [ 141 ] = xx [ 99 ] * xx [ 93 ] ; xx [ 149 ] = xx [
3 ] * xx [ 141 ] ; xx [ 150 ] = xx [ 29 ] + xx [ 123 ] ; xx [ 123 ] = xx [
126 ] / xx [ 61 ] ; xx [ 151 ] = xx [ 150 ] - xx [ 126 ] * xx [ 123 ] ; xx [
152 ] = xx [ 99 ] * xx [ 99 ] ; xx [ 153 ] = xx [ 3 ] * xx [ 152 ] - xx [ 1 ]
; xx [ 154 ] = xx [ 123 ] * xx [ 133 ] ; xx [ 155 ] = xx [ 124 ] - xx [ 154 ]
; xx [ 156 ] = xx [ 149 ] * xx [ 151 ] + xx [ 153 ] * xx [ 155 ] ; xx [ 157 ]
= xx [ 35 ] - xx [ 154 ] ; xx [ 154 ] = xx [ 133 ] / xx [ 61 ] ; xx [ 158 ] =
xx [ 128 ] - xx [ 154 ] * xx [ 133 ] ; xx [ 159 ] = xx [ 149 ] * xx [ 157 ] +
xx [ 158 ] * xx [ 153 ] ; xx [ 160 ] = xx [ 156 ] * xx [ 149 ] + xx [ 159 ] *
xx [ 153 ] ; xx [ 161 ] = xx [ 29 ] + xx [ 160 ] ; xx [ 162 ] = ( xx [ 152 ]
+ xx [ 93 ] * xx [ 93 ] ) * xx [ 3 ] - xx [ 1 ] ; xx [ 152 ] = xx [ 127 ] -
xx [ 154 ] * xx [ 41 ] ; xx [ 163 ] = xx [ 113 ] + xx [ 123 ] * xx [ 41 ] ;
xx [ 164 ] = xx [ 162 ] * ( xx [ 153 ] * xx [ 152 ] - xx [ 163 ] * xx [ 149 ]
) ; xx [ 165 ] = xx [ 151 ] * xx [ 153 ] - xx [ 149 ] * xx [ 155 ] ; xx [ 151
] = xx [ 153 ] * xx [ 157 ] - xx [ 149 ] * xx [ 158 ] ; xx [ 155 ] = xx [ 149
] * xx [ 165 ] + xx [ 153 ] * xx [ 151 ] ; xx [ 157 ] = xx [ 160 ] * xx [ 132
] - xx [ 155 ] * xx [ 135 ] ; xx [ 158 ] = xx [ 164 ] - xx [ 157 ] ; xx [ 160
] = xx [ 161 ] * xx [ 7 ] - xx [ 158 ] ; xx [ 166 ] = xx [ 41 ] / xx [ 61 ] ;
xx [ 167 ] = ( xx [ 163 ] * xx [ 153 ] + xx [ 149 ] * xx [ 152 ] ) * xx [ 162
] ; xx [ 152 ] = xx [ 132 ] * xx [ 164 ] + xx [ 167 ] * xx [ 135 ] ; xx [ 163
] = xx [ 156 ] * xx [ 153 ] - xx [ 159 ] * xx [ 149 ] ; xx [ 156 ] = xx [ 153
] * xx [ 165 ] - xx [ 149 ] * xx [ 151 ] ; xx [ 151 ] = xx [ 132 ] * xx [ 163
] - xx [ 135 ] * xx [ 156 ] ; xx [ 159 ] = xx [ 26 ] + ( xx [ 130 ] - xx [
166 ] * xx [ 41 ] ) * xx [ 162 ] * xx [ 162 ] - xx [ 152 ] - xx [ 152 ] + xx
[ 157 ] * xx [ 132 ] - xx [ 151 ] * xx [ 135 ] ; xx [ 130 ] = xx [ 7 ] * xx [
158 ] - xx [ 159 ] ; xx [ 152 ] = xx [ 7 ] * xx [ 160 ] - xx [ 130 ] ; ii [ 0
] = factorSymmetricPosDef ( xx + 152 , 1 , xx + 157 ) ; if ( ii [ 0 ] != 0 )
{ return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint2' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 157 ] = ( xx [ 90 ] + xx [ 137 ] * xx [ 7 ] + xx [ 148
] ) / xx [ 152 ] ; xx [ 162 ] = xx [ 23 ] * xx [ 94 ] ; xx [ 164 ] = xx [ 44
] - xx [ 3 ] * xx [ 162 ] * xx [ 94 ] ; xx [ 165 ] = xx [ 137 ] - xx [ 157 ]
* xx [ 160 ] ; xx [ 137 ] = xx [ 167 ] + xx [ 151 ] ; xx [ 151 ] = xx [ 137 ]
+ xx [ 7 ] * xx [ 163 ] ; xx [ 167 ] = xx [ 131 ] - xx [ 151 ] * xx [ 157 ] ;
xx [ 131 ] = xx [ 94 ] * xx [ 167 ] ; xx [ 168 ] = xx [ 94 ] * xx [ 165 ] ;
xx [ 169 ] = xx [ 165 ] + xx [ 3 ] * ( xx [ 4 ] * xx [ 131 ] - xx [ 168 ] *
xx [ 94 ] ) ; xx [ 165 ] = xx [ 3 ] * xx [ 4 ] * xx [ 162 ] ; xx [ 162 ] = xx
[ 112 ] + xx [ 3 ] * ( xx [ 138 ] - xx [ 134 ] ) + xx [ 165 ] ; xx [ 112 ] =
xx [ 97 ] * xx [ 162 ] ; xx [ 134 ] = xx [ 136 ] - xx [ 142 ] + xx [ 164 ] ;
xx [ 136 ] = xx [ 134 ] * xx [ 97 ] ; xx [ 138 ] = xx [ 7 ] * xx [ 97 ] ; xx
[ 142 ] = xx [ 144 ] * xx [ 94 ] - xx [ 4 ] * xx [ 145 ] ; xx [ 170 ] = xx [
144 ] * xx [ 4 ] + xx [ 94 ] * xx [ 145 ] ; xx [ 144 ] = xx [ 97 ] * xx [ 142
] - xx [ 170 ] * xx [ 95 ] ; xx [ 145 ] = xx [ 23 ] * xx [ 144 ] ; xx [ 171 ]
= xx [ 3 ] * ( xx [ 112 ] * xx [ 97 ] - xx [ 95 ] * xx [ 136 ] ) - xx [ 162 ]
- xx [ 3 ] * xx [ 95 ] * xx [ 138 ] - xx [ 3 ] * ( xx [ 95 ] * xx [ 142 ] +
xx [ 170 ] * xx [ 97 ] ) * xx [ 145 ] ; xx [ 142 ] = xx [ 4 ] * xx [ 94 ] ;
xx [ 162 ] = xx [ 3 ] * xx [ 142 ] ; xx [ 170 ] = xx [ 29 ] + xx [ 156 ] ; xx
[ 156 ] = xx [ 151 ] / xx [ 152 ] ; xx [ 172 ] = xx [ 170 ] - xx [ 151 ] * xx
[ 156 ] ; xx [ 173 ] = xx [ 4 ] * xx [ 4 ] ; xx [ 174 ] = xx [ 3 ] * xx [ 173
] - xx [ 1 ] ; xx [ 175 ] = xx [ 156 ] * xx [ 160 ] ; xx [ 176 ] = xx [ 163 ]
- xx [ 175 ] ; xx [ 177 ] = xx [ 162 ] * xx [ 172 ] + xx [ 174 ] * xx [ 176 ]
; xx [ 178 ] = xx [ 155 ] - xx [ 175 ] ; xx [ 175 ] = xx [ 160 ] / xx [ 152 ]
; xx [ 179 ] = xx [ 161 ] - xx [ 175 ] * xx [ 160 ] ; xx [ 180 ] = xx [ 162 ]
* xx [ 178 ] + xx [ 179 ] * xx [ 174 ] ; xx [ 181 ] = xx [ 177 ] * xx [ 162 ]
+ xx [ 180 ] * xx [ 174 ] ; xx [ 182 ] = ( xx [ 173 ] + xx [ 94 ] * xx [ 94 ]
) * xx [ 3 ] - xx [ 1 ] ; xx [ 173 ] = xx [ 158 ] - xx [ 175 ] * xx [ 130 ] ;
xx [ 183 ] = xx [ 137 ] + xx [ 156 ] * xx [ 130 ] ; xx [ 184 ] = xx [ 182 ] *
( xx [ 174 ] * xx [ 173 ] - xx [ 183 ] * xx [ 162 ] ) ; xx [ 185 ] = xx [ 172
] * xx [ 174 ] - xx [ 162 ] * xx [ 176 ] ; xx [ 172 ] = xx [ 178 ] * xx [ 174
] - xx [ 162 ] * xx [ 179 ] ; xx [ 176 ] = xx [ 162 ] * xx [ 185 ] + xx [ 174
] * xx [ 172 ] ; xx [ 178 ] = xx [ 181 ] * xx [ 164 ] - xx [ 176 ] * xx [ 165
] ; xx [ 179 ] = xx [ 184 ] - xx [ 178 ] ; xx [ 186 ] = xx [ 130 ] / xx [ 152
] ; xx [ 187 ] = ( xx [ 183 ] * xx [ 174 ] + xx [ 162 ] * xx [ 173 ] ) * xx [
182 ] ; xx [ 173 ] = xx [ 164 ] * xx [ 184 ] + xx [ 187 ] * xx [ 165 ] ; xx [
183 ] = xx [ 164 ] * ( xx [ 177 ] * xx [ 174 ] - xx [ 180 ] * xx [ 162 ] ) -
xx [ 165 ] * ( xx [ 174 ] * xx [ 185 ] - xx [ 162 ] * xx [ 172 ] ) ; xx [ 172
] = xx [ 7 ] * ( ( xx [ 29 ] + xx [ 181 ] ) * xx [ 7 ] - xx [ 179 ] ) - ( xx
[ 7 ] * xx [ 179 ] - ( ( xx [ 159 ] - xx [ 186 ] * xx [ 130 ] ) * xx [ 182 ]
* xx [ 182 ] - xx [ 173 ] - xx [ 173 ] + xx [ 178 ] * xx [ 164 ] - xx [ 183 ]
* xx [ 165 ] ) ) + xx [ 26 ] ; ii [ 0 ] = factorSymmetricPosDef ( xx + 172 ,
1 , xx + 159 ) ; if ( ii [ 0 ] != 0 ) { return sm_ssci_recordRunTimeError (
"physmod:sm:core:compiler:mechanism:mechanism:degenerateMassFoll" ,
 "'MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint1' has a degenerate mass distribution on its follower side."
, neDiagMgr ) ; } xx [ 159 ] = ( xx [ 90 ] + xx [ 157 ] * xx [ 130 ] + xx [
164 ] * xx [ 169 ] - xx [ 165 ] * ( xx [ 167 ] - ( xx [ 4 ] * xx [ 168 ] + xx
[ 131 ] * xx [ 94 ] ) * xx [ 3 ] ) + xx [ 169 ] * xx [ 7 ] + xx [ 171 ] ) /
xx [ 172 ] ; xx [ 90 ] = xx [ 159 ] * xx [ 165 ] ; xx [ 131 ] = xx [ 90 ] *
xx [ 94 ] ; xx [ 167 ] = xx [ 7 ] * xx [ 159 ] + xx [ 164 ] * xx [ 159 ] ; xx
[ 168 ] = xx [ 167 ] * xx [ 94 ] ; xx [ 169 ] = xx [ 90 ] - xx [ 3 ] * ( xx [
131 ] * xx [ 94 ] + xx [ 4 ] * xx [ 168 ] ) ; xx [ 90 ] = ( xx [ 168 ] * xx [
94 ] - xx [ 4 ] * xx [ 131 ] ) * xx [ 3 ] - xx [ 167 ] ; xx [ 131 ] = xx [
157 ] + xx [ 159 ] * xx [ 186 ] + xx [ 156 ] * xx [ 169 ] + xx [ 175 ] * xx [
90 ] ; xx [ 157 ] = xx [ 159 ] + xx [ 131 ] ; xx [ 167 ] = xx [ 169 ] + xx [
157 ] * xx [ 135 ] ; xx [ 168 ] = xx [ 90 ] - xx [ 131 ] * xx [ 7 ] - xx [
157 ] * xx [ 132 ] ; xx [ 90 ] = xx [ 93 ] * xx [ 168 ] ; xx [ 169 ] = xx [
93 ] * xx [ 167 ] ; xx [ 173 ] = xx [ 167 ] + xx [ 3 ] * ( xx [ 99 ] * xx [
90 ] - xx [ 169 ] * xx [ 93 ] ) ; xx [ 167 ] = xx [ 168 ] - ( xx [ 99 ] * xx
[ 169 ] + xx [ 90 ] * xx [ 93 ] ) * xx [ 3 ] ; xx [ 90 ] = xx [ 129 ] + xx [
157 ] * xx [ 166 ] + xx [ 173 ] * xx [ 123 ] + xx [ 154 ] * xx [ 167 ] ; xx [
129 ] = xx [ 157 ] + xx [ 90 ] ; xx [ 157 ] = xx [ 167 ] - xx [ 90 ] * xx [ 7
] - xx [ 129 ] * xx [ 91 ] ; xx [ 167 ] = xx [ 89 ] - ( xx [ 3 ] * xx [ 87 ]
* xx [ 88 ] + ( xx [ 53 ] + xx [ 51 ] ) * xx [ 3 ] ) + xx [ 23 ] ; xx [ 51 ]
= xx [ 3 ] * xx [ 25 ] * xx [ 24 ] ; xx [ 25 ] = xx [ 17 ] + xx [ 3 ] * ( xx
[ 21 ] * xx [ 20 ] - xx [ 11 ] * xx [ 5 ] ) - xx [ 51 ] + xx [ 23 ] ; xx [ 5
] = xx [ 25 ] / xx [ 10 ] ; xx [ 21 ] = xx [ 2 ] * xx [ 5 ] ; xx [ 53 ] = xx
[ 21 ] * xx [ 8 ] ; xx [ 87 ] = xx [ 3 ] * xx [ 6 ] * xx [ 53 ] ; xx [ 88 ] =
xx [ 21 ] - xx [ 3 ] * xx [ 53 ] * xx [ 8 ] ; xx [ 21 ] = xx [ 26 ] * xx [ 5
] + xx [ 46 ] * xx [ 87 ] - xx [ 45 ] * xx [ 88 ] ; xx [ 53 ] = xx [ 54 ] - (
xx [ 51 ] + ( xx [ 59 ] * xx [ 16 ] + xx [ 19 ] * xx [ 58 ] ) * xx [ 3 ] ) +
xx [ 23 ] ; xx [ 51 ] = ( xx [ 53 ] - ( xx [ 21 ] + xx [ 7 ] * xx [ 88 ] ) )
/ xx [ 47 ] ; xx [ 54 ] = xx [ 87 ] + xx [ 38 ] * xx [ 51 ] ; xx [ 58 ] = xx
[ 88 ] + xx [ 50 ] * xx [ 51 ] ; xx [ 59 ] = xx [ 58 ] * xx [ 12 ] ; xx [ 87
] = xx [ 54 ] * xx [ 12 ] ; xx [ 88 ] = xx [ 54 ] + xx [ 3 ] * ( xx [ 15 ] *
xx [ 59 ] - xx [ 87 ] * xx [ 12 ] ) ; xx [ 54 ] = xx [ 58 ] - ( xx [ 15 ] *
xx [ 87 ] + xx [ 59 ] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 58 ] = ( xx [ 167 ] - (
xx [ 21 ] + xx [ 42 ] * xx [ 51 ] + xx [ 88 ] * xx [ 83 ] - xx [ 75 ] * xx [
54 ] + xx [ 7 ] * xx [ 54 ] ) ) / xx [ 74 ] ; xx [ 21 ] = xx [ 54 ] + xx [ 82
] * xx [ 58 ] ; xx [ 54 ] = ( xx [ 21 ] - ( xx [ 18 ] * xx [ 18 ] * xx [ 21 ]
- ( xx [ 88 ] + xx [ 62 ] * xx [ 58 ] ) * xx [ 18 ] * xx [ 14 ] ) * xx [ 3 ]
) / xx [ 72 ] ; xx [ 21 ] = xx [ 18 ] * xx [ 54 ] ; xx [ 59 ] = xx [ 3 ] * xx
[ 18 ] * xx [ 21 ] - xx [ 54 ] ; xx [ 87 ] = xx [ 3 ] * xx [ 21 ] * xx [ 14 ]
; xx [ 21 ] = xx [ 58 ] - ( xx [ 66 ] * xx [ 59 ] - xx [ 78 ] * xx [ 87 ] ) ;
xx [ 58 ] = xx [ 83 ] * xx [ 21 ] - xx [ 87 ] ; xx [ 87 ] = xx [ 59 ] + xx [
7 ] * xx [ 21 ] - xx [ 75 ] * xx [ 21 ] ; xx [ 59 ] = xx [ 12 ] * xx [ 87 ] ;
xx [ 88 ] = xx [ 58 ] * xx [ 12 ] ; xx [ 89 ] = xx [ 58 ] - ( xx [ 15 ] * xx
[ 59 ] + xx [ 88 ] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 58 ] = xx [ 87 ] + xx [ 3
] * ( xx [ 15 ] * xx [ 88 ] - xx [ 59 ] * xx [ 12 ] ) ; xx [ 59 ] = xx [ 51 ]
- ( xx [ 56 ] * xx [ 21 ] + xx [ 69 ] * xx [ 89 ] + xx [ 58 ] * xx [ 63 ] ) ;
xx [ 51 ] = xx [ 21 ] + xx [ 59 ] ; xx [ 87 ] = xx [ 58 ] + xx [ 7 ] * xx [
59 ] - xx [ 45 ] * xx [ 51 ] ; xx [ 58 ] = xx [ 5 ] - ( xx [ 51 ] * xx [ 79 ]
+ ( xx [ 87 ] + xx [ 3 ] * ( xx [ 6 ] * ( xx [ 89 ] + xx [ 51 ] * xx [ 46 ] )
* xx [ 8 ] - xx [ 8 ] * xx [ 87 ] * xx [ 8 ] ) ) * xx [ 43 ] ) ; xx [ 5 ] =
xx [ 3 ] * xx [ 108 ] * xx [ 107 ] ; xx [ 51 ] = xx [ 101 ] + xx [ 3 ] * ( xx
[ 105 ] * xx [ 103 ] - xx [ 92 ] * xx [ 102 ] ) - xx [ 5 ] + xx [ 23 ] ; xx [
87 ] = xx [ 51 ] / xx [ 10 ] ; xx [ 88 ] = xx [ 2 ] * xx [ 87 ] ; xx [ 89 ] =
xx [ 88 ] * xx [ 86 ] ; xx [ 102 ] = xx [ 3 ] * xx [ 89 ] * xx [ 86 ] - xx [
88 ] ; xx [ 88 ] = xx [ 3 ] * xx [ 104 ] * xx [ 89 ] ; xx [ 89 ] = xx [ 91 ]
* xx [ 102 ] - xx [ 111 ] * xx [ 88 ] - xx [ 26 ] * xx [ 87 ] ; xx [ 105 ] =
xx [ 115 ] - xx [ 3 ] * ( xx [ 121 ] * xx [ 96 ] + xx [ 98 ] * xx [ 117 ] ) -
xx [ 5 ] + xx [ 23 ] ; xx [ 5 ] = ( xx [ 89 ] + xx [ 7 ] * xx [ 102 ] + xx [
105 ] ) / xx [ 61 ] ; xx [ 108 ] = xx [ 102 ] - xx [ 5 ] * xx [ 133 ] ; xx [
102 ] = xx [ 88 ] - xx [ 126 ] * xx [ 5 ] ; xx [ 88 ] = xx [ 93 ] * xx [ 102
] ; xx [ 115 ] = xx [ 93 ] * xx [ 108 ] ; xx [ 117 ] = xx [ 108 ] + xx [ 3 ]
* ( xx [ 99 ] * xx [ 88 ] - xx [ 115 ] * xx [ 93 ] ) ; xx [ 108 ] = xx [ 102
] - ( xx [ 99 ] * xx [ 115 ] + xx [ 88 ] * xx [ 93 ] ) * xx [ 3 ] ; xx [ 88 ]
= xx [ 89 ] + xx [ 5 ] * xx [ 41 ] + xx [ 132 ] * xx [ 117 ] - xx [ 135 ] *
xx [ 108 ] ; xx [ 89 ] = xx [ 120 ] + xx [ 3 ] * ( xx [ 95 ] * xx [ 143 ] -
xx [ 139 ] * xx [ 97 ] ) - xx [ 3 ] * xx [ 146 ] * xx [ 147 ] + xx [ 23 ] ;
xx [ 102 ] = ( xx [ 88 ] + xx [ 117 ] * xx [ 7 ] + xx [ 89 ] ) / xx [ 152 ] ;
xx [ 115 ] = xx [ 117 ] - xx [ 102 ] * xx [ 160 ] ; xx [ 117 ] = xx [ 108 ] -
xx [ 151 ] * xx [ 102 ] ; xx [ 108 ] = xx [ 94 ] * xx [ 117 ] ; xx [ 120 ] =
xx [ 94 ] * xx [ 115 ] ; xx [ 121 ] = xx [ 115 ] + xx [ 3 ] * ( xx [ 4 ] * xx
[ 108 ] - xx [ 120 ] * xx [ 94 ] ) ; xx [ 115 ] = xx [ 134 ] - ( ( xx [ 95 ]
* xx [ 112 ] + xx [ 136 ] * xx [ 97 ] ) * xx [ 3 ] + xx [ 3 ] * xx [ 138 ] *
xx [ 97 ] ) - xx [ 3 ] * xx [ 145 ] * xx [ 144 ] + xx [ 44 ] ; xx [ 44 ] = (
xx [ 88 ] + xx [ 102 ] * xx [ 130 ] + xx [ 164 ] * xx [ 121 ] - xx [ 165 ] *
( xx [ 117 ] - ( xx [ 4 ] * xx [ 120 ] + xx [ 108 ] * xx [ 94 ] ) * xx [ 3 ]
) + xx [ 121 ] * xx [ 7 ] + xx [ 115 ] ) / xx [ 172 ] ; xx [ 88 ] = xx [ 44 ]
* xx [ 165 ] ; xx [ 108 ] = xx [ 88 ] * xx [ 94 ] ; xx [ 112 ] = xx [ 7 ] *
xx [ 44 ] + xx [ 164 ] * xx [ 44 ] ; xx [ 117 ] = xx [ 112 ] * xx [ 94 ] ; xx
[ 120 ] = xx [ 88 ] - xx [ 3 ] * ( xx [ 108 ] * xx [ 94 ] + xx [ 4 ] * xx [
117 ] ) ; xx [ 88 ] = ( xx [ 117 ] * xx [ 94 ] - xx [ 4 ] * xx [ 108 ] ) * xx
[ 3 ] - xx [ 112 ] ; xx [ 108 ] = xx [ 102 ] + xx [ 44 ] * xx [ 186 ] + xx [
156 ] * xx [ 120 ] + xx [ 175 ] * xx [ 88 ] ; xx [ 102 ] = xx [ 44 ] + xx [
108 ] ; xx [ 112 ] = xx [ 120 ] + xx [ 102 ] * xx [ 135 ] ; xx [ 117 ] = xx [
88 ] - xx [ 108 ] * xx [ 7 ] - xx [ 102 ] * xx [ 132 ] ; xx [ 88 ] = xx [ 93
] * xx [ 117 ] ; xx [ 120 ] = xx [ 93 ] * xx [ 112 ] ; xx [ 121 ] = xx [ 112
] + xx [ 3 ] * ( xx [ 99 ] * xx [ 88 ] - xx [ 120 ] * xx [ 93 ] ) ; xx [ 112
] = xx [ 117 ] - ( xx [ 99 ] * xx [ 120 ] + xx [ 88 ] * xx [ 93 ] ) * xx [ 3
] ; xx [ 88 ] = xx [ 5 ] + xx [ 102 ] * xx [ 166 ] + xx [ 121 ] * xx [ 123 ]
+ xx [ 154 ] * xx [ 112 ] ; xx [ 5 ] = xx [ 102 ] + xx [ 88 ] ; xx [ 102 ] =
xx [ 112 ] - xx [ 88 ] * xx [ 7 ] - xx [ 5 ] * xx [ 91 ] ; xx [ 112 ] = xx [
87 ] + xx [ 43 ] * ( xx [ 102 ] - ( xx [ 104 ] * xx [ 86 ] * ( xx [ 121 ] +
xx [ 5 ] * xx [ 111 ] ) + xx [ 86 ] * xx [ 102 ] * xx [ 86 ] ) * xx [ 3 ] ) -
xx [ 5 ] * xx [ 79 ] ; xx [ 5 ] = xx [ 57 ] * xx [ 21 ] - xx [ 54 ] + xx [ 60
] * xx [ 59 ] + xx [ 58 ] * xx [ 27 ] + xx [ 44 ] * xx [ 171 ] + xx [ 108 ] *
xx [ 148 ] - xx [ 88 ] * xx [ 122 ] + xx [ 112 ] * xx [ 110 ] ; xx [ 143 ] =
xx [ 70 ] + xx [ 57 ] * xx [ 30 ] + xx [ 60 ] * xx [ 76 ] + ( xx [ 28 ] - (
xx [ 49 ] * xx [ 79 ] + ( xx [ 81 ] + xx [ 3 ] * ( xx [ 6 ] * ( xx [ 84 ] +
xx [ 49 ] * xx [ 46 ] ) * xx [ 8 ] - xx [ 8 ] * xx [ 81 ] * xx [ 8 ] ) ) * xx
[ 43 ] ) ) * xx [ 27 ] + xx [ 159 ] * xx [ 171 ] + xx [ 131 ] * xx [ 148 ] -
xx [ 90 ] * xx [ 122 ] + ( xx [ 100 ] + xx [ 43 ] * ( xx [ 157 ] - ( xx [ 104
] * xx [ 86 ] * ( xx [ 173 ] + xx [ 129 ] * xx [ 111 ] ) + xx [ 86 ] * xx [
157 ] * xx [ 86 ] ) * xx [ 3 ] ) - xx [ 129 ] * xx [ 79 ] ) * xx [ 110 ] ; xx
[ 144 ] = xx [ 5 ] ; xx [ 145 ] = xx [ 5 ] ; xx [ 146 ] = xx [ 167 ] * xx [
21 ] + xx [ 53 ] * xx [ 59 ] + xx [ 25 ] * xx [ 58 ] + xx [ 115 ] * xx [ 44 ]
+ xx [ 108 ] * xx [ 89 ] + xx [ 88 ] * xx [ 105 ] + xx [ 112 ] * xx [ 51 ] ;
xx [ 5 ] = state [ 3 ] + state [ 5 ] ; xx [ 21 ] = xx [ 5 ] * xx [ 5 ] * xx [
46 ] ; xx [ 28 ] = xx [ 21 ] * xx [ 8 ] ; xx [ 30 ] = xx [ 5 ] * xx [ 45 ] *
xx [ 5 ] ; xx [ 44 ] = xx [ 30 ] * xx [ 8 ] ; xx [ 49 ] = xx [ 29 ] * ( xx [
3 ] * ( xx [ 28 ] * xx [ 8 ] - xx [ 6 ] * xx [ 44 ] ) - xx [ 21 ] ) - input [
26 ] ; xx [ 21 ] = 4.010704565915762e-6 ; xx [ 54 ] = 4.010704565915763e-7 ;
xx [ 58 ] = xx [ 5 ] + state [ 7 ] ; xx [ 59 ] = state [ 9 ] + state [ 11 ] ;
xx [ 70 ] = xx [ 59 ] + state [ 13 ] ; xx [ 76 ] = xx [ 70 ] + state [ 15 ] ;
xx [ 81 ] = xx [ 58 ] + xx [ 76 ] ; xx [ 84 ] = xx [ 106 ] * xx [ 24 ] - xx [
22 ] * xx [ 107 ] ; xx [ 87 ] = xx [ 22 ] * xx [ 106 ] + xx [ 24 ] * xx [ 107
] ; xx [ 88 ] = state [ 16 ] + pm_math_canonicalAngle ( xx [ 3 ] * atan2 (
sqrt ( xx [ 84 ] * xx [ 84 ] ) , fabs ( - xx [ 87 ] ) ) * ( ( xx [ 87 ] * xx
[ 84 ] ) < 0.0 ? - 1.0 : + 1.0 ) - state [ 16 ] ) ; xx [ 84 ] = xx [ 81 ] *
xx [ 54 ] - xx [ 88 ] * xx [ 21 ] ; xx [ 87 ] = xx [ 84 ] - input [ 28 ] ; xx
[ 90 ] = xx [ 21 ] * state [ 6 ] + xx [ 54 ] * state [ 7 ] + xx [ 87 ] + xx [
7 ] * xx [ 49 ] ; xx [ 100 ] = xx [ 90 ] / xx [ 10 ] ; xx [ 102 ] = xx [ 49 ]
- xx [ 2 ] * xx [ 100 ] ; xx [ 108 ] = ( xx [ 7 ] * ( xx [ 5 ] + xx [ 58 ] )
* state [ 7 ] + ( xx [ 6 ] * xx [ 28 ] + xx [ 44 ] * xx [ 8 ] ) * xx [ 3 ] -
xx [ 30 ] ) * xx [ 29 ] - input [ 24 ] ; xx [ 28 ] = xx [ 8 ] * xx [ 108 ] ;
xx [ 30 ] = xx [ 6 ] * xx [ 28 ] ; xx [ 44 ] = xx [ 8 ] * xx [ 102 ] ; xx [
58 ] = xx [ 102 ] - ( xx [ 30 ] + xx [ 44 ] * xx [ 8 ] ) * xx [ 3 ] ; xx [
102 ] = xx [ 83 ] * state [ 3 ] * state [ 3 ] ; xx [ 112 ] = xx [ 102 ] * xx
[ 12 ] ; xx [ 117 ] = xx [ 75 ] * state [ 3 ] * state [ 3 ] ; xx [ 120 ] = xx
[ 117 ] * xx [ 12 ] ; xx [ 121 ] = xx [ 7 ] * ( state [ 3 ] + xx [ 5 ] ) *
state [ 5 ] + ( xx [ 15 ] * xx [ 112 ] + xx [ 120 ] * xx [ 12 ] ) * xx [ 3 ]
- xx [ 117 ] ; xx [ 117 ] = xx [ 3 ] * ( xx [ 112 ] * xx [ 12 ] - xx [ 15 ] *
xx [ 120 ] ) - xx [ 102 ] ; xx [ 102 ] = xx [ 121 ] * xx [ 48 ] + xx [ 40 ] *
xx [ 117 ] ; xx [ 40 ] = xx [ 58 ] - input [ 32 ] + xx [ 102 ] ; xx [ 48 ] =
xx [ 21 ] * state [ 4 ] + xx [ 54 ] * state [ 5 ] ; xx [ 112 ] = xx [ 28 ] *
xx [ 8 ] ; xx [ 28 ] = xx [ 108 ] + xx [ 3 ] * ( xx [ 6 ] * xx [ 44 ] - xx [
112 ] ) ; xx [ 44 ] = xx [ 121 ] * xx [ 31 ] + xx [ 117 ] * xx [ 39 ] ; xx [
31 ] = xx [ 87 ] - xx [ 26 ] * xx [ 100 ] + xx [ 28 ] * xx [ 46 ] - xx [ 45 ]
* xx [ 58 ] - input [ 34 ] - xx [ 44 ] ; xx [ 39 ] = ( xx [ 48 ] + xx [ 31 ]
+ xx [ 40 ] * xx [ 7 ] ) / xx [ 47 ] ; xx [ 58 ] = xx [ 40 ] - xx [ 50 ] * xx
[ 39 ] ; xx [ 40 ] = xx [ 77 ] * xx [ 121 ] + xx [ 117 ] * xx [ 36 ] ; xx [
36 ] = xx [ 28 ] - input [ 30 ] + xx [ 40 ] - xx [ 38 ] * xx [ 39 ] ; xx [ 28
] = xx [ 12 ] * xx [ 36 ] ; xx [ 77 ] = xx [ 12 ] * xx [ 58 ] ; xx [ 120 ] =
xx [ 58 ] - ( xx [ 15 ] * xx [ 28 ] + xx [ 77 ] * xx [ 12 ] ) * xx [ 3 ] ; xx
[ 58 ] = state [ 3 ] * state [ 3 ] ; xx [ 129 ] = xx [ 7 ] * xx [ 58 ] ; xx [
131 ] = xx [ 129 ] * xx [ 34 ] ; xx [ 34 ] = xx [ 120 ] - input [ 38 ] + xx [
131 ] ; xx [ 134 ] = xx [ 21 ] * state [ 2 ] + xx [ 54 ] * state [ 3 ] ; xx [
136 ] = xx [ 36 ] + xx [ 3 ] * ( xx [ 15 ] * xx [ 77 ] - xx [ 28 ] * xx [ 12
] ) ; xx [ 28 ] = xx [ 129 ] * xx [ 64 ] ; xx [ 36 ] = ( xx [ 134 ] + xx [ 31
] - xx [ 42 ] * xx [ 39 ] + xx [ 136 ] * xx [ 83 ] - xx [ 75 ] * xx [ 120 ] -
input [ 40 ] - xx [ 28 ] + xx [ 34 ] * xx [ 7 ] ) / xx [ 74 ] ; xx [ 31 ] =
xx [ 34 ] - xx [ 82 ] * xx [ 36 ] ; xx [ 34 ] = xx [ 80 ] * xx [ 129 ] ; xx [
64 ] = ( xx [ 31 ] - ( xx [ 18 ] * xx [ 18 ] * xx [ 31 ] - xx [ 18 ] * ( xx [
136 ] - input [ 36 ] + xx [ 34 ] - xx [ 62 ] * xx [ 36 ] ) * xx [ 14 ] ) * xx
[ 3 ] ) / xx [ 72 ] ; xx [ 31 ] = 1.0 ; xx [ 77 ] = xx [ 65 ] * state [ 3 ] *
state [ 3 ] ; xx [ 65 ] = xx [ 31 ] * xx [ 77 ] ; xx [ 80 ] =
2.220446049250313e-16 ; xx [ 120 ] = xx [ 3 ] * xx [ 73 ] * state [ 3 ] *
state [ 3 ] ; xx [ 73 ] = xx [ 80 ] * xx [ 120 ] ; xx [ 136 ] = xx [ 65 ] -
xx [ 73 ] ; xx [ 138 ] = xx [ 31 ] * xx [ 120 ] ; xx [ 31 ] = xx [ 80 ] * xx
[ 77 ] ; xx [ 77 ] = xx [ 138 ] + xx [ 31 ] ; xx [ 80 ] = xx [ 3 ] * xx [ 19
] * xx [ 16 ] ; xx [ 120 ] = xx [ 3 ] * xx [ 67 ] * state [ 5 ] * state [ 5 ]
; xx [ 67 ] = xx [ 71 ] * state [ 5 ] * state [ 5 ] ; xx [ 139 ] = xx [ 3 ] *
xx [ 16 ] * xx [ 16 ] - xx [ 1 ] ; xx [ 147 ] = xx [ 80 ] * xx [ 120 ] - xx [
67 ] * xx [ 139 ] ; xx [ 157 ] = xx [ 3 ] * xx [ 20 ] * xx [ 20 ] - xx [ 1 ]
; xx [ 159 ] = xx [ 3 ] * xx [ 157 ] * state [ 3 ] * state [ 5 ] ; xx [ 168 ]
= xx [ 136 ] * xx [ 68 ] - xx [ 77 ] * xx [ 71 ] + xx [ 147 ] - xx [ 159 ] ;
xx [ 169 ] = 4.0 ; xx [ 173 ] = xx [ 11 ] * xx [ 20 ] ; xx [ 177 ] = xx [ 169
] * xx [ 173 ] * state [ 3 ] * state [ 5 ] ; xx [ 178 ] = xx [ 120 ] * xx [
139 ] ; xx [ 120 ] = xx [ 80 ] * xx [ 67 ] ; xx [ 67 ] = xx [ 177 ] - ( xx [
68 ] * xx [ 77 ] + xx [ 136 ] * xx [ 71 ] + xx [ 178 ] + xx [ 120 ] ) ; xx [
80 ] = 1.0e-3 ; xx [ 139 ] = xx [ 18 ] * xx [ 80 ] ; xx [ 179 ] = xx [ 58 ] *
( xx [ 80 ] - xx [ 3 ] * xx [ 18 ] * xx [ 139 ] ) ; xx [ 180 ] = xx [ 13 ] *
xx [ 13 ] * xx [ 179 ] ; xx [ 181 ] = xx [ 3 ] * xx [ 139 ] * xx [ 14 ] * xx
[ 58 ] ; xx [ 58 ] = xx [ 13 ] * xx [ 13 ] * xx [ 181 ] ; xx [ 13 ] = state [
5 ] * state [ 5 ] ; xx [ 139 ] = xx [ 80 ] * xx [ 12 ] ; xx [ 182 ] = xx [ 13
] * ( xx [ 80 ] - xx [ 3 ] * xx [ 139 ] * xx [ 12 ] ) ; xx [ 184 ] = xx [ 3 ]
* xx [ 15 ] * xx [ 139 ] * xx [ 13 ] ; xx [ 13 ] = xx [ 19 ] * xx [ 184 ] ;
xx [ 139 ] = xx [ 19 ] * xx [ 182 ] ; xx [ 185 ] = ( xx [ 7 ] - xx [ 52 ] ) *
state [ 5 ] * state [ 3 ] ; xx [ 52 ] = xx [ 3 ] * xx [ 55 ] * state [ 5 ] *
state [ 3 ] ; xx [ 55 ] = xx [ 19 ] * xx [ 52 ] ; xx [ 188 ] = xx [ 19 ] * xx
[ 185 ] ; xx [ 189 ] = state [ 7 ] * state [ 7 ] ; xx [ 190 ] = xx [ 80 ] *
xx [ 8 ] ; xx [ 191 ] = xx [ 189 ] * ( xx [ 80 ] - xx [ 3 ] * xx [ 190 ] * xx
[ 8 ] ) ; xx [ 192 ] = xx [ 3 ] * xx [ 6 ] * xx [ 190 ] * xx [ 189 ] ; xx [
189 ] = xx [ 11 ] * xx [ 192 ] ; xx [ 190 ] = xx [ 11 ] * xx [ 191 ] ; xx [
193 ] = xx [ 5 ] * xx [ 17 ] * state [ 7 ] ; xx [ 17 ] = xx [ 5 ] * xx [ 3 ]
* xx [ 9 ] * state [ 7 ] ; xx [ 9 ] = xx [ 11 ] * xx [ 17 ] ; xx [ 194 ] = xx
[ 11 ] * xx [ 193 ] ; xx [ 195 ] = xx [ 37 ] * state [ 7 ] * state [ 7 ] ; xx
[ 196 ] = xx [ 3 ] * xx [ 173 ] ; xx [ 173 ] = xx [ 3 ] * xx [ 32 ] * state [
7 ] * state [ 7 ] ; xx [ 32 ] = xx [ 95 ] * xx [ 97 ] ; xx [ 197 ] = xx [ 3 ]
* xx [ 32 ] * state [ 9 ] * state [ 9 ] ; xx [ 198 ] = xx [ 3 ] * xx [ 95 ] *
xx [ 95 ] - xx [ 1 ] ; xx [ 199 ] = xx [ 198 ] * state [ 9 ] * state [ 9 ] ;
xx [ 200 ] = xx [ 3 ] * xx [ 32 ] ; xx [ 32 ] = xx [ 3 ] * xx [ 142 ] * state
[ 11 ] * state [ 11 ] ; xx [ 142 ] = xx [ 174 ] * state [ 11 ] * state [ 11 ]
; xx [ 201 ] = xx [ 3 ] * xx [ 96 ] * xx [ 96 ] - xx [ 1 ] ; xx [ 202 ] = xx
[ 162 ] * xx [ 197 ] - xx [ 199 ] * xx [ 174 ] + xx [ 200 ] * xx [ 32 ] - xx
[ 142 ] * xx [ 198 ] - xx [ 3 ] * xx [ 201 ] * state [ 9 ] * state [ 11 ] ;
xx [ 203 ] = xx [ 202 ] * xx [ 153 ] ; xx [ 204 ] = xx [ 162 ] * xx [ 199 ] ;
xx [ 162 ] = xx [ 197 ] * xx [ 174 ] ; xx [ 174 ] = xx [ 32 ] * xx [ 198 ] ;
xx [ 32 ] = xx [ 200 ] * xx [ 142 ] ; xx [ 142 ] = xx [ 98 ] * xx [ 96 ] ; xx
[ 198 ] = xx [ 169 ] * xx [ 142 ] * state [ 9 ] * state [ 11 ] ; xx [ 169 ] =
xx [ 204 ] + xx [ 162 ] + xx [ 174 ] + xx [ 32 ] + xx [ 198 ] ; xx [ 200 ] =
xx [ 3 ] * xx [ 142 ] ; xx [ 142 ] = xx [ 3 ] * xx [ 141 ] * state [ 13 ] *
state [ 13 ] ; xx [ 141 ] = xx [ 153 ] * state [ 13 ] * state [ 13 ] ; xx [
205 ] = xx [ 200 ] * xx [ 142 ] - xx [ 141 ] * xx [ 201 ] ; xx [ 206 ] = xx [
3 ] * xx [ 103 ] * xx [ 103 ] - xx [ 1 ] ; xx [ 207 ] = xx [ 3 ] * xx [ 59 ]
* xx [ 206 ] * state [ 13 ] ; xx [ 208 ] = xx [ 203 ] + xx [ 169 ] * xx [ 149
] + xx [ 205 ] - xx [ 207 ] ; xx [ 209 ] = xx [ 149 ] * xx [ 202 ] ; xx [ 210
] = xx [ 142 ] * xx [ 201 ] ; xx [ 142 ] = xx [ 200 ] * xx [ 141 ] ; xx [ 141
] = xx [ 3 ] * xx [ 92 ] * xx [ 103 ] ; xx [ 200 ] = xx [ 3 ] * xx [ 59 ] *
xx [ 141 ] * state [ 13 ] ; xx [ 201 ] = xx [ 169 ] * xx [ 153 ] - xx [ 209 ]
+ xx [ 210 ] + xx [ 142 ] - xx [ 200 ] ; xx [ 211 ] = xx [ 125 ] * state [ 15
] * state [ 15 ] ; xx [ 212 ] = xx [ 3 ] * xx [ 109 ] * state [ 15 ] * state
[ 15 ] ; xx [ 109 ] = state [ 9 ] * state [ 9 ] ; xx [ 213 ] = xx [ 80 ] * xx
[ 97 ] ; xx [ 214 ] = state [ 11 ] * state [ 11 ] ; xx [ 215 ] = xx [ 80 ] *
xx [ 94 ] ; xx [ 216 ] = xx [ 214 ] * ( xx [ 3 ] * xx [ 215 ] * xx [ 94 ] -
xx [ 80 ] ) ; xx [ 217 ] = xx [ 216 ] * xx [ 97 ] ; xx [ 218 ] = xx [ 3 ] *
xx [ 4 ] * xx [ 215 ] * xx [ 214 ] ; xx [ 214 ] = xx [ 218 ] * xx [ 97 ] ; xx
[ 215 ] = xx [ 3 ] * xx [ 140 ] * state [ 11 ] * state [ 9 ] ; xx [ 140 ] =
xx [ 215 ] * xx [ 97 ] ; xx [ 219 ] = ( xx [ 7 ] - xx [ 116 ] ) * state [ 11
] * state [ 9 ] ; xx [ 116 ] = xx [ 219 ] * xx [ 97 ] ; xx [ 220 ] = state [
13 ] * state [ 13 ] ; xx [ 221 ] = xx [ 80 ] * xx [ 93 ] ; xx [ 222 ] = xx [
220 ] * ( xx [ 3 ] * xx [ 221 ] * xx [ 93 ] - xx [ 80 ] ) ; xx [ 223 ] = xx [
98 ] * xx [ 222 ] ; xx [ 224 ] = xx [ 3 ] * xx [ 99 ] * xx [ 221 ] * xx [ 220
] ; xx [ 220 ] = xx [ 98 ] * xx [ 224 ] ; xx [ 221 ] = xx [ 59 ] * xx [ 3 ] *
xx [ 118 ] * state [ 13 ] ; xx [ 118 ] = xx [ 98 ] * xx [ 221 ] ; xx [ 225 ]
= xx [ 59 ] * ( xx [ 7 ] - xx [ 114 ] ) * state [ 13 ] ; xx [ 114 ] = xx [ 98
] * xx [ 225 ] ; xx [ 226 ] = state [ 15 ] * state [ 15 ] ; xx [ 227 ] = xx [
80 ] * xx [ 86 ] ; xx [ 228 ] = xx [ 226 ] * ( xx [ 3 ] * xx [ 227 ] * xx [
86 ] - xx [ 80 ] ) ; xx [ 229 ] = xx [ 3 ] * xx [ 104 ] * xx [ 227 ] * xx [
226 ] ; xx [ 226 ] = xx [ 92 ] * xx [ 229 ] ; xx [ 227 ] = xx [ 92 ] * xx [
228 ] ; xx [ 230 ] = xx [ 70 ] * xx [ 101 ] * state [ 15 ] ; xx [ 101 ] = xx
[ 92 ] * xx [ 230 ] ; xx [ 231 ] = xx [ 70 ] * xx [ 3 ] * xx [ 85 ] * state [
15 ] ; xx [ 85 ] = xx [ 92 ] * xx [ 231 ] ; xx [ 232 ] = xx [ 18 ] * xx [ 64
] ; xx [ 233 ] = xx [ 3 ] * xx [ 18 ] * xx [ 232 ] - xx [ 64 ] ; xx [ 234 ] =
xx [ 3 ] * xx [ 232 ] * xx [ 14 ] ; xx [ 232 ] = xx [ 36 ] + xx [ 66 ] * xx [
233 ] - xx [ 78 ] * xx [ 234 ] ; xx [ 36 ] = xx [ 129 ] - xx [ 234 ] - xx [
232 ] * xx [ 83 ] ; xx [ 234 ] = xx [ 233 ] - xx [ 232 ] * xx [ 7 ] + xx [
232 ] * xx [ 75 ] ; xx [ 233 ] = xx [ 234 ] * xx [ 12 ] ; xx [ 235 ] = xx [
12 ] * xx [ 36 ] ; xx [ 236 ] = xx [ 36 ] - ( xx [ 15 ] * xx [ 233 ] + xx [
235 ] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 36 ] = xx [ 234 ] + xx [ 3 ] * ( xx [
15 ] * xx [ 235 ] - xx [ 233 ] * xx [ 12 ] ) ; xx [ 233 ] = xx [ 39 ] + xx [
69 ] * xx [ 236 ] + xx [ 36 ] * xx [ 63 ] - xx [ 232 ] * xx [ 56 ] ; xx [ 39
] = xx [ 232 ] + xx [ 233 ] ; xx [ 234 ] = xx [ 36 ] - xx [ 233 ] * xx [ 7 ]
+ xx [ 117 ] + xx [ 39 ] * xx [ 45 ] ; xx [ 36 ] = xx [ 100 ] + ( xx [ 234 ]
+ xx [ 3 ] * ( xx [ 6 ] * xx [ 8 ] * ( xx [ 236 ] + xx [ 121 ] - xx [ 39 ] *
xx [ 46 ] ) - xx [ 234 ] * xx [ 8 ] * xx [ 8 ] ) ) * xx [ 43 ] - xx [ 39 ] *
xx [ 79 ] ; xx [ 39 ] = xx [ 21 ] * state [ 8 ] + xx [ 54 ] * state [ 9 ] ;
xx [ 100 ] = xx [ 70 ] * xx [ 70 ] * xx [ 91 ] ; xx [ 234 ] = xx [ 100 ] * xx
[ 86 ] ; xx [ 235 ] = xx [ 70 ] * xx [ 70 ] * xx [ 111 ] ; xx [ 236 ] = xx [
235 ] * xx [ 86 ] ; xx [ 237 ] = ( xx [ 3 ] * ( xx [ 104 ] * xx [ 234 ] + xx
[ 236 ] * xx [ 86 ] ) - xx [ 235 ] ) * xx [ 29 ] - input [ 20 ] ; xx [ 235 ]
= input [ 22 ] + xx [ 84 ] ; xx [ 84 ] = xx [ 21 ] * state [ 14 ] + xx [ 54 ]
* state [ 15 ] + xx [ 235 ] + xx [ 7 ] * xx [ 237 ] ; xx [ 238 ] = xx [ 84 ]
/ xx [ 10 ] ; xx [ 239 ] = xx [ 237 ] - xx [ 2 ] * xx [ 238 ] ; xx [ 240 ] =
xx [ 29 ] * ( ( xx [ 234 ] * xx [ 86 ] - xx [ 104 ] * xx [ 236 ] ) * xx [ 3 ]
- xx [ 100 ] - xx [ 7 ] * ( xx [ 70 ] + xx [ 76 ] ) * state [ 15 ] ) - input
[ 18 ] ; xx [ 29 ] = xx [ 86 ] * xx [ 240 ] ; xx [ 76 ] = xx [ 104 ] * xx [
29 ] ; xx [ 100 ] = xx [ 86 ] * xx [ 239 ] ; xx [ 234 ] = xx [ 239 ] + xx [ 3
] * ( xx [ 76 ] - xx [ 100 ] * xx [ 86 ] ) ; xx [ 236 ] = xx [ 59 ] * xx [
132 ] * xx [ 59 ] ; xx [ 239 ] = xx [ 236 ] * xx [ 93 ] ; xx [ 241 ] = xx [
59 ] * xx [ 59 ] * xx [ 135 ] ; xx [ 242 ] = xx [ 241 ] * xx [ 93 ] ; xx [
243 ] = ( xx [ 239 ] * xx [ 93 ] - xx [ 99 ] * xx [ 242 ] ) * xx [ 3 ] - xx [
236 ] - xx [ 7 ] * ( xx [ 59 ] + xx [ 70 ] ) * state [ 13 ] ; xx [ 236 ] = xx
[ 3 ] * ( xx [ 99 ] * xx [ 239 ] + xx [ 242 ] * xx [ 93 ] ) - xx [ 241 ] ; xx
[ 239 ] = xx [ 35 ] * xx [ 243 ] + xx [ 128 ] * xx [ 236 ] ; xx [ 35 ] = xx [
234 ] - input [ 14 ] + xx [ 239 ] ; xx [ 128 ] = xx [ 21 ] * state [ 12 ] +
xx [ 54 ] * state [ 13 ] ; xx [ 241 ] = xx [ 236 ] * xx [ 127 ] - xx [ 113 ]
* xx [ 243 ] ; xx [ 113 ] = xx [ 29 ] * xx [ 86 ] ; xx [ 29 ] = xx [ 240 ] -
( xx [ 104 ] * xx [ 100 ] + xx [ 113 ] ) * xx [ 3 ] ; xx [ 100 ] = xx [ 241 ]
- ( input [ 16 ] + xx [ 235 ] - xx [ 26 ] * xx [ 238 ] + xx [ 91 ] * xx [ 234
] - xx [ 111 ] * xx [ 29 ] ) ; xx [ 127 ] = ( xx [ 128 ] + xx [ 35 ] * xx [ 7
] - xx [ 100 ] ) / xx [ 61 ] ; xx [ 234 ] = xx [ 35 ] - xx [ 127 ] * xx [ 133
] ; xx [ 35 ] = xx [ 150 ] * xx [ 243 ] + xx [ 236 ] * xx [ 124 ] ; xx [ 124
] = xx [ 29 ] - input [ 12 ] + xx [ 35 ] - xx [ 126 ] * xx [ 127 ] ; xx [ 29
] = xx [ 93 ] * xx [ 124 ] ; xx [ 150 ] = xx [ 93 ] * xx [ 234 ] ; xx [ 242 ]
= xx [ 234 ] + xx [ 3 ] * ( xx [ 99 ] * xx [ 29 ] - xx [ 150 ] * xx [ 93 ] )
; xx [ 234 ] = xx [ 164 ] * state [ 9 ] * state [ 9 ] ; xx [ 244 ] = xx [ 234
] * xx [ 94 ] ; xx [ 245 ] = xx [ 165 ] * state [ 9 ] * state [ 9 ] ; xx [
246 ] = xx [ 245 ] * xx [ 94 ] ; xx [ 247 ] = ( xx [ 244 ] * xx [ 94 ] - xx [
4 ] * xx [ 246 ] ) * xx [ 3 ] - xx [ 234 ] - xx [ 7 ] * ( state [ 9 ] + xx [
59 ] ) * state [ 11 ] ; xx [ 59 ] = xx [ 3 ] * ( xx [ 4 ] * xx [ 244 ] + xx [
246 ] * xx [ 94 ] ) - xx [ 245 ] ; xx [ 234 ] = xx [ 155 ] * xx [ 247 ] + xx
[ 161 ] * xx [ 59 ] ; xx [ 155 ] = xx [ 242 ] - input [ 8 ] + xx [ 234 ] ; xx
[ 161 ] = xx [ 21 ] * state [ 10 ] + xx [ 54 ] * state [ 11 ] ; xx [ 21 ] =
xx [ 124 ] - ( xx [ 99 ] * xx [ 150 ] + xx [ 29 ] * xx [ 93 ] ) * xx [ 3 ] ;
xx [ 29 ] = xx [ 59 ] * xx [ 158 ] - xx [ 137 ] * xx [ 247 ] ; xx [ 54 ] = xx
[ 100 ] - xx [ 127 ] * xx [ 41 ] - ( xx [ 132 ] * xx [ 242 ] - xx [ 135 ] *
xx [ 21 ] ) - input [ 10 ] + xx [ 29 ] ; xx [ 100 ] = ( xx [ 161 ] + xx [ 155
] * xx [ 7 ] - xx [ 54 ] ) / xx [ 152 ] ; xx [ 124 ] = xx [ 155 ] - xx [ 100
] * xx [ 160 ] ; xx [ 137 ] = xx [ 170 ] * xx [ 247 ] + xx [ 59 ] * xx [ 163
] ; xx [ 150 ] = xx [ 21 ] - input [ 6 ] + xx [ 137 ] - xx [ 151 ] * xx [ 100
] ; xx [ 21 ] = xx [ 94 ] * xx [ 150 ] ; xx [ 155 ] = xx [ 94 ] * xx [ 124 ]
; xx [ 158 ] = xx [ 124 ] + xx [ 3 ] * ( xx [ 4 ] * xx [ 21 ] - xx [ 155 ] *
xx [ 94 ] ) ; xx [ 124 ] = xx [ 7 ] * xx [ 109 ] ; xx [ 163 ] = xx [ 176 ] *
xx [ 124 ] ; xx [ 170 ] = ( xx [ 187 ] + xx [ 183 ] ) * xx [ 124 ] ; xx [ 176
] = ( xx [ 39 ] + xx [ 7 ] * ( xx [ 158 ] - input [ 2 ] - xx [ 163 ] ) - ( xx
[ 54 ] - xx [ 100 ] * xx [ 130 ] - ( xx [ 164 ] * xx [ 158 ] - xx [ 165 ] * (
xx [ 150 ] - ( xx [ 4 ] * xx [ 155 ] + xx [ 21 ] * xx [ 94 ] ) * xx [ 3 ] ) )
- input [ 4 ] + xx [ 170 ] ) ) / xx [ 172 ] ; xx [ 21 ] = xx [ 124 ] - xx [
176 ] * xx [ 165 ] ; xx [ 54 ] = xx [ 21 ] * xx [ 94 ] ; xx [ 150 ] = xx [ 7
] * xx [ 176 ] + xx [ 164 ] * xx [ 176 ] ; xx [ 155 ] = xx [ 150 ] * xx [ 94
] ; xx [ 158 ] = xx [ 3 ] * ( xx [ 54 ] * xx [ 94 ] - xx [ 4 ] * xx [ 155 ] )
- xx [ 21 ] ; xx [ 21 ] = ( xx [ 4 ] * xx [ 54 ] + xx [ 155 ] * xx [ 94 ] ) *
xx [ 3 ] - xx [ 150 ] ; xx [ 54 ] = xx [ 100 ] + xx [ 176 ] * xx [ 186 ] + xx
[ 156 ] * xx [ 158 ] + xx [ 175 ] * xx [ 21 ] ; xx [ 100 ] = xx [ 176 ] + xx
[ 54 ] ; xx [ 150 ] = xx [ 158 ] + xx [ 247 ] + xx [ 100 ] * xx [ 135 ] ; xx
[ 155 ] = xx [ 21 ] - xx [ 54 ] * xx [ 7 ] + xx [ 59 ] - xx [ 100 ] * xx [
132 ] ; xx [ 21 ] = xx [ 93 ] * xx [ 155 ] ; xx [ 158 ] = xx [ 93 ] * xx [
150 ] ; xx [ 183 ] = xx [ 150 ] + xx [ 3 ] * ( xx [ 99 ] * xx [ 21 ] - xx [
158 ] * xx [ 93 ] ) ; xx [ 150 ] = xx [ 155 ] - ( xx [ 99 ] * xx [ 158 ] + xx
[ 21 ] * xx [ 93 ] ) * xx [ 3 ] ; xx [ 21 ] = xx [ 127 ] + xx [ 100 ] * xx [
166 ] + xx [ 183 ] * xx [ 123 ] + xx [ 154 ] * xx [ 150 ] ; xx [ 127 ] = xx [
100 ] + xx [ 21 ] ; xx [ 100 ] = xx [ 150 ] - xx [ 21 ] * xx [ 7 ] + xx [ 236
] - xx [ 127 ] * xx [ 91 ] ; xx [ 150 ] = xx [ 238 ] + xx [ 43 ] * ( xx [ 100
] - ( xx [ 104 ] * xx [ 86 ] * ( xx [ 183 ] + xx [ 243 ] + xx [ 127 ] * xx [
111 ] ) + xx [ 86 ] * xx [ 100 ] * xx [ 86 ] ) * xx [ 3 ] ) - xx [ 127 ] * xx
[ 79 ] ; xx [ 100 ] = xx [ 65 ] - xx [ 73 ] ; xx [ 65 ] = xx [ 138 ] + xx [
31 ] ; xx [ 31 ] = xx [ 100 ] * xx [ 71 ] + xx [ 65 ] * xx [ 68 ] + xx [ 120
] + xx [ 178 ] - xx [ 177 ] ; xx [ 73 ] = xx [ 68 ] * xx [ 100 ] - xx [ 65 ]
* xx [ 71 ] + xx [ 147 ] - xx [ 159 ] ; xx [ 68 ] = xx [ 162 ] + xx [ 204 ] +
xx [ 32 ] + xx [ 174 ] + xx [ 198 ] ; xx [ 32 ] = xx [ 209 ] - xx [ 68 ] * xx
[ 153 ] - ( xx [ 142 ] + xx [ 210 ] ) + xx [ 200 ] ; xx [ 71 ] = xx [ 68 ] *
xx [ 149 ] + xx [ 203 ] + xx [ 205 ] - xx [ 207 ] ; xx [ 158 ] = xx [ 64 ] -
( xx [ 45 ] * xx [ 168 ] + xx [ 67 ] * xx [ 46 ] + xx [ 179 ] - ( xx [ 180 ]
- xx [ 58 ] ) * xx [ 3 ] - ( xx [ 75 ] * xx [ 77 ] + xx [ 136 ] * xx [ 83 ] )
+ xx [ 182 ] - ( xx [ 13 ] * xx [ 16 ] + xx [ 19 ] * xx [ 139 ] ) * xx [ 3 ]
+ xx [ 3 ] * ( xx [ 185 ] - ( xx [ 55 ] * xx [ 16 ] + xx [ 19 ] * xx [ 188 ]
) * xx [ 3 ] ) + xx [ 191 ] + xx [ 3 ] * ( xx [ 189 ] * xx [ 20 ] - xx [ 11 ]
* xx [ 190 ] ) + ( xx [ 193 ] + xx [ 3 ] * ( xx [ 9 ] * xx [ 20 ] - xx [ 11 ]
* xx [ 194 ] ) ) * xx [ 3 ] - xx [ 23 ] * ( xx [ 168 ] * xx [ 37 ] - xx [ 67
] * xx [ 33 ] - ( xx [ 195 ] * xx [ 157 ] + xx [ 196 ] * xx [ 173 ] ) - xx [
3 ] * xx [ 5 ] * ( xx [ 3 ] * xx [ 22 ] * xx [ 22 ] - xx [ 1 ] ) * state [ 7
] ) - ( xx [ 23 ] * ( xx [ 208 ] * xx [ 125 ] + xx [ 119 ] * xx [ 201 ] - (
xx [ 211 ] * xx [ 206 ] + xx [ 141 ] * xx [ 212 ] ) - xx [ 3 ] * xx [ 70 ] *
( xx [ 3 ] * xx [ 106 ] * xx [ 106 ] - xx [ 1 ] ) * state [ 15 ] ) + xx [ 91
] * xx [ 208 ] + xx [ 111 ] * xx [ 201 ] + xx [ 132 ] * xx [ 202 ] + xx [ 169
] * xx [ 135 ] + xx [ 109 ] * ( xx [ 3 ] * xx [ 213 ] * xx [ 97 ] - xx [ 80 ]
) - ( xx [ 164 ] * xx [ 199 ] - xx [ 165 ] * xx [ 197 ] ) + xx [ 216 ] - ( xx
[ 217 ] * xx [ 97 ] - xx [ 95 ] * xx [ 214 ] ) * xx [ 3 ] + xx [ 3 ] * ( ( xx
[ 95 ] * xx [ 140 ] + xx [ 116 ] * xx [ 97 ] ) * xx [ 3 ] - xx [ 219 ] ) + xx
[ 222 ] - ( xx [ 98 ] * xx [ 223 ] - xx [ 220 ] * xx [ 96 ] ) * xx [ 3 ] + xx
[ 3 ] * ( ( xx [ 118 ] * xx [ 96 ] + xx [ 98 ] * xx [ 114 ] ) * xx [ 3 ] - xx
[ 225 ] ) + xx [ 228 ] - xx [ 3 ] * ( xx [ 226 ] * xx [ 103 ] + xx [ 92 ] *
xx [ 227 ] ) + xx [ 3 ] * ( xx [ 3 ] * ( xx [ 92 ] * xx [ 101 ] - xx [ 85 ] *
xx [ 103 ] ) - xx [ 230 ] ) ) ) + xx [ 232 ] * xx [ 57 ] + xx [ 233 ] * xx [
60 ] + xx [ 36 ] * xx [ 27 ] - xx [ 176 ] * xx [ 171 ] - xx [ 54 ] * xx [ 148
] + xx [ 21 ] * xx [ 122 ] - xx [ 150 ] * xx [ 110 ] ; xx [ 159 ] = xx [ 232
] * xx [ 167 ] - ( xx [ 45 ] * xx [ 31 ] + xx [ 46 ] * xx [ 73 ] + xx [ 75 ]
* xx [ 100 ] - xx [ 65 ] * xx [ 83 ] + xx [ 181 ] - xx [ 3 ] * ( xx [ 58 ] +
xx [ 180 ] ) + xx [ 3 ] * ( xx [ 19 ] * xx [ 13 ] - xx [ 139 ] * xx [ 16 ] )
- xx [ 184 ] + xx [ 3 ] * ( xx [ 3 ] * ( xx [ 19 ] * xx [ 55 ] - xx [ 188 ] *
xx [ 16 ] ) - xx [ 52 ] ) + ( xx [ 190 ] * xx [ 20 ] + xx [ 11 ] * xx [ 189 ]
) * xx [ 3 ] - xx [ 192 ] + xx [ 3 ] * ( ( xx [ 194 ] * xx [ 20 ] + xx [ 11 ]
* xx [ 9 ] ) * xx [ 3 ] - xx [ 17 ] ) - xx [ 23 ] * ( xx [ 31 ] * xx [ 37 ] -
xx [ 33 ] * xx [ 73 ] + xx [ 173 ] * xx [ 157 ] - xx [ 196 ] * xx [ 195 ] -
xx [ 3 ] * xx [ 5 ] * xx [ 3 ] * xx [ 22 ] * xx [ 24 ] * state [ 7 ] ) - ( xx
[ 23 ] * ( xx [ 32 ] * xx [ 125 ] + xx [ 119 ] * xx [ 71 ] + xx [ 141 ] * xx
[ 211 ] - xx [ 212 ] * xx [ 206 ] - xx [ 3 ] * xx [ 70 ] * xx [ 3 ] * xx [
106 ] * xx [ 107 ] * state [ 15 ] ) + xx [ 91 ] * xx [ 32 ] + xx [ 111 ] * xx
[ 71 ] + xx [ 3 ] * ( xx [ 95 ] * xx [ 217 ] + xx [ 214 ] * xx [ 97 ] ) - xx
[ 218 ] - ( xx [ 165 ] * xx [ 199 ] + xx [ 164 ] * xx [ 197 ] + xx [ 3 ] * xx
[ 95 ] * xx [ 213 ] * xx [ 109 ] ) + xx [ 3 ] * ( xx [ 3 ] * ( xx [ 140 ] *
xx [ 97 ] - xx [ 95 ] * xx [ 116 ] ) - xx [ 215 ] ) - ( xx [ 68 ] * xx [ 132
] - xx [ 135 ] * xx [ 202 ] ) + xx [ 3 ] * ( xx [ 223 ] * xx [ 96 ] + xx [ 98
] * xx [ 220 ] ) - xx [ 224 ] + xx [ 3 ] * ( xx [ 3 ] * ( xx [ 98 ] * xx [
118 ] - xx [ 114 ] * xx [ 96 ] ) - xx [ 221 ] ) - ( xx [ 229 ] + ( xx [ 227 ]
* xx [ 103 ] - xx [ 92 ] * xx [ 226 ] ) * xx [ 3 ] ) + xx [ 3 ] * ( ( xx [
101 ] * xx [ 103 ] + xx [ 92 ] * xx [ 85 ] ) * xx [ 3 ] - xx [ 231 ] ) ) ) +
xx [ 233 ] * xx [ 53 ] + xx [ 36 ] * xx [ 25 ] - xx [ 115 ] * xx [ 176 ] - xx
[ 54 ] * xx [ 89 ] - xx [ 21 ] * xx [ 105 ] - xx [ 150 ] * xx [ 51 ] ; memcpy
( xx + 19 , xx + 143 , 4 * sizeof ( double ) ) ; factorAndSolveSymmetric ( xx
+ 19 , 2 , xx + 23 , ii + 0 , xx + 158 , xx + 16 , xx + 176 ) ; xx [ 1 ] = (
xx [ 16 ] * xx [ 27 ] + xx [ 25 ] * xx [ 17 ] - xx [ 90 ] ) / xx [ 10 ] ; xx
[ 5 ] = xx [ 49 ] + xx [ 2 ] * xx [ 1 ] ; xx [ 9 ] = xx [ 5 ] * xx [ 8 ] ; xx
[ 11 ] = xx [ 5 ] - ( xx [ 30 ] + xx [ 9 ] * xx [ 8 ] ) * xx [ 3 ] ; xx [ 5 ]
= xx [ 11 ] - input [ 32 ] + xx [ 102 ] ; xx [ 13 ] = xx [ 108 ] + xx [ 3 ] *
( xx [ 6 ] * xx [ 9 ] - xx [ 112 ] ) ; xx [ 9 ] = xx [ 87 ] + xx [ 26 ] * xx
[ 1 ] + xx [ 13 ] * xx [ 46 ] - xx [ 45 ] * xx [ 11 ] - input [ 34 ] - xx [
44 ] ; xx [ 11 ] = ( xx [ 16 ] * xx [ 60 ] + xx [ 53 ] * xx [ 17 ] - ( xx [
48 ] + xx [ 9 ] + xx [ 5 ] * xx [ 7 ] ) ) / xx [ 47 ] ; xx [ 19 ] = xx [ 5 ]
+ xx [ 50 ] * xx [ 11 ] ; xx [ 5 ] = xx [ 13 ] - input [ 30 ] + xx [ 40 ] +
xx [ 38 ] * xx [ 11 ] ; xx [ 13 ] = xx [ 5 ] * xx [ 12 ] ; xx [ 20 ] = xx [
19 ] * xx [ 12 ] ; xx [ 21 ] = xx [ 19 ] - ( xx [ 15 ] * xx [ 13 ] + xx [ 20
] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 19 ] = xx [ 21 ] - input [ 38 ] + xx [ 131
] ; xx [ 22 ] = xx [ 5 ] + xx [ 3 ] * ( xx [ 15 ] * xx [ 20 ] - xx [ 13 ] *
xx [ 12 ] ) ; xx [ 5 ] = ( xx [ 57 ] * xx [ 16 ] + xx [ 167 ] * xx [ 17 ] - (
xx [ 134 ] + xx [ 9 ] + xx [ 42 ] * xx [ 11 ] + xx [ 22 ] * xx [ 83 ] - xx [
75 ] * xx [ 21 ] - input [ 40 ] - xx [ 28 ] + xx [ 19 ] * xx [ 7 ] ) ) / xx [
74 ] ; xx [ 9 ] = xx [ 19 ] + xx [ 82 ] * xx [ 5 ] ; xx [ 13 ] = ( xx [ 16 ]
- ( xx [ 9 ] - ( xx [ 18 ] * xx [ 9 ] * xx [ 18 ] - ( xx [ 22 ] - input [ 36
] + xx [ 34 ] + xx [ 62 ] * xx [ 5 ] ) * xx [ 18 ] * xx [ 14 ] ) * xx [ 3 ] )
) / xx [ 72 ] ; xx [ 9 ] = xx [ 18 ] * xx [ 13 ] ; xx [ 19 ] = xx [ 13 ] - xx
[ 3 ] * xx [ 18 ] * xx [ 9 ] ; xx [ 18 ] = xx [ 3 ] * xx [ 9 ] * xx [ 14 ] ;
xx [ 9 ] = xx [ 5 ] - ( xx [ 66 ] * xx [ 19 ] + xx [ 78 ] * xx [ 18 ] ) ; xx
[ 5 ] = xx [ 129 ] + xx [ 18 ] + xx [ 83 ] * xx [ 9 ] ; xx [ 14 ] = xx [ 19 ]
+ xx [ 7 ] * xx [ 9 ] - xx [ 75 ] * xx [ 9 ] ; xx [ 18 ] = xx [ 12 ] * xx [
14 ] ; xx [ 19 ] = xx [ 5 ] * xx [ 12 ] ; xx [ 20 ] = xx [ 5 ] - ( xx [ 15 ]
* xx [ 18 ] + xx [ 19 ] * xx [ 12 ] ) * xx [ 3 ] ; xx [ 5 ] = xx [ 14 ] + xx
[ 3 ] * ( xx [ 15 ] * xx [ 19 ] - xx [ 18 ] * xx [ 12 ] ) ; xx [ 12 ] = xx [
11 ] - ( xx [ 56 ] * xx [ 9 ] + xx [ 69 ] * xx [ 20 ] + xx [ 5 ] * xx [ 63 ]
) ; xx [ 11 ] = xx [ 9 ] + xx [ 12 ] ; xx [ 14 ] = xx [ 5 ] + xx [ 7 ] * xx [
12 ] + xx [ 117 ] - xx [ 45 ] * xx [ 11 ] ; xx [ 5 ] = xx [ 1 ] - ( xx [ 11 ]
* xx [ 79 ] + ( xx [ 14 ] + xx [ 3 ] * ( xx [ 6 ] * ( xx [ 20 ] + xx [ 121 ]
+ xx [ 11 ] * xx [ 46 ] ) * xx [ 8 ] - xx [ 8 ] * xx [ 14 ] * xx [ 8 ] ) ) *
xx [ 43 ] ) ; xx [ 1 ] = ( xx [ 84 ] + xx [ 16 ] * xx [ 110 ] + xx [ 51 ] *
xx [ 17 ] ) / xx [ 10 ] ; xx [ 6 ] = xx [ 237 ] - xx [ 2 ] * xx [ 1 ] ; xx [
2 ] = xx [ 86 ] * xx [ 6 ] ; xx [ 8 ] = xx [ 6 ] + xx [ 3 ] * ( xx [ 76 ] -
xx [ 2 ] * xx [ 86 ] ) ; xx [ 6 ] = xx [ 8 ] - input [ 14 ] + xx [ 239 ] ; xx
[ 10 ] = xx [ 240 ] - ( xx [ 104 ] * xx [ 2 ] + xx [ 113 ] ) * xx [ 3 ] ; xx
[ 2 ] = xx [ 241 ] - ( input [ 16 ] + xx [ 235 ] - xx [ 26 ] * xx [ 1 ] + xx
[ 91 ] * xx [ 8 ] - xx [ 111 ] * xx [ 10 ] ) ; xx [ 8 ] = ( xx [ 128 ] + xx [
6 ] * xx [ 7 ] - xx [ 2 ] + xx [ 105 ] * xx [ 17 ] - xx [ 16 ] * xx [ 122 ] )
/ xx [ 61 ] ; xx [ 14 ] = xx [ 6 ] - xx [ 8 ] * xx [ 133 ] ; xx [ 6 ] = xx [
10 ] - input [ 12 ] + xx [ 35 ] - xx [ 126 ] * xx [ 8 ] ; xx [ 10 ] = xx [ 93
] * xx [ 6 ] ; xx [ 15 ] = xx [ 93 ] * xx [ 14 ] ; xx [ 18 ] = xx [ 14 ] + xx
[ 3 ] * ( xx [ 99 ] * xx [ 10 ] - xx [ 15 ] * xx [ 93 ] ) ; xx [ 14 ] = xx [
18 ] - input [ 8 ] + xx [ 234 ] ; xx [ 19 ] = xx [ 6 ] - ( xx [ 99 ] * xx [
15 ] + xx [ 10 ] * xx [ 93 ] ) * xx [ 3 ] ; xx [ 6 ] = xx [ 2 ] - xx [ 8 ] *
xx [ 41 ] - ( xx [ 132 ] * xx [ 18 ] - xx [ 135 ] * xx [ 19 ] ) - input [ 10
] + xx [ 29 ] ; xx [ 2 ] = ( xx [ 161 ] + xx [ 14 ] * xx [ 7 ] - xx [ 6 ] +
xx [ 16 ] * xx [ 148 ] + xx [ 89 ] * xx [ 17 ] ) / xx [ 152 ] ; xx [ 10 ] =
xx [ 14 ] - xx [ 2 ] * xx [ 160 ] ; xx [ 14 ] = xx [ 19 ] - input [ 6 ] + xx
[ 137 ] - xx [ 151 ] * xx [ 2 ] ; xx [ 15 ] = xx [ 94 ] * xx [ 14 ] ; xx [ 18
] = xx [ 94 ] * xx [ 10 ] ; xx [ 19 ] = xx [ 10 ] + xx [ 3 ] * ( xx [ 4 ] *
xx [ 15 ] - xx [ 18 ] * xx [ 94 ] ) ; xx [ 10 ] = ( xx [ 39 ] + xx [ 7 ] * (
xx [ 19 ] - input [ 2 ] - xx [ 163 ] ) - ( xx [ 6 ] - xx [ 2 ] * xx [ 130 ] -
( xx [ 164 ] * xx [ 19 ] - xx [ 165 ] * ( xx [ 14 ] - ( xx [ 4 ] * xx [ 18 ]
+ xx [ 15 ] * xx [ 94 ] ) * xx [ 3 ] ) ) - input [ 4 ] + xx [ 170 ] ) + xx [
16 ] * xx [ 171 ] + xx [ 115 ] * xx [ 17 ] ) / xx [ 172 ] ; xx [ 6 ] = xx [
124 ] - xx [ 10 ] * xx [ 165 ] ; xx [ 14 ] = xx [ 6 ] * xx [ 94 ] ; xx [ 15 ]
= xx [ 7 ] * xx [ 10 ] + xx [ 164 ] * xx [ 10 ] ; xx [ 16 ] = xx [ 15 ] * xx
[ 94 ] ; xx [ 17 ] = xx [ 3 ] * ( xx [ 14 ] * xx [ 94 ] - xx [ 4 ] * xx [ 16
] ) - xx [ 6 ] ; xx [ 6 ] = ( xx [ 4 ] * xx [ 14 ] + xx [ 16 ] * xx [ 94 ] )
* xx [ 3 ] - xx [ 15 ] ; xx [ 4 ] = xx [ 2 ] + xx [ 10 ] * xx [ 186 ] + xx [
156 ] * xx [ 17 ] + xx [ 175 ] * xx [ 6 ] ; xx [ 2 ] = xx [ 10 ] + xx [ 4 ] ;
xx [ 14 ] = xx [ 17 ] + xx [ 247 ] + xx [ 2 ] * xx [ 135 ] ; xx [ 15 ] = xx [
6 ] - xx [ 4 ] * xx [ 7 ] + xx [ 59 ] - xx [ 2 ] * xx [ 132 ] ; xx [ 6 ] = xx
[ 93 ] * xx [ 15 ] ; xx [ 16 ] = xx [ 93 ] * xx [ 14 ] ; xx [ 17 ] = xx [ 14
] + xx [ 3 ] * ( xx [ 99 ] * xx [ 6 ] - xx [ 16 ] * xx [ 93 ] ) ; xx [ 14 ] =
xx [ 15 ] - ( xx [ 99 ] * xx [ 16 ] + xx [ 6 ] * xx [ 93 ] ) * xx [ 3 ] ; xx
[ 6 ] = xx [ 8 ] + xx [ 2 ] * xx [ 166 ] + xx [ 17 ] * xx [ 123 ] + xx [ 154
] * xx [ 14 ] ; xx [ 8 ] = xx [ 2 ] + xx [ 6 ] ; xx [ 2 ] = xx [ 14 ] - xx [
6 ] * xx [ 7 ] + xx [ 236 ] - xx [ 8 ] * xx [ 91 ] ; xx [ 7 ] = xx [ 1 ] + xx
[ 43 ] * ( xx [ 2 ] - ( xx [ 104 ] * xx [ 86 ] * ( xx [ 17 ] + xx [ 243 ] +
xx [ 8 ] * xx [ 111 ] ) + xx [ 86 ] * xx [ 2 ] * xx [ 86 ] ) * xx [ 3 ] ) -
xx [ 8 ] * xx [ 79 ] ; logVector [ 0 ] = state [ 0 ] ; logVector [ 1 ] =
state [ 1 ] ; logVector [ 2 ] = xx [ 0 ] * state [ 2 ] ; logVector [ 3 ] = xx
[ 0 ] * state [ 3 ] ; logVector [ 4 ] = xx [ 0 ] * state [ 4 ] ; logVector [
5 ] = xx [ 0 ] * state [ 5 ] ; logVector [ 6 ] = xx [ 0 ] * state [ 6 ] ;
logVector [ 7 ] = xx [ 0 ] * state [ 7 ] ; logVector [ 8 ] = xx [ 0 ] * state
[ 8 ] ; logVector [ 9 ] = xx [ 0 ] * state [ 9 ] ; logVector [ 10 ] = xx [ 0
] * state [ 10 ] ; logVector [ 11 ] = xx [ 0 ] * state [ 11 ] ; logVector [
12 ] = xx [ 0 ] * state [ 12 ] ; logVector [ 13 ] = xx [ 0 ] * state [ 13 ] ;
logVector [ 14 ] = xx [ 0 ] * state [ 14 ] ; logVector [ 15 ] = xx [ 0 ] *
state [ 15 ] ; logVector [ 16 ] = xx [ 13 ] ; logVector [ 17 ] = xx [ 0 ] *
xx [ 9 ] ; logVector [ 18 ] = xx [ 0 ] * xx [ 12 ] ; logVector [ 19 ] = xx [
0 ] * xx [ 5 ] ; logVector [ 20 ] = - ( xx [ 0 ] * xx [ 10 ] ) ; logVector [
21 ] = - ( xx [ 4 ] * xx [ 0 ] ) ; logVector [ 22 ] = - ( xx [ 6 ] * xx [ 0 ]
) ; logVector [ 23 ] = - ( xx [ 7 ] * xx [ 0 ] ) ; logVector [ 24 ] = xx [ 88
] * xx [ 0 ] ; logVector [ 25 ] = - ( xx [ 81 ] * xx [ 0 ] ) ; logVector [ 26
] = - ( xx [ 0 ] * ( xx [ 11 ] + xx [ 5 ] - ( xx [ 8 ] + xx [ 7 ] ) ) ) ;
errorResult [ 0 ] = 0.0 ; return NULL ; }
