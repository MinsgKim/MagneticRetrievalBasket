#include <math.h>
#include <string.h>
#include "pm_std.h"
#include "sm_std.h"
#include "ne_std.h"
#include "ne_dae.h"
#include "sm_ssci_run_time_errors.h"
#include "sm_RuntimeDerivedValuesBundle.h"
#include "MagneticBasket_Simscape_Optimizer_dda62cd9_1_geometries.h"
PmfMessageId MagneticBasket_Simscape_Optimizer_dda62cd9_1_compOutputs ( const
RuntimeDerivedValuesBundle * rtdv , const double * state , const int *
modeVector , const double * input , const double * inputDot , const double *
inputDdot , const double * discreteState , double * output ,
NeuDiagnosticManager * neDiagMgr ) { const double * rtdvd = rtdv -> mDoubles
. mValues ; const int * rtdvi = rtdv -> mInts . mValues ; double xx [ 81 ] ;
( void ) rtdvd ; ( void ) rtdvi ; ( void ) modeVector ; ( void ) input ; (
void ) inputDot ; ( void ) inputDdot ; ( void ) discreteState ; ( void )
neDiagMgr ; xx [ 0 ] = 2.0 ; xx [ 1 ] = 1.0e-3 ; xx [ 2 ] = 0.5 ; xx [ 3 ] =
xx [ 2 ] * state [ 8 ] ; xx [ 4 ] = sin ( xx [ 3 ] ) ; xx [ 5 ] = xx [ 1 ] *
xx [ 4 ] ; xx [ 6 ] = xx [ 0 ] * xx [ 5 ] * xx [ 4 ] - xx [ 1 ] ; xx [ 7 ] =
0.045 ; xx [ 8 ] = 0.0 ; xx [ 9 ] = cos ( xx [ 3 ] ) ; xx [ 3 ] =
1.500000000000001e-3 + xx [ 0 ] * xx [ 9 ] * xx [ 5 ] ; xx [ 5 ] =
0.9657342571406726 ; xx [ 10 ] = 0.2595329354531968 ; xx [ 11 ] = xx [ 5 ] *
xx [ 9 ] + xx [ 10 ] * xx [ 4 ] ; xx [ 12 ] = xx [ 11 ] * xx [ 11 ] ; xx [ 13
] = 1.0 ; xx [ 14 ] = xx [ 0 ] * xx [ 12 ] - xx [ 13 ] ; xx [ 15 ] = xx [ 5 ]
* xx [ 4 ] - xx [ 10 ] * xx [ 9 ] ; xx [ 5 ] = xx [ 0 ] * xx [ 11 ] * xx [ 15
] ; xx [ 10 ] = xx [ 9 ] * xx [ 9 ] ; xx [ 11 ] = xx [ 0 ] * xx [ 10 ] - xx [
13 ] ; xx [ 16 ] = xx [ 0 ] * xx [ 9 ] * xx [ 4 ] ; xx [ 17 ] = 2.0e-3 ; xx [
18 ] = xx [ 2 ] * state [ 10 ] ; xx [ 19 ] = sin ( xx [ 18 ] ) ; xx [ 20 ] =
xx [ 1 ] * xx [ 19 ] ; xx [ 21 ] = xx [ 17 ] - xx [ 0 ] * xx [ 20 ] * xx [ 19
] ; xx [ 22 ] = cos ( xx [ 18 ] ) ; xx [ 18 ] = xx [ 0 ] * xx [ 22 ] * xx [
20 ] ; xx [ 20 ] = xx [ 18 ] * xx [ 4 ] ; xx [ 23 ] = xx [ 21 ] * xx [ 4 ] ;
xx [ 24 ] = xx [ 21 ] - xx [ 0 ] * ( xx [ 9 ] * xx [ 20 ] + xx [ 23 ] * xx [
4 ] ) - xx [ 6 ] ; xx [ 21 ] = ( xx [ 9 ] * xx [ 23 ] - xx [ 20 ] * xx [ 4 ]
) * xx [ 0 ] + xx [ 18 ] + xx [ 3 ] ; xx [ 18 ] = 0.9258049549196806 ; xx [
20 ] = xx [ 9 ] * xx [ 22 ] - xx [ 4 ] * xx [ 19 ] ; xx [ 23 ] = xx [ 9 ] *
xx [ 19 ] + xx [ 22 ] * xx [ 4 ] ; xx [ 9 ] = 0.3780015680472346 ; xx [ 19 ]
= xx [ 18 ] * xx [ 20 ] + xx [ 23 ] * xx [ 9 ] ; xx [ 22 ] = xx [ 19 ] * xx [
19 ] ; xx [ 25 ] = xx [ 0 ] * xx [ 22 ] - xx [ 13 ] ; xx [ 26 ] = xx [ 9 ] *
xx [ 20 ] - xx [ 23 ] * xx [ 18 ] ; xx [ 9 ] = xx [ 0 ] * xx [ 19 ] * xx [ 26
] ; xx [ 18 ] = xx [ 20 ] * xx [ 20 ] ; xx [ 19 ] = xx [ 0 ] * xx [ 18 ] - xx
[ 13 ] ; xx [ 27 ] = xx [ 0 ] * xx [ 23 ] * xx [ 20 ] ; xx [ 28 ] = xx [ 2 ]
* state [ 12 ] ; xx [ 29 ] = sin ( xx [ 28 ] ) ; xx [ 30 ] = xx [ 1 ] * xx [
29 ] ; xx [ 31 ] = xx [ 17 ] - xx [ 0 ] * xx [ 30 ] * xx [ 29 ] ; xx [ 32 ] =
cos ( xx [ 28 ] ) ; xx [ 28 ] = xx [ 0 ] * xx [ 32 ] * xx [ 30 ] ; xx [ 30 ]
= xx [ 23 ] * xx [ 28 ] ; xx [ 33 ] = xx [ 23 ] * xx [ 31 ] ; xx [ 34 ] = xx
[ 31 ] - xx [ 0 ] * ( xx [ 30 ] * xx [ 20 ] + xx [ 23 ] * xx [ 33 ] ) + xx [
24 ] ; xx [ 31 ] = ( xx [ 33 ] * xx [ 20 ] - xx [ 23 ] * xx [ 30 ] ) * xx [ 0
] + xx [ 28 ] + xx [ 21 ] ; xx [ 28 ] = 0.9301380840958209 ; xx [ 30 ] = xx [
23 ] * xx [ 29 ] - xx [ 32 ] * xx [ 20 ] ; xx [ 33 ] = xx [ 29 ] * xx [ 20 ]
+ xx [ 23 ] * xx [ 32 ] ; xx [ 20 ] = 0.3672099460997151 ; xx [ 29 ] = xx [
28 ] * xx [ 30 ] - xx [ 33 ] * xx [ 20 ] ; xx [ 32 ] = xx [ 29 ] * xx [ 29 ]
; xx [ 35 ] = xx [ 0 ] * xx [ 32 ] - xx [ 13 ] ; xx [ 36 ] = xx [ 20 ] * xx [
30 ] + xx [ 33 ] * xx [ 28 ] ; xx [ 20 ] = xx [ 0 ] * xx [ 36 ] * xx [ 29 ] ;
xx [ 28 ] = xx [ 30 ] * xx [ 30 ] ; xx [ 29 ] = xx [ 0 ] * xx [ 28 ] - xx [
13 ] ; xx [ 37 ] = xx [ 0 ] * xx [ 33 ] * xx [ 30 ] ; xx [ 38 ] = xx [ 2 ] *
state [ 14 ] ; xx [ 39 ] = sin ( xx [ 38 ] ) ; xx [ 40 ] = xx [ 1 ] * xx [ 39
] ; xx [ 41 ] = xx [ 17 ] - xx [ 0 ] * xx [ 40 ] * xx [ 39 ] ; xx [ 42 ] = xx
[ 33 ] * xx [ 41 ] ; xx [ 43 ] = cos ( xx [ 38 ] ) ; xx [ 38 ] = xx [ 0 ] *
xx [ 43 ] * xx [ 40 ] ; xx [ 40 ] = xx [ 33 ] * xx [ 38 ] ; xx [ 44 ] = xx [
43 ] * xx [ 30 ] + xx [ 33 ] * xx [ 39 ] ; xx [ 45 ] = 0.9766452861856919 ;
xx [ 46 ] = 0.2148580577294411 ; xx [ 47 ] = xx [ 39 ] * xx [ 30 ] - xx [ 33
] * xx [ 43 ] ; xx [ 39 ] = xx [ 44 ] * xx [ 45 ] + xx [ 46 ] * xx [ 47 ] ;
xx [ 43 ] = xx [ 39 ] * xx [ 39 ] ; xx [ 48 ] = xx [ 0 ] * xx [ 43 ] - xx [
13 ] ; xx [ 49 ] = xx [ 45 ] * xx [ 47 ] - xx [ 44 ] * xx [ 46 ] ; xx [ 45 ]
= xx [ 0 ] * xx [ 39 ] * xx [ 49 ] ; xx [ 39 ] = xx [ 44 ] * xx [ 44 ] ; xx [
46 ] = xx [ 0 ] * xx [ 39 ] - xx [ 13 ] ; xx [ 50 ] = xx [ 0 ] * xx [ 44 ] *
xx [ 47 ] ; xx [ 44 ] = xx [ 2 ] * state [ 6 ] ; xx [ 51 ] = sin ( xx [ 44 ]
) ; xx [ 52 ] = xx [ 1 ] * xx [ 51 ] ; xx [ 53 ] = xx [ 0 ] * xx [ 52 ] * xx
[ 51 ] - xx [ 17 ] ; xx [ 54 ] = xx [ 2 ] * state [ 4 ] ; xx [ 55 ] = sin (
xx [ 54 ] ) ; xx [ 56 ] = 0.7071067811865476 ; xx [ 57 ] = xx [ 2 ] * state [
2 ] ; xx [ 2 ] = xx [ 56 ] * cos ( xx [ 57 ] ) ; xx [ 58 ] = xx [ 56 ] * sin
( xx [ 57 ] ) ; xx [ 57 ] = xx [ 2 ] + xx [ 58 ] ; xx [ 59 ] = xx [ 56 ] * xx
[ 57 ] ; xx [ 60 ] = xx [ 2 ] - xx [ 58 ] ; xx [ 2 ] = xx [ 60 ] * xx [ 56 ]
; xx [ 58 ] = xx [ 59 ] + xx [ 2 ] ; xx [ 61 ] = xx [ 59 ] - xx [ 2 ] ; xx [
2 ] = cos ( xx [ 54 ] ) ; xx [ 54 ] = xx [ 55 ] * xx [ 58 ] + xx [ 61 ] * xx
[ 2 ] ; xx [ 59 ] = cos ( xx [ 44 ] ) ; xx [ 44 ] = xx [ 0 ] * xx [ 59 ] * xx
[ 52 ] ; xx [ 52 ] = xx [ 54 ] * xx [ 44 ] ; xx [ 62 ] = xx [ 61 ] * xx [ 55
] - xx [ 2 ] * xx [ 58 ] ; xx [ 63 ] = xx [ 53 ] * xx [ 54 ] ; xx [ 64 ] = xx
[ 1 ] * xx [ 55 ] ; xx [ 65 ] = xx [ 0 ] * xx [ 64 ] * xx [ 55 ] - xx [ 17 ]
; xx [ 17 ] = xx [ 0 ] * xx [ 2 ] * xx [ 64 ] ; xx [ 2 ] = xx [ 61 ] * xx [
17 ] ; xx [ 55 ] = xx [ 65 ] * xx [ 61 ] ; xx [ 64 ] = xx [ 60 ] * xx [ 1 ] ;
xx [ 1 ] = xx [ 0 ] * xx [ 60 ] * xx [ 64 ] - 2.500000000000001e-3 ; xx [ 60
] = xx [ 0 ] * xx [ 64 ] * xx [ 57 ] ; xx [ 57 ] = xx [ 56 ] * xx [ 56 ] * xx
[ 60 ] ; xx [ 64 ] = xx [ 56 ] * xx [ 1 ] * xx [ 56 ] ; xx [ 56 ] = xx [ 1 ]
- xx [ 0 ] * ( xx [ 57 ] + xx [ 64 ] ) + state [ 0 ] + 0.014 ; xx [ 1 ] = xx
[ 65 ] + xx [ 0 ] * ( xx [ 2 ] * xx [ 58 ] - xx [ 61 ] * xx [ 55 ] ) + xx [
56 ] ; xx [ 65 ] = xx [ 60 ] + ( xx [ 64 ] - xx [ 57 ] ) * xx [ 0 ] ; xx [ 57
] = xx [ 17 ] - ( xx [ 55 ] * xx [ 58 ] + xx [ 61 ] * xx [ 2 ] ) * xx [ 0 ] -
xx [ 65 ] ; xx [ 2 ] = xx [ 59 ] * xx [ 62 ] + xx [ 54 ] * xx [ 51 ] ; xx [
17 ] = 0.1585578743085479 ; xx [ 55 ] = 0.9873496850127389 ; xx [ 60 ] = xx [
54 ] * xx [ 59 ] - xx [ 51 ] * xx [ 62 ] ; xx [ 51 ] = xx [ 2 ] * xx [ 17 ] +
xx [ 55 ] * xx [ 60 ] ; xx [ 59 ] = xx [ 51 ] * xx [ 51 ] ; xx [ 64 ] = xx [
0 ] * xx [ 59 ] - xx [ 13 ] ; xx [ 66 ] = xx [ 17 ] * xx [ 60 ] - xx [ 2 ] *
xx [ 55 ] ; xx [ 17 ] = xx [ 0 ] * xx [ 51 ] * xx [ 66 ] ; xx [ 51 ] = xx [ 2
] * xx [ 2 ] ; xx [ 55 ] = xx [ 0 ] * xx [ 51 ] - xx [ 13 ] ; xx [ 67 ] = xx
[ 0 ] * xx [ 2 ] * xx [ 60 ] ; xx [ 2 ] = 0.3290701976886963 ; xx [ 68 ] =
0.9443054616982379 ; xx [ 69 ] = xx [ 2 ] * xx [ 62 ] + xx [ 54 ] * xx [ 68 ]
; xx [ 70 ] = xx [ 69 ] * xx [ 69 ] ; xx [ 71 ] = xx [ 0 ] * xx [ 70 ] - xx [
13 ] ; xx [ 72 ] = xx [ 68 ] * xx [ 62 ] - xx [ 54 ] * xx [ 2 ] ; xx [ 2 ] =
xx [ 0 ] * xx [ 69 ] * xx [ 72 ] ; xx [ 68 ] = xx [ 62 ] * xx [ 62 ] ; xx [
69 ] = xx [ 0 ] * xx [ 68 ] - xx [ 13 ] ; xx [ 73 ] = xx [ 0 ] * xx [ 54 ] *
xx [ 62 ] ; xx [ 74 ] = 0.2849406335136161 ; xx [ 75 ] = 0.9585451660578437 ;
xx [ 76 ] = xx [ 74 ] * xx [ 58 ] - xx [ 61 ] * xx [ 75 ] ; xx [ 77 ] = xx [
76 ] * xx [ 76 ] ; xx [ 78 ] = xx [ 0 ] * xx [ 77 ] - xx [ 13 ] ; xx [ 79 ] =
xx [ 75 ] * xx [ 58 ] + xx [ 61 ] * xx [ 74 ] ; xx [ 74 ] = xx [ 0 ] * xx [
79 ] * xx [ 76 ] ; xx [ 75 ] = xx [ 58 ] * xx [ 58 ] ; xx [ 76 ] = xx [ 0 ] *
xx [ 75 ] - xx [ 13 ] ; xx [ 80 ] = xx [ 0 ] * xx [ 61 ] * xx [ 58 ] ; output
[ 0 ] = state [ 8 ] ; output [ 1 ] = state [ 10 ] ; output [ 2 ] = state [ 12
] ; output [ 3 ] = state [ 14 ] ; output [ 4 ] = state [ 16 ] ; output [ 5 ]
= state [ 6 ] ; output [ 6 ] = state [ 4 ] ; output [ 7 ] = - ( xx [ 6 ] + xx
[ 7 ] ) ; output [ 8 ] = xx [ 8 ] ; output [ 9 ] = xx [ 3 ] ; output [ 10 ] =
xx [ 14 ] ; output [ 11 ] = xx [ 8 ] ; output [ 12 ] = xx [ 5 ] ; output [ 13
] = xx [ 8 ] ; output [ 14 ] = ( xx [ 12 ] + xx [ 15 ] * xx [ 15 ] ) * xx [ 0
] - xx [ 13 ] ; output [ 15 ] = xx [ 8 ] ; output [ 16 ] = - xx [ 5 ] ;
output [ 17 ] = xx [ 8 ] ; output [ 18 ] = xx [ 14 ] ; output [ 19 ] = xx [
11 ] ; output [ 20 ] = xx [ 8 ] ; output [ 21 ] = xx [ 16 ] ; output [ 22 ] =
xx [ 8 ] ; output [ 23 ] = ( xx [ 10 ] + xx [ 4 ] * xx [ 4 ] ) * xx [ 0 ] -
xx [ 13 ] ; output [ 24 ] = xx [ 8 ] ; output [ 25 ] = - xx [ 16 ] ; output [
26 ] = xx [ 8 ] ; output [ 27 ] = xx [ 11 ] ; output [ 28 ] = xx [ 24 ] - xx
[ 7 ] ; output [ 29 ] = xx [ 8 ] ; output [ 30 ] = xx [ 21 ] ; output [ 31 ]
= xx [ 25 ] ; output [ 32 ] = xx [ 8 ] ; output [ 33 ] = - xx [ 9 ] ; output
[ 34 ] = xx [ 8 ] ; output [ 35 ] = ( xx [ 22 ] + xx [ 26 ] * xx [ 26 ] ) *
xx [ 0 ] - xx [ 13 ] ; output [ 36 ] = xx [ 8 ] ; output [ 37 ] = xx [ 9 ] ;
output [ 38 ] = xx [ 8 ] ; output [ 39 ] = xx [ 25 ] ; output [ 40 ] = xx [
19 ] ; output [ 41 ] = xx [ 8 ] ; output [ 42 ] = xx [ 27 ] ; output [ 43 ] =
xx [ 8 ] ; output [ 44 ] = ( xx [ 18 ] + xx [ 23 ] * xx [ 23 ] ) * xx [ 0 ] -
xx [ 13 ] ; output [ 45 ] = xx [ 8 ] ; output [ 46 ] = - xx [ 27 ] ; output [
47 ] = xx [ 8 ] ; output [ 48 ] = xx [ 19 ] ; output [ 49 ] = xx [ 34 ] - xx
[ 7 ] ; output [ 50 ] = xx [ 8 ] ; output [ 51 ] = xx [ 31 ] ; output [ 52 ]
= xx [ 35 ] ; output [ 53 ] = xx [ 8 ] ; output [ 54 ] = - xx [ 20 ] ; output
[ 55 ] = xx [ 8 ] ; output [ 56 ] = ( xx [ 32 ] + xx [ 36 ] * xx [ 36 ] ) *
xx [ 0 ] - xx [ 13 ] ; output [ 57 ] = xx [ 8 ] ; output [ 58 ] = xx [ 20 ] ;
output [ 59 ] = xx [ 8 ] ; output [ 60 ] = xx [ 35 ] ; output [ 61 ] = xx [
29 ] ; output [ 62 ] = xx [ 8 ] ; output [ 63 ] = - xx [ 37 ] ; output [ 64 ]
= xx [ 8 ] ; output [ 65 ] = ( xx [ 28 ] + xx [ 33 ] * xx [ 33 ] ) * xx [ 0 ]
- xx [ 13 ] ; output [ 66 ] = xx [ 8 ] ; output [ 67 ] = xx [ 37 ] ; output [
68 ] = xx [ 8 ] ; output [ 69 ] = xx [ 29 ] ; output [ 70 ] = xx [ 41 ] - (
xx [ 33 ] * xx [ 42 ] - xx [ 40 ] * xx [ 30 ] ) * xx [ 0 ] + xx [ 34 ] - xx [
7 ] ; output [ 71 ] = xx [ 8 ] ; output [ 72 ] = xx [ 38 ] - xx [ 0 ] * ( xx
[ 33 ] * xx [ 40 ] + xx [ 42 ] * xx [ 30 ] ) + xx [ 31 ] ; output [ 73 ] = xx
[ 48 ] ; output [ 74 ] = xx [ 8 ] ; output [ 75 ] = xx [ 45 ] ; output [ 76 ]
= xx [ 8 ] ; output [ 77 ] = ( xx [ 43 ] + xx [ 49 ] * xx [ 49 ] ) * xx [ 0 ]
- xx [ 13 ] ; output [ 78 ] = xx [ 8 ] ; output [ 79 ] = - xx [ 45 ] ; output
[ 80 ] = xx [ 8 ] ; output [ 81 ] = xx [ 48 ] ; output [ 82 ] = xx [ 46 ] ;
output [ 83 ] = xx [ 8 ] ; output [ 84 ] = xx [ 50 ] ; output [ 85 ] = xx [ 8
] ; output [ 86 ] = ( xx [ 39 ] + xx [ 47 ] * xx [ 47 ] ) * xx [ 0 ] - xx [
13 ] ; output [ 87 ] = xx [ 8 ] ; output [ 88 ] = - xx [ 50 ] ; output [ 89 ]
= xx [ 8 ] ; output [ 90 ] = xx [ 46 ] ; output [ 91 ] = xx [ 53 ] - ( xx [
52 ] * xx [ 62 ] + xx [ 54 ] * xx [ 63 ] ) * xx [ 0 ] + xx [ 1 ] - xx [ 7 ] ;
output [ 92 ] = xx [ 8 ] ; output [ 93 ] = xx [ 44 ] + xx [ 0 ] * ( xx [ 63 ]
* xx [ 62 ] - xx [ 54 ] * xx [ 52 ] ) + xx [ 57 ] ; output [ 94 ] = xx [ 64 ]
; output [ 95 ] = xx [ 8 ] ; output [ 96 ] = xx [ 17 ] ; output [ 97 ] = xx [
8 ] ; output [ 98 ] = ( xx [ 59 ] + xx [ 66 ] * xx [ 66 ] ) * xx [ 0 ] - xx [
13 ] ; output [ 99 ] = xx [ 8 ] ; output [ 100 ] = - xx [ 17 ] ; output [ 101
] = xx [ 8 ] ; output [ 102 ] = xx [ 64 ] ; output [ 103 ] = xx [ 55 ] ;
output [ 104 ] = xx [ 8 ] ; output [ 105 ] = xx [ 67 ] ; output [ 106 ] = xx
[ 8 ] ; output [ 107 ] = ( xx [ 51 ] + xx [ 60 ] * xx [ 60 ] ) * xx [ 0 ] -
xx [ 13 ] ; output [ 108 ] = xx [ 8 ] ; output [ 109 ] = - xx [ 67 ] ; output
[ 110 ] = xx [ 8 ] ; output [ 111 ] = xx [ 55 ] ; output [ 112 ] = xx [ 1 ] -
xx [ 7 ] ; output [ 113 ] = xx [ 8 ] ; output [ 114 ] = xx [ 57 ] ; output [
115 ] = xx [ 71 ] ; output [ 116 ] = xx [ 8 ] ; output [ 117 ] = - xx [ 2 ] ;
output [ 118 ] = xx [ 8 ] ; output [ 119 ] = ( xx [ 70 ] + xx [ 72 ] * xx [
72 ] ) * xx [ 0 ] - xx [ 13 ] ; output [ 120 ] = xx [ 8 ] ; output [ 121 ] =
xx [ 2 ] ; output [ 122 ] = xx [ 8 ] ; output [ 123 ] = xx [ 71 ] ; output [
124 ] = xx [ 69 ] ; output [ 125 ] = xx [ 8 ] ; output [ 126 ] = xx [ 73 ] ;
output [ 127 ] = xx [ 8 ] ; output [ 128 ] = ( xx [ 68 ] + xx [ 54 ] * xx [
54 ] ) * xx [ 0 ] - xx [ 13 ] ; output [ 129 ] = xx [ 8 ] ; output [ 130 ] =
- xx [ 73 ] ; output [ 131 ] = xx [ 8 ] ; output [ 132 ] = xx [ 69 ] ; output
[ 133 ] = xx [ 56 ] - xx [ 7 ] ; output [ 134 ] = xx [ 8 ] ; output [ 135 ] =
- xx [ 65 ] ; output [ 136 ] = xx [ 78 ] ; output [ 137 ] = xx [ 8 ] ; output
[ 138 ] = - xx [ 74 ] ; output [ 139 ] = xx [ 8 ] ; output [ 140 ] = ( xx [
77 ] + xx [ 79 ] * xx [ 79 ] ) * xx [ 0 ] - xx [ 13 ] ; output [ 141 ] = xx [
8 ] ; output [ 142 ] = xx [ 74 ] ; output [ 143 ] = xx [ 8 ] ; output [ 144 ]
= xx [ 78 ] ; output [ 145 ] = xx [ 76 ] ; output [ 146 ] = xx [ 8 ] ; output
[ 147 ] = - xx [ 80 ] ; output [ 148 ] = xx [ 8 ] ; output [ 149 ] = ( xx [
75 ] + xx [ 61 ] * xx [ 61 ] ) * xx [ 0 ] - xx [ 13 ] ; output [ 150 ] = xx [
8 ] ; output [ 151 ] = xx [ 80 ] ; output [ 152 ] = xx [ 8 ] ; output [ 153 ]
= xx [ 76 ] ; return NULL ; }
