#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "MagneticBasket_Simscape_Optimizer_capi_host.h"
#define sizeof(s) ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el) ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s) (s)
#ifndef SS_UINT64
#define SS_UINT64 19
#endif
#ifndef SS_INT64
#define SS_INT64 20
#endif
#else
#include "builtin_typeid_types.h"
#include "MagneticBasket_Simscape_Optimizer.h"
#include "MagneticBasket_Simscape_Optimizer_capi.h"
#include "MagneticBasket_Simscape_Optimizer_private.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s)               ((NULL))
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif
static const rtwCAPI_Signals rtBlockSignals [ ] = { { 0 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Gain" ) , TARGET_STRING ( "" ) , 0 , 0 , 0
, 0 , 0 } , { 1 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Subsystem1/Gain" ) , TARGET_STRING ( "" )
, 0 , 0 , 0 , 0 , 0 } , { 2 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_10_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 3 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_10_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 4 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_10_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 5 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_11_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 6 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_11_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 7 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_11_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 8 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_12_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 9 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_12_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 10 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_12_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 11 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_13_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 12 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_13_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 13 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_13_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 14 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_14_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 15 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_14_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 16 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_14_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 17 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_1_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 18 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_1_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 19 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_1_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 20 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_2_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 21 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_2_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 22 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_2_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 23 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_3_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 24 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_3_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 25 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_3_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 26 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_4_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 27 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_4_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 28 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_4_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 29 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_5_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 30 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_5_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 31 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_5_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 32 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_6_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 33 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_6_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 34 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_6_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 35 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_7_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 36 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_7_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 37 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_7_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 38 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_8_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 39 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_8_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 40 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_8_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 41 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_9_1_1"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 42 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_9_1_2"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 43 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/INPUT_9_1_3"
) , TARGET_STRING ( "" ) , 0 , 0 , 1 , 0 , 0 } , { 44 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/OUTPUT_1_0"
) , TARGET_STRING ( "" ) , 0 , 0 , 2 , 0 , 0 } , { 45 , 0 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/STATE_1" ) ,
TARGET_STRING ( "" ) , 0 , 0 , 3 , 0 , 0 } , { 46 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 47 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 48 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 49 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Gain3"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 50 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 51 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 52 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 53 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 54 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 55 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 56 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Indexing Needed"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 57 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 58 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 59 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 60 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 61 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 62 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 63 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Indexing Needed"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 64 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 65 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 66 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 67 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 68 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 69 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 70 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Indexing Needed"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 71 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 72 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 73 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 74 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 75 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 76 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 77 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Indexing Needed"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 78 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 79 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 80 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 81 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 82 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 83 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 84 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Indexing Needed"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 85 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 86 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 87 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 88 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Magnet DipoleMoment"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 4 , 0 , 1 } , { 89 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 90 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Gain2"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 91 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Indexing Needed"
) , TARGET_STRING ( "world에서 본 뱡향" ) , 0 , 0 , 5 , 0 , 0 } , { 92 , 0 ,
TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 93 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 94 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Submatrix1"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 95 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 96 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 97 , 4 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 98 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 99 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 100 , 8 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 101 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 102 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 103 , 12 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 104 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 105 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 106 , 16 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 107 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 108 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 109 , 20 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 110 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 111 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 112 , 24 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 113 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 114 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 115 , 28 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 116 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 117 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 118 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 119 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 120 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 121 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 122 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 123 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 124 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 125 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 126 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 127 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 128 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 129 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 130 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 131 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 132 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 133 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 134 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 135 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 136 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 137 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 138 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 139 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 140 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 141 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 142 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 143 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 144 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 145 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 146 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 147 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 148 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 149 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 150 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 151 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 152 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 153 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 154 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 155 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 156 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 157 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 158 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 159 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 160 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 161 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 162 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 163 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 164 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 165 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 166 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 167 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 168 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 169 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 170 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 171 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 172 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 173 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 174 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 175 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 176 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 177 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 178 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 179 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 180 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 181 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 182 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 183 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 184 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 185 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 186 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 187 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 188 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 189 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 190 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 191 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 192 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 193 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 194 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 195 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 196 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 197 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 198 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 199 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 200 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 201 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 202 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 203 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 204 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 205 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 206 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 207 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 208 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 209 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 210 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 211 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 212 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "5 * dot(m1,r_hat) * dot(m2,r_hat)" ) , 0 , 0 , 7 , 0 , 0
} , { 213 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Gain3"
) , TARGET_STRING ( "3*myu_0/(4*pi*r^4)" ) , 0 , 0 , 7 , 0 , 0 } , { 214 , 0
, TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 215 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Product1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 216 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Product4"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 217 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Product6"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 218 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 219 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Add2"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 220 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Gain"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 221 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 222 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Transpose"
) , TARGET_STRING ( "" ) , 0 , 0 , 8 , 0 , 0 } , { 223 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 224 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/MatrixMultiply1"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 225 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 4 , 0 , 0 } , { 226 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Add"
) , TARGET_STRING ( "" ) , 0 , 0 , 6 , 0 , 0 } , { 227 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/PS-Simulink Converter1/EVAL_KEY/RESHAPE"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 228 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 229 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 230 , 1 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 231 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 232 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 233 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 234 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 235 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 236 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 237 , 2 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 238 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 239 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 240 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 241 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 242 , 3 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 243 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 244 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 245 , 5 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 246 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 247 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 248 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 249 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 250 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 251 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 252 , 6 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 253 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 254 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 255 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 256 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 257 , 7 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 258 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 259 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 260 , 9 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 261 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 262 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 263 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 264 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 265 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 266 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 267 , 10 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 268 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 269 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 270 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 271 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 272 , 11 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 273 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 274 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 275 , 13 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 276 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 277 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 278 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 279 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 280 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 281 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 282 , 14 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 283 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 284 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 285 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 286 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 287 , 15 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 288 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 289 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 290 , 17 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 291 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 292 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 293 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 294 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 295 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 296 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 297 , 18 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 298 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 299 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 300 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 301 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 302 , 19 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 303 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 304 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 305 , 21 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 306 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 307 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 308 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 309 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 310 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 311 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 312 , 22 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 313 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 314 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 315 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 316 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 317 , 23 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 318 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 319 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 320 , 25 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 321 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 322 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 323 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Element Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 10 , 0 , 0 } , { 324 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Cross Product/Sum"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 325 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 326 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 327 , 26 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 328 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Sum of Elements"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 329 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Switch"
) , TARGET_STRING ( "" ) , 0 , 0 , 9 , 0 , 0 } , { 330 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Math Function1"
) , TARGET_STRING ( "" ) , 0 , 0 , 7 , 0 , 0 } , { 331 , 0 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Divide"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 332 , 27 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Product"
) , TARGET_STRING ( "" ) , 0 , 0 , 5 , 0 , 0 } , { 0 , 0 , ( NULL ) , ( NULL
) , 0 , 0 , 0 , 0 , 0 } } ; static const rtwCAPI_BlockParameters
rtBlockParameters [ ] = { { 333 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Gain" ) , TARGET_STRING ( "Gain" ) , 0 , 7
, 0 } , { 334 , TARGET_STRING (
"MagneticBasket_Simscape_Optimizer/Subsystem1/Gain" ) , TARGET_STRING (
"Gain" ) , 0 , 7 , 0 } , { 335 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 336 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 337 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 338 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 339 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 340 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 341 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 342 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 343 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 344 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 345 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 346 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 347 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 348 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Magnet DipoleMoment"
) , TARGET_STRING ( "Value" ) , 0 , 11 , 0 } , { 349 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 350 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 351 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 352 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 353 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 354 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 355 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 356 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 357 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 358 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 359 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 360 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 361 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 362 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 363 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 364 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 365 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 366 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 367 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 368 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 369 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 370 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 371 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 372 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 373 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 374 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 375 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 376 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 377 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 378 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 379 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 380 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 381 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 382 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 383 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 384 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 385 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 386 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 387 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 388 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 389 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 390 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 391 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 392 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 393 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 394 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 395 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 396 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 397 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 398 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 399 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 400 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 401 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 402 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 403 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 404 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 405 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 406 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Gain"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 407 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 408 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1"
) , TARGET_STRING ( "maxzero" ) , 0 , 7 , 0 } , { 409 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 410 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Constant1"
) , TARGET_STRING ( "Value" ) , 0 , 6 , 0 } , { 411 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Gain1"
) , TARGET_STRING ( "Gain" ) , 0 , 7 , 0 } , { 412 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 413 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 414 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator1/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 415 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 416 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 417 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator2/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 418 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 419 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 420 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator3/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 421 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 422 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 423 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator4/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 424 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 425 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 426 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator5/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 427 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 428 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 429 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator6/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 430 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Force/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 431 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 432 , TARGET_STRING (
 "MagneticBasket_Simscape_Optimizer/Subsystem1/Magnetic Torque//Force Generator7/Calculate Magnetic Force//Torque/Magnetic Torque/Normalize Vector1/Constant"
) , TARGET_STRING ( "Value" ) , 0 , 7 , 0 } , { 0 , ( NULL ) , ( NULL ) , 0 ,
0 , 0 } } ; static int_T rt_LoggedStateIdxList [ ] = { - 1 } ; static const
rtwCAPI_Signals rtRootInputs [ ] = { { 0 , 0 , ( NULL ) , ( NULL ) , 0 , 0 ,
0 , 0 , 0 } } ; static const rtwCAPI_Signals rtRootOutputs [ ] = { { 0 , 0 ,
( NULL ) , ( NULL ) , 0 , 0 , 0 , 0 , 0 } } ; static const
rtwCAPI_ModelParameters rtModelParameters [ ] = { { 433 , TARGET_STRING (
"Fixed" ) , 1 , 7 , 0 } , { 434 , TARGET_STRING ( "SimParams" ) , 2 , 7 , 0 }
, { 435 , TARGET_STRING ( "x" ) , 0 , 14 , 0 } , { 0 , ( NULL ) , 0 , 0 , 0 }
} ;
#ifndef HOST_CAPI_BUILD
static void * rtDataAddrMap [ ] = { & rtB . o3cnonqf2y [ 0 ] , & rtB .
gdleddk2bu [ 0 ] , & rtB . g1yo2p3pxt [ 0 ] , & rtB . az4eo3zh5o [ 0 ] , &
rtB . ch1tbnhyon [ 0 ] , & rtB . eoseoj2f3b [ 0 ] , & rtB . aigpe1li2p [ 0 ]
, & rtB . pbn2nzk2k3 [ 0 ] , & rtB . hffhwq1imf [ 0 ] , & rtB . gbxuu4wtte [
0 ] , & rtB . k5nlme3kns [ 0 ] , & rtB . j52cfj2d1w [ 0 ] , & rtB .
ftvz2ekipg [ 0 ] , & rtB . awgskjjmsb [ 0 ] , & rtB . fl4rvtkarn [ 0 ] , &
rtB . jp22vgjwwm [ 0 ] , & rtB . pzpjv5s5n2 [ 0 ] , & rtB . c513o5pg5r [ 0 ]
, & rtB . pfmi5qvzhf [ 0 ] , & rtB . gi35bo1nw4 [ 0 ] , & rtB . iii5xyl5wn [
0 ] , & rtB . mjqwg1vruy [ 0 ] , & rtB . eo0yhyzcnn [ 0 ] , & rtB .
c5dvbg4axm [ 0 ] , & rtB . bxbu5xamaf [ 0 ] , & rtB . foqtyz1zpt [ 0 ] , &
rtB . gy3u0julyl [ 0 ] , & rtB . hytxiurw0s [ 0 ] , & rtB . nubcraskrw [ 0 ]
, & rtB . lp0cxari03 [ 0 ] , & rtB . phx5lsta5w [ 0 ] , & rtB . hoap0zr3eu [
0 ] , & rtB . enp5f3s002 [ 0 ] , & rtB . nsl31ekm0f [ 0 ] , & rtB .
hrqvwki2jl [ 0 ] , & rtB . dypcby02yn [ 0 ] , & rtB . a100glwqkz [ 0 ] , &
rtB . c3twv1nzdj [ 0 ] , & rtB . iurh4h1cla [ 0 ] , & rtB . har0a2wrmi [ 0 ]
, & rtB . ib2ieuf2an [ 0 ] , & rtB . ojpovrodvr [ 0 ] , & rtB . fcgmpow5uz [
0 ] , & rtB . jylntbxfpj [ 0 ] , & rtB . nfyrk0gene [ 0 ] , & rtB .
aay54cmgia [ 0 ] , & rtB . nsmqj1yv3h [ 0 ] , & rtB . j4tq1cqvap [ 0 ] , &
rtB . i3fsmpf0l4 [ 0 ] , & rtB . ecnsyysan1 [ 0 ] , & rtB . dg2tz20qh3 [ 0 ]
, & rtB . ciemy5f0t3 [ 0 ] , & rtB . bjyrxtyjst [ 0 ] , & rtB . mdrxplc5qx [
0 ] , & rtB . dakqbtmvna [ 0 ] , & rtB . cwwkrda1qf [ 0 ] , & rtB .
hjkncigkim [ 0 ] , & rtB . leoubbnbua [ 0 ] , & rtB . iwp0nxhxfg [ 0 ] , &
rtB . mrse3gkxz4 [ 0 ] , & rtB . idcqc0uibn [ 0 ] , & rtB . dndiwk2b1s [ 0 ]
, & rtB . nalsmrc4um [ 0 ] , & rtB . ol0mge5uht [ 0 ] , & rtB . igdgu1b2hk [
0 ] , & rtB . bcueeadxdi [ 0 ] , & rtB . j54ujtxdmj [ 0 ] , & rtB .
lkwtwvkutl [ 0 ] , & rtB . cbmxpbd2iu [ 0 ] , & rtB . b0dwpzqyrm [ 0 ] , &
rtB . dw12uizmd2 [ 0 ] , & rtB . cpzvyphms1 [ 0 ] , & rtB . ley50ghvfj [ 0 ]
, & rtB . gvavv3jbs0 [ 0 ] , & rtB . dunjuwexfg [ 0 ] , & rtB . l01qgkxm4q [
0 ] , & rtB . lxihdocxga [ 0 ] , & rtB . f34tk41pk1 [ 0 ] , & rtB .
obk5xjqd1g [ 0 ] , & rtB . emdpdwkzg3 [ 0 ] , & rtB . l1wqaovs4e [ 0 ] , &
rtB . oimhqzeh2t [ 0 ] , & rtB . g3vpjafupq [ 0 ] , & rtB . j4fldanieq [ 0 ]
, & rtB . fzbgmnan4p [ 0 ] , & rtB . fzqka3pwrj [ 0 ] , & rtB . e3buvmif1c [
0 ] , & rtB . babxk0cwue [ 0 ] , & rtB . bcudxny0c0 [ 0 ] , & rtB .
j2tnneegmr [ 0 ] , & rtB . o2znllyvoq [ 0 ] , & rtB . gmivsxqevf [ 0 ] , &
rtB . hnzvyzmoy4 [ 0 ] , & rtB . exqp2uohe5 [ 0 ] , & rtB . fmqgelgm0i [ 0 ]
, & rtB . oywwaak55z , & rtB . n3uyejqyiv [ 0 ] , & rtB . itxt4ytwah [ 0 ] ,
& rtB . jknqpcyts5 , & rtB . fo44k43jik [ 0 ] , & rtB . p5hpxiz5u3 [ 0 ] , &
rtB . n2cw123fsj , & rtB . gtc4x0vmrj [ 0 ] , & rtB . inhs5zftab [ 0 ] , &
rtB . fh04srz5g0 , & rtB . mkzi43umox [ 0 ] , & rtB . itx5qkpva4 [ 0 ] , &
rtB . d4it13dar5 , & rtB . dl4tl2eodk [ 0 ] , & rtB . jltrdwyzmw [ 0 ] , &
rtB . efmfmwthqv , & rtB . ivaqfjshfj [ 0 ] , & rtB . l02fn5ndrj [ 0 ] , &
rtB . jgkbpstejl , & rtB . evftq3sk0q [ 0 ] , & rtB . ibott2adsg [ 0 ] , &
rtB . jmgra3bzio , & rtB . p1t00nt2c4 , & rtB . dm0nsvocoi [ 0 ] , & rtB .
gfwkvsxkxh [ 0 ] , & rtB . pi3p3kites [ 0 ] , & rtB . dxk5dr2guu [ 0 ] , &
rtB . pzoqxgrplf [ 0 ] , & rtB . n43m55vtpu , & rtB . mnl5wynqkw , & rtB .
j2caeelflf [ 0 ] , & rtB . axlew0gacz [ 0 ] , & rtB . ffk0uso1jv [ 0 ] , &
rtB . nhpvoyt4rs [ 0 ] , & rtB . n5zmerun4k [ 0 ] , & rtB . gzdxak0rgb [ 0 ]
, & rtB . d11srgkxv1 [ 0 ] , & rtB . nj4zxwlnaz , & rtB . ojzjawyp2p , & rtB
. mswvoafvyk [ 0 ] , & rtB . ljdiomnr0e [ 0 ] , & rtB . eecvbm2tyj [ 0 ] , &
rtB . oexwmjcbqo [ 0 ] , & rtB . ga3a4xdbbc [ 0 ] , & rtB . b4eypdlrbg , &
rtB . ollqza00rh , & rtB . cviiadazje [ 0 ] , & rtB . a5iwna5ofw [ 0 ] , &
rtB . n5qkyzyc0p [ 0 ] , & rtB . p4bi154apu [ 0 ] , & rtB . kypr3ufomq [ 0 ]
, & rtB . hrdkh50xl2 [ 0 ] , & rtB . bv20qr1sq1 [ 0 ] , & rtB . hv30h5tsxp ,
& rtB . gmo4itlisv , & rtB . pl32irduyc [ 0 ] , & rtB . l0eiwonrmi [ 0 ] , &
rtB . nexxpkogi1 [ 0 ] , & rtB . lgifpnih0e [ 0 ] , & rtB . hrtkouuc3i [ 0 ]
, & rtB . iqlb144ibi , & rtB . fd02dku0jj , & rtB . o0sl3javac [ 0 ] , & rtB
. pctih4a2gv [ 0 ] , & rtB . d15mjte2mm [ 0 ] , & rtB . kucavif4r4 [ 0 ] , &
rtB . j3daqzj5si [ 0 ] , & rtB . lk2g50sct1 [ 0 ] , & rtB . fusdzrgc4y [ 0 ]
, & rtB . encqrxuu3y , & rtB . pxg5b2slaw , & rtB . nbrrv2d3zg [ 0 ] , & rtB
. eegu2ofo43 [ 0 ] , & rtB . h4h1ivsmcr [ 0 ] , & rtB . cuxk4j3qys [ 0 ] , &
rtB . hacw1qj3aj [ 0 ] , & rtB . bdozk1fpny , & rtB . f1arg5nwa3 , & rtB .
i1lnvldtbv [ 0 ] , & rtB . c402o0md5n [ 0 ] , & rtB . hwppgfjqwy [ 0 ] , &
rtB . ort4y5zlaf [ 0 ] , & rtB . eckq3h2tzb [ 0 ] , & rtB . ahw1dydhgo [ 0 ]
, & rtB . bdswxd3b20 [ 0 ] , & rtB . j1ydxe2usi , & rtB . ezvzwoavxd , & rtB
. lyhkcvfx5d [ 0 ] , & rtB . mw4hhpohry [ 0 ] , & rtB . msqcnwfh2s [ 0 ] , &
rtB . enoycngmy3 [ 0 ] , & rtB . bhsuwz2d2l [ 0 ] , & rtB . ceb1rzrjkf , &
rtB . lffwcdbemm , & rtB . bqjv3qidt0 [ 0 ] , & rtB . odlamrju54 [ 0 ] , &
rtB . p3ru4cddgi [ 0 ] , & rtB . gxha2ko2vz [ 0 ] , & rtB . chig3jpwbn [ 0 ]
, & rtB . lkoo4fzueb [ 0 ] , & rtB . kup2makpml [ 0 ] , & rtB . cab5adcq1r ,
& rtB . cdgxtvlobx , & rtB . f4nnctaymu [ 0 ] , & rtB . m4bi5zptkh [ 0 ] , &
rtB . mrdivqgmvm [ 0 ] , & rtB . h1zdr1aani [ 0 ] , & rtB . mcasp4wuog [ 0 ]
, & rtB . avcgwoxwhd , & rtB . k3vxq1m344 , & rtB . aaw4gtn4yc [ 0 ] , & rtB
. d2uj0bf0jq [ 0 ] , & rtB . joimmmlrrj [ 0 ] , & rtB . le0lvwzojd [ 0 ] , &
rtB . e3hmodlujb [ 0 ] , & rtB . jnz14mmcr3 [ 0 ] , & rtB . oflmpiaez1 [ 0 ]
, & rtB . aizglk3yjb , & rtB . bbtvepo4p4 , & rtB . azxdea4e00 [ 0 ] , & rtB
. jqxy4sgz1r [ 0 ] , & rtB . gdhtjuaiqp [ 0 ] , & rtB . lperucgciu [ 0 ] , &
rtB . afxh4mt3gu [ 0 ] , & rtB . avjauevoe4 , & rtB . j2h4dmgm1o , & rtB .
lrkvtt341w [ 0 ] , & rtB . li0pbutf3d [ 0 ] , & rtB . nzexd0tfyn [ 0 ] , &
rtB . nwco0rdsmq [ 0 ] , & rtB . kdzk5onxbh [ 0 ] , & rtB . kptzapguuy [ 0 ]
, & rtB . cxwqk2amda [ 0 ] , & rtB . j0nzmiwajs , & rtB . kerf1yzlz0 [ 0 ] ,
& rtB . nuca2ykly0 [ 0 ] , & rtB . lkk1p3hkce , & rtB . p4b33uxmxl [ 0 ] , &
rtB . i210rcby1e [ 0 ] , & rtB . njlotg3eyk [ 0 ] , & rtB . i1wwep0quu , &
rtB . dvlkqx1wrk [ 0 ] , & rtB . kdrism1ogn [ 0 ] , & rtB . jd0epngxde , &
rtB . jytdis5y4s [ 0 ] , & rtB . d4sl1kvckp , & rtB . k0e5tb5t4r [ 0 ] , &
rtB . nl3crkyc5i [ 0 ] , & rtB . ltesgv5ex5 , & rtB . fi20a5zfv5 [ 0 ] , &
rtB . npob1uuqas [ 0 ] , & rtB . j22ln0sz4g , & rtB . e4k3bfcn1z [ 0 ] , &
rtB . ljklcjruel [ 0 ] , & rtB . giyavk0qxh [ 0 ] , & rtB . eded4se0e2 , &
rtB . huyyr23tvs [ 0 ] , & rtB . jqlckttaam [ 0 ] , & rtB . eb5u3a2x55 , &
rtB . epqt4glrc2 [ 0 ] , & rtB . bkxjvl2wwl , & rtB . jcdegmljub [ 0 ] , &
rtB . h2ica3jpzg [ 0 ] , & rtB . ij2rjc05nz , & rtB . ndcqoyp0vx [ 0 ] , &
rtB . o4wni3515w [ 0 ] , & rtB . gmj04jp0xd , & rtB . pyxwiqcuay [ 0 ] , &
rtB . e30x0ooogi [ 0 ] , & rtB . dgv155wtql [ 0 ] , & rtB . lvldddhyjb , &
rtB . f2e0m2twoo [ 0 ] , & rtB . cschdiniut [ 0 ] , & rtB . lntlwwhs4g , &
rtB . j1izfy2apg [ 0 ] , & rtB . idwieepnv4 , & rtB . h2vbgw30tf [ 0 ] , &
rtB . hrudfjbz3i [ 0 ] , & rtB . my44vuapi0 , & rtB . cg3eiq4i4s [ 0 ] , &
rtB . fvtfdqea2r [ 0 ] , & rtB . fbdf2tllgu , & rtB . k3faf5atq0 [ 0 ] , &
rtB . bsaqv0u04c [ 0 ] , & rtB . alzuotumkg [ 0 ] , & rtB . hlnjhscgik , &
rtB . cafqksoz00 [ 0 ] , & rtB . mdxnlbsfvm [ 0 ] , & rtB . oc1jrzdnhe , &
rtB . fkzq14ohjb [ 0 ] , & rtB . jqzltcnumt , & rtB . etawd0whn4 [ 0 ] , &
rtB . e4tk4xaa2e [ 0 ] , & rtB . ntt11gc034 , & rtB . bpewjnssvt [ 0 ] , &
rtB . p5ihfgrspp [ 0 ] , & rtB . nfak5a1vn4 , & rtB . hvd2sf2hag [ 0 ] , &
rtB . fcsl1vvaio [ 0 ] , & rtB . pu15hoqemb [ 0 ] , & rtB . b4b10et3s3 , &
rtB . dv2svmenmy [ 0 ] , & rtB . lebemiabnp [ 0 ] , & rtB . p00fanezjt , &
rtB . flzta24bjr [ 0 ] , & rtB . mp054zeo2w , & rtB . pfosbwya4i [ 0 ] , &
rtB . fc15cwanga [ 0 ] , & rtB . ieyrepmb05 , & rtB . favu2vq25u [ 0 ] , &
rtB . bjhxaihh05 [ 0 ] , & rtB . lbmtuwi01h , & rtB . lyjnnzstpa [ 0 ] , &
rtB . b3yujzbge2 [ 0 ] , & rtB . gf5n1xgexk [ 0 ] , & rtB . cjrpkdpszw , &
rtB . ksrswfmi3x [ 0 ] , & rtB . d1zqfcz4zo [ 0 ] , & rtB . egvr4ozujh , &
rtB . logpc2hca5 [ 0 ] , & rtB . btylsma43j , & rtB . my424itvki [ 0 ] , &
rtB . emllridil0 [ 0 ] , & rtB . c00wpmk13d , & rtB . gg2oxfw35t [ 0 ] , &
rtB . joumqp0zpt [ 0 ] , & rtB . kq5inrtfnq , & rtB . huacygyj0o [ 0 ] , &
rtB . oiy31a51av [ 0 ] , & rtB . athb42mz2d [ 0 ] , & rtB . ifl2ongx1d , &
rtB . glnk1qsnt1 [ 0 ] , & rtB . eejq0hbv30 [ 0 ] , & rtB . gtu1mufolx , &
rtB . a3qlbfbxxo [ 0 ] , & rtB . e0gxuicxzo , & rtB . kbih5351xy [ 0 ] , &
rtB . ovchpiv2ws [ 0 ] , & rtP . Gain_Gain_iby0aikld1 , & rtP . Gain_Gain , &
rtP . NormalizeVector_maxzero , & rtP . MagnetDipoleMoment_Value [ 0 ] , &
rtP . NormalizeVector_maxzero_hach1vuaw5 , & rtP .
MagnetDipoleMoment_Value_hus3icj2zx [ 0 ] , & rtP .
NormalizeVector_maxzero_jvlgwkux2b , & rtP .
MagnetDipoleMoment_Value_m4h5ot0d2d [ 0 ] , & rtP .
NormalizeVector_maxzero_levzkkbbsj , & rtP .
MagnetDipoleMoment_Value_am3fjqdljc [ 0 ] , & rtP .
NormalizeVector_maxzero_anlse3xi4h , & rtP .
MagnetDipoleMoment_Value_hqiqmgkprz [ 0 ] , & rtP .
NormalizeVector_maxzero_ap3h53mgxx , & rtP .
MagnetDipoleMoment_Value_clffg1mnnq [ 0 ] , & rtP .
NormalizeVector_maxzero_fqwaltgljo , & rtP .
MagnetDipoleMoment_Value_af4j3ouf2p [ 0 ] , & rtP . Constant_Value_lu1rr1djht
, & rtP . Constant_Value_kebaex1duz , & rtP . Constant_Value_mqmlcee2pu , &
rtP . Constant_Value_iinhm22lbs , & rtP . Constant_Value_e34xedeirm , & rtP .
Constant_Value_gimcsalf5s , & rtP . Constant_Value_nfejvnig4t , & rtP .
NormalizeVector_maxzero_jpxndxsk21 , & rtP . Constant_Value , & rtP .
Gain_Gain_o3ussx2pey , & rtP . NormalizeVector_maxzero_ab1bjfgz5c , & rtP .
NormalizeVector1_maxzero , & rtP . Constant_Value_ffpsqzlbdy , & rtP .
Constant1_Value [ 0 ] , & rtP . Gain1_Gain , & rtP .
NormalizeVector_maxzero_pjti1ch2qi , & rtP . Constant_Value_hrhqyti5xa , &
rtP . Gain_Gain_ch2juginai , & rtP . NormalizeVector_maxzero_nzh4nf0nvr , &
rtP . NormalizeVector1_maxzero_hfyr12og40 , & rtP . Constant_Value_ccdbgnkv2m
, & rtP . Constant1_Value_mat3dcwilq [ 0 ] , & rtP . Gain1_Gain_inzuw1ftmc ,
& rtP . NormalizeVector_maxzero_jceokergvn , & rtP .
Constant_Value_crcwx2inuv , & rtP . Gain_Gain_cnz04octod , & rtP .
NormalizeVector_maxzero_k3zfgywmre , & rtP .
NormalizeVector1_maxzero_edxcmuydlr , & rtP . Constant_Value_pxe2jyk15f , &
rtP . Constant1_Value_obh4dra0ft [ 0 ] , & rtP . Gain1_Gain_iijf52nd0p , &
rtP . NormalizeVector_maxzero_hbfzo5qsh2 , & rtP . Constant_Value_mulpe1a5eh
, & rtP . Gain_Gain_mxfyehzrbo , & rtP . NormalizeVector_maxzero_kt51ehrppa ,
& rtP . NormalizeVector1_maxzero_os1qknjcxt , & rtP .
Constant_Value_kyzvzubrmq , & rtP . Constant1_Value_htaifltrf1 [ 0 ] , & rtP
. Gain1_Gain_kgqmwgbxmo , & rtP . NormalizeVector_maxzero_nfiirapz00 , & rtP
. Constant_Value_o54qqmia1e , & rtP . Gain_Gain_ch2qj2nra0 , & rtP .
NormalizeVector_maxzero_dfin15pekf , & rtP .
NormalizeVector1_maxzero_pbowfqqdx5 , & rtP . Constant_Value_c50rewkplm , &
rtP . Constant1_Value_ajigyq4ohn [ 0 ] , & rtP . Gain1_Gain_pklggewqk2 , &
rtP . NormalizeVector_maxzero_ihmtmp1vz1 , & rtP . Constant_Value_l5p2q1mni4
, & rtP . Gain_Gain_kzgbyxe2u1 , & rtP . NormalizeVector_maxzero_kytveykcqy ,
& rtP . NormalizeVector1_maxzero_byir5pwee4 , & rtP .
Constant_Value_pik4yifivl , & rtP . Constant1_Value_dagdlneepy [ 0 ] , & rtP
. Gain1_Gain_mt5ivn2ky5 , & rtP . NormalizeVector_maxzero_dnlmszn2ri , & rtP
. Constant_Value_h2tkrf03re , & rtP . Gain_Gain_gnnygtgrhf , & rtP .
NormalizeVector_maxzero_hkbi1qbiut , & rtP .
NormalizeVector1_maxzero_j5bp1hk5bd , & rtP . Constant_Value_nu3mnt2521 , &
rtP . Constant1_Value_k2adofcs3i [ 0 ] , & rtP . Gain1_Gain_fzh3lvxzcm , &
rtP . Constant_Value_foidw3bu40 , & rtP . Constant_Value_aldt0rz0ce , & rtP .
Constant_Value_hsgo25r4sk , & rtP . Constant_Value_bgbfcrykss , & rtP .
Constant_Value_ag31g0q0wg , & rtP . Constant_Value_lpo3xrcqjh , & rtP .
Constant_Value_p4qexwk3wd , & rtP . Constant_Value_a4fqp5ywdt , & rtP .
Constant_Value_lwhq5q2fck , & rtP . Constant_Value_f4bgccfy2v , & rtP .
Constant_Value_hooldt5mbh , & rtP . Constant_Value_fyrymlyvwp , & rtP .
Constant_Value_fkwr40m5ih , & rtP . Constant_Value_gxy2qdgtgj , & rtP .
Constant_Value_oqtflizduh , & rtP . Constant_Value_cuvlenn5ze , & rtP .
Constant_Value_ndmgh4xmik , & rtP . Constant_Value_kqizw1abe3 , & rtP .
Constant_Value_oc3ljvgxhv , & rtP . Constant_Value_m0tx3vya3h , & rtP .
Constant_Value_jgubgdhktv , & rtP . Fixed , & rtP . SimParams , & rtP . x [ 0
] , } ; static int32_T * rtVarDimsAddrMap [ ] = { ( NULL ) } ;
#endif
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap [ ] = { { "double" ,
"real_T" , 0 , 0 , sizeof ( real_T ) , ( uint8_T ) SS_DOUBLE , 0 , 0 , 0 } ,
{ "struct" , "struct_C8RyKn0xmbw5jBmaos4hYE" , 12 , 1 , sizeof (
struct_C8RyKn0xmbw5jBmaos4hYE ) , ( uint8_T ) SS_STRUCT , 0 , 0 , 0 } , {
"struct" , "struct_0rpQxwwyqCy2sr7lBYKBtB" , 10 , 13 , sizeof (
struct_0rpQxwwyqCy2sr7lBYKBtB ) , ( uint8_T ) SS_STRUCT , 0 , 0 , 0 } } ;
#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif
static TARGET_CONST rtwCAPI_ElementMap rtElementMap [ ] = { { ( NULL ) , 0 ,
0 , 0 , 0 } , { "PlateLength" , rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE ,
PlateLength ) , 0 , 12 , 0 } , { "PlateWidth" , rt_offsetof (
struct_C8RyKn0xmbw5jBmaos4hYE , PlateWidth ) , 0 , 12 , 0 } , { "PlateVolume"
, rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE , PlateVolume ) , 0 , 12 , 0 }
, { "num_links" , rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE , num_links ) ,
0 , 12 , 0 } , { "MagnetHomeConfig" , rt_offsetof (
struct_C8RyKn0xmbw5jBmaos4hYE , MagnetHomeConfig ) , 0 , 13 , 0 } , {
"MagnetDirection" , rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE ,
MagnetDirection ) , 0 , 13 , 0 } , { "Br" , rt_offsetof (
struct_C8RyKn0xmbw5jBmaos4hYE , Br ) , 0 , 12 , 0 } , { "Volume" ,
rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE , Volume ) , 0 , 12 , 0 } , {
"myu" , rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE , myu ) , 0 , 12 , 0 } ,
{ "MagnetMagnitude" , rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE ,
MagnetMagnitude ) , 0 , 12 , 0 } , { "MangetDipoleMoment" , rt_offsetof (
struct_C8RyKn0xmbw5jBmaos4hYE , MangetDipoleMoment ) , 0 , 13 , 0 } , { "Kt"
, rt_offsetof ( struct_C8RyKn0xmbw5jBmaos4hYE , Kt ) , 0 , 12 , 0 } , {
"gravity_onoff" , rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB , gravity_onoff
) , 0 , 12 , 0 } , { "stiffness" , rt_offsetof (
struct_0rpQxwwyqCy2sr7lBYKBtB , stiffness ) , 0 , 12 , 0 } , { "damping" ,
rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB , damping ) , 0 , 12 , 0 } , {
"density" , rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB , density ) , 0 , 12
, 0 } , { "Opacity" , rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB , Opacity )
, 0 , 12 , 0 } , { "arrow_Opacity" , rt_offsetof (
struct_0rpQxwwyqCy2sr7lBYKBtB , arrow_Opacity ) , 0 , 12 , 0 } , {
"arrow_scale" , rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB , arrow_scale ) ,
0 , 12 , 0 } , { "force_onoff" , rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB
, force_onoff ) , 0 , 12 , 0 } , { "torque_onoff" , rt_offsetof (
struct_0rpQxwwyqCy2sr7lBYKBtB , torque_onoff ) , 0 , 12 , 0 } , {
"SimulationTime" , rt_offsetof ( struct_0rpQxwwyqCy2sr7lBYKBtB ,
SimulationTime ) , 0 , 12 , 0 } } ; static const rtwCAPI_DimensionMap
rtDimensionMap [ ] = { { rtwCAPI_VECTOR , 0 , 2 , 0 } , { rtwCAPI_VECTOR , 2
, 2 , 0 } , { rtwCAPI_VECTOR , 4 , 2 , 0 } , { rtwCAPI_VECTOR , 6 , 2 , 0 } ,
{ rtwCAPI_VECTOR , 8 , 2 , 0 } , { rtwCAPI_MATRIX_COL_MAJOR , 8 , 2 , 0 } , {
rtwCAPI_MATRIX_COL_MAJOR , 10 , 2 , 0 } , { rtwCAPI_SCALAR , 12 , 2 , 0 } , {
rtwCAPI_MATRIX_COL_MAJOR , 14 , 2 , 0 } , { rtwCAPI_MATRIX_COL_MAJOR , 2 , 2
, 0 } , { rtwCAPI_MATRIX_COL_MAJOR , 16 , 2 , 0 } , { rtwCAPI_VECTOR , 14 , 2
, 0 } , { rtwCAPI_MATRIX_COL_MAJOR , 12 , 2 , 0 } , {
rtwCAPI_MATRIX_COL_MAJOR , 18 , 2 , 0 } , { rtwCAPI_VECTOR , 20 , 2 , 0 } } ;
static const uint_T rtDimensionArray [ ] = { 7 , 1 , 4 , 1 , 154 , 1 , 18 , 1
, 3 , 1 , 3 , 3 , 1 , 1 , 1 , 3 , 6 , 1 , 1 , 2 , 1 , 14 } ; static const
real_T rtcapiStoredFloats [ ] = { 0.0 } ; static const rtwCAPI_FixPtMap
rtFixPtMap [ ] = { { ( NULL ) , ( NULL ) , rtwCAPI_FIX_RESERVED , 0 , 0 , (
boolean_T ) 0 } , } ; static const rtwCAPI_SampleTimeMap rtSampleTimeMap [ ]
= { { ( const void * ) & rtcapiStoredFloats [ 0 ] , ( const void * ) &
rtcapiStoredFloats [ 0 ] , ( int8_T ) 0 , ( uint8_T ) 0 } , { ( NULL ) , (
NULL ) , 1 , 0 } } ; static rtwCAPI_ModelMappingStaticInfo mmiStatic = { {
rtBlockSignals , 333 , rtRootInputs , 0 , rtRootOutputs , 0 } , {
rtBlockParameters , 100 , rtModelParameters , 3 } , { ( NULL ) , 0 } , {
rtDataTypeMap , rtDimensionMap , rtFixPtMap , rtElementMap , rtSampleTimeMap
, rtDimensionArray } , "float" , { 343361676U , 2573916448U , 664120902U ,
3813395425U } , ( NULL ) , 0 , ( boolean_T ) 0 , rt_LoggedStateIdxList } ;
const rtwCAPI_ModelMappingStaticInfo *
MagneticBasket_Simscape_Optimizer_GetCAPIStaticMap ( void ) { return &
mmiStatic ; }
#ifndef HOST_CAPI_BUILD
void MagneticBasket_Simscape_Optimizer_InitializeDataMapInfo ( void ) {
rtwCAPI_SetVersion ( ( * rt_dataMapInfoPtr ) . mmi , 1 ) ;
rtwCAPI_SetStaticMap ( ( * rt_dataMapInfoPtr ) . mmi , & mmiStatic ) ;
rtwCAPI_SetLoggingStaticMap ( ( * rt_dataMapInfoPtr ) . mmi , ( NULL ) ) ;
rtwCAPI_SetDataAddressMap ( ( * rt_dataMapInfoPtr ) . mmi , rtDataAddrMap ) ;
rtwCAPI_SetVarDimsAddressMap ( ( * rt_dataMapInfoPtr ) . mmi ,
rtVarDimsAddrMap ) ; rtwCAPI_SetInstanceLoggingInfo ( ( * rt_dataMapInfoPtr )
. mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArray ( ( * rt_dataMapInfoPtr ) . mmi
, ( NULL ) ) ; rtwCAPI_SetChildMMIArrayLen ( ( * rt_dataMapInfoPtr ) . mmi ,
0 ) ; }
#else
#ifdef __cplusplus
extern "C" {
#endif
void MagneticBasket_Simscape_Optimizer_host_InitializeDataMapInfo (
MagneticBasket_Simscape_Optimizer_host_DataMapInfo_T * dataMap , const char *
path ) { rtwCAPI_SetVersion ( dataMap -> mmi , 1 ) ; rtwCAPI_SetStaticMap (
dataMap -> mmi , & mmiStatic ) ; rtwCAPI_SetDataAddressMap ( dataMap -> mmi ,
( NULL ) ) ; rtwCAPI_SetVarDimsAddressMap ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetPath ( dataMap -> mmi , path ) ; rtwCAPI_SetFullPath ( dataMap ->
mmi , ( NULL ) ) ; rtwCAPI_SetChildMMIArray ( dataMap -> mmi , ( NULL ) ) ;
rtwCAPI_SetChildMMIArrayLen ( dataMap -> mmi , 0 ) ; }
#ifdef __cplusplus
}
#endif
#endif
