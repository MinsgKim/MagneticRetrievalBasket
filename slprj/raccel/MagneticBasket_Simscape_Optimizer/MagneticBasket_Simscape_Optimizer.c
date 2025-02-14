#include "MagneticBasket_Simscape_Optimizer.h"
#include "rtwtypes.h"
#include <string.h>
#include "mwmathutil.h"
#include <stddef.h>
#include "MagneticBasket_Simscape_Optimizer_private.h"
#include "rt_logging_mmi.h"
#include "MagneticBasket_Simscape_Optimizer_capi.h"
#include "MagneticBasket_Simscape_Optimizer_dt.h"
extern void * CreateDiagnosticAsVoidPtr_wrapper ( const char * id , int nargs
, ... ) ; extern ssExecutionInfo gblExecutionInfo ; RTWExtModeInfo *
gblRTWExtModeInfo = NULL ; void raccelForceExtModeShutdown ( boolean_T
extModeStartPktReceived ) { if ( ! extModeStartPktReceived ) { boolean_T
stopRequested = false ; rtExtModeWaitForStartPkt ( gblRTWExtModeInfo , 1 , &
stopRequested ) ; } rtExtModeShutdown ( 1 ) ; }
#include "slsv_diagnostic_codegen_c_api.h"
#include "slsa_sim_engine.h"
#ifdef RSIM_WITH_SOLVER_MULTITASKING
boolean_T gbl_raccel_isMultitasking = 1 ;
#else
boolean_T gbl_raccel_isMultitasking = 0 ;
#endif
boolean_T gbl_raccel_tid01eq = 0 ; int_T gbl_raccel_NumST = 2 ; const char_T
* gbl_raccel_Version = "23.2 (R2023b) 01-Aug-2023" ; void
raccel_setup_MMIStateLog ( SimStruct * S ) {
#ifdef UseMMIDataLogging
rt_FillStateSigInfoFromMMI ( ssGetRTWLogInfo ( S ) , & ssGetErrorStatus ( S )
) ;
#else
UNUSED_PARAMETER ( S ) ;
#endif
} static DataMapInfo rt_dataMapInfo ; DataMapInfo * rt_dataMapInfoPtr = &
rt_dataMapInfo ; rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr = & (
rt_dataMapInfo . mmi ) ; int_T enableFcnCallFlag [ ] = { 1 , 1 } ; const char
* raccelLoadInputsAndAperiodicHitTimes ( SimStruct * S , const char *
inportFileName , int * matFileFormat ) { return rt_RAccelReadInportsMatFile (
S , inportFileName , matFileFormat ) ; }
#include "simstruc.h"
#include "fixedpoint.h"
#include "slsa_sim_engine.h"
#include "simtarget/slSimTgtSLExecSimBridge.h"
B rtB ; X rtX ; DW rtDW ; static SimStruct model_S ; SimStruct * const rtS =
& model_S ; void MdlStart ( void ) { CXPtMax * _rtXPerturbMax ; CXPtMin *
_rtXPerturbMin ; NeModelParameters modelParameters ; NeModelParameters
modelParameters_p ; NeslSimulationData * simulationData ; NeslSimulator * tmp
; NeuDiagnosticManager * diagnosticManager ; NeuDiagnosticTree *
diagnosticTree ; NeuDiagnosticTree * diagnosticTree_e ; NeuDiagnosticTree *
diagnosticTree_p ; char * msg ; char * msg_e ; char * msg_p ; real_T tmp_m [
168 ] ; real_T time ; real_T tmp_e ; int32_T tmp_i ; int_T tmp_g [ 43 ] ;
boolean_T tmp_p ; boolean_T val ; { bool externalInputIsInDatasetFormat =
false ; void * pISigstreamManager = rt_GetISigstreamManager ( rtS ) ;
rtwISigstreamManagerGetInputIsInDatasetFormat ( pISigstreamManager , &
externalInputIsInDatasetFormat ) ; if ( externalInputIsInDatasetFormat ) { }
} _rtXPerturbMax = ( ( CXPtMax * ) ssGetJacobianPerturbationBoundsMaxVec (
rtS ) ) ; _rtXPerturbMin = ( ( CXPtMin * )
ssGetJacobianPerturbationBoundsMinVec ( rtS ) ) ; { { { bool
isStreamoutAlreadyRegistered = false ; { sdiSignalSourceInfoU srcInfo ;
sdiLabelU loggedName = sdiGetLabelFromChars ( "Subsystem1" ) ; sdiLabelU
origSigName = sdiGetLabelFromChars ( "" ) ; sdiLabelU propName =
sdiGetLabelFromChars ( "Subsystem1" ) ; sdiLabelU blockPath =
sdiGetLabelFromChars ( "MagneticBasket_Simscape_Optimizer/To Workspace" ) ;
sdiLabelU blockSID = sdiGetLabelFromChars ( "" ) ; sdiLabelU subPath =
sdiGetLabelFromChars ( "" ) ; sdiDims sigDims ; sdiLabelU sigName =
sdiGetLabelFromChars ( "Subsystem1" ) ; sdiAsyncRepoDataTypeHandle hDT =
sdiAsyncRepoGetBuiltInDataTypeHandle ( DATA_TYPE_DOUBLE ) ; { sdiComplexity
sigComplexity = REAL ; sdiSampleTimeContinuity stCont =
SAMPLE_TIME_CONTINUOUS ; int_T sigDimsArray [ 1 ] = { 7 } ; sigDims . nDims =
1 ; sigDims . dimensions = sigDimsArray ; srcInfo . numBlockPathElems = 1 ;
srcInfo . fullBlockPath = ( sdiFullBlkPathU ) & blockPath ; srcInfo . SID = (
sdiSignalIDU ) & blockSID ; srcInfo . subPath = subPath ; srcInfo . portIndex
= 0 + 1 ; srcInfo . signalName = sigName ; srcInfo . sigSourceUUID = 0 ; rtDW
. l4r0r31ubh . AQHandles = sdiStartAsyncioQueueCreation ( hDT , & srcInfo ,
rt_dataMapInfo . mmi . InstanceMap . fullPath ,
"30d6c7e4-0241-434a-89fe-82fb47456870" , sigComplexity , & sigDims ,
DIMENSIONS_MODE_FIXED , stCont , "" ) ; sdiCompleteAsyncioQueueCreation (
rtDW . l4r0r31ubh . AQHandles , hDT , & srcInfo ) ; if ( rtDW . l4r0r31ubh .
AQHandles ) { sdiSetSignalSampleTimeString ( rtDW . l4r0r31ubh . AQHandles ,
"&#xC5F0;&#xC18D;" , 0.0 , ssGetTFinal ( rtS ) ) ; sdiSetSignalRefRate ( rtDW
. l4r0r31ubh . AQHandles , 0.0 ) ; sdiSetRunStartTime ( rtDW . l4r0r31ubh .
AQHandles , ssGetTaskTime ( rtS , 0 ) ) ; sdiAsyncRepoSetSignalExportSettings
( rtDW . l4r0r31ubh . AQHandles , 1 , 0 ) ; sdiAsyncRepoSetSignalExportName (
rtDW . l4r0r31ubh . AQHandles , loggedName , origSigName , propName ) ;
sdiAsyncRepoSetBlockPathDomain ( rtDW . l4r0r31ubh . AQHandles ) ; }
sdiFreeLabel ( sigName ) ; sdiFreeLabel ( loggedName ) ; sdiFreeLabel (
origSigName ) ; sdiFreeLabel ( propName ) ; sdiFreeLabel ( blockPath ) ;
sdiFreeLabel ( blockSID ) ; sdiFreeLabel ( subPath ) ; } } if ( !
isStreamoutAlreadyRegistered ) { { sdiLabelU varName = sdiGetLabelFromChars (
"Config" ) ; sdiRegisterWksVariable ( rtDW . l4r0r31ubh . AQHandles , varName
, "timeseries" ) ; sdiFreeLabel ( varName ) ; } } } } } tmp =
nesl_lease_simulator (
"MagneticBasket_Simscape_Optimizer/Solver Configuration_1" , 0 , 0 ) ; rtDW .
ctycytv1sy = ( void * ) tmp ; tmp_p = pointer_is_null ( rtDW . ctycytv1sy ) ;
if ( tmp_p ) { MagneticBasket_Simscape_Optimizer_dda62cd9_1_gateway ( ) ; tmp
= nesl_lease_simulator (
"MagneticBasket_Simscape_Optimizer/Solver Configuration_1" , 0 , 0 ) ; rtDW .
ctycytv1sy = ( void * ) tmp ; } slsaSaveRawMemoryForSimTargetOP ( rtS ,
"MagneticBasket_Simscape_Optimizer/Solver Configuration_100" , ( void * * ) (
& rtDW . ctycytv1sy ) , 0U * sizeof ( real_T ) , nesl_save_simdata ,
nesl_restore_simdata ) ; simulationData = nesl_create_simulation_data ( ) ;
rtDW . m2wgpmvjvy = ( void * ) simulationData ; diagnosticManager =
rtw_create_diagnostics ( ) ; rtDW . oqiguspxfn = ( void * ) diagnosticManager
; modelParameters . mSolverType = NE_SOLVER_TYPE_DAE ; modelParameters .
mSolverAbsTol = 0.001 ; modelParameters . mSolverRelTol = 0.001 ;
modelParameters . mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_MAYBE ;
modelParameters . mStartTime = 0.0 ; modelParameters . mLoadInitialState =
false ; modelParameters . mUseSimState = false ; modelParameters .
mLinTrimCompile = false ; modelParameters . mLoggingMode = SSC_LOGGING_OFF ;
modelParameters . mRTWModifiedTimeStamp = 6.61247773E+8 ; modelParameters .
mZcDisabled = false ; modelParameters . mUseModelRefSolver = false ;
modelParameters . mTargetFPGAHIL = false ; tmp_e = 0.001 ; modelParameters .
mSolverTolerance = tmp_e ; tmp_e = 0.0 ; modelParameters . mFixedStepSize =
tmp_e ; tmp_p = true ; modelParameters . mVariableStepSolver = tmp_p ; tmp_p
= false ; modelParameters . mIsUsingODEN = tmp_p ; tmp_p =
slIsRapidAcceleratorSimulating ( ) ; val = ssGetGlobalInitialStatesAvailable
( rtS ) ; if ( tmp_p ) { val = ( val && ssIsFirstInitCond ( rtS ) ) ; }
modelParameters . mLoadInitialState = val ; modelParameters . mZcDisabled =
false ; diagnosticManager = ( NeuDiagnosticManager * ) rtDW . oqiguspxfn ;
diagnosticTree = neu_diagnostic_manager_get_initial_tree ( diagnosticManager
) ; tmp_i = nesl_initialize_simulator ( ( NeslSimulator * ) rtDW . ctycytv1sy
, & modelParameters , diagnosticManager ) ; if ( tmp_i != 0 ) { tmp_p =
error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if ( tmp_p ) { msg =
rtw_diagnostics_msg ( diagnosticTree ) ; ssSetErrorStatus ( rtS , msg ) ; } }
simulationData = ( NeslSimulationData * ) rtDW . m2wgpmvjvy ; time = ssGetT (
rtS ) ; simulationData -> mData -> mTime . mN = 1 ; simulationData -> mData
-> mTime . mX = & time ; simulationData -> mData -> mContStates . mN = 18 ;
simulationData -> mData -> mContStates . mX = & rtX . eouclruk0d [ 0 ] ;
simulationData -> mData -> mDiscStates . mN = 0 ; simulationData -> mData ->
mDiscStates . mX = & rtDW . exk4rrasy3 ; simulationData -> mData ->
mModeVector . mN = 0 ; simulationData -> mData -> mModeVector . mX = & rtDW .
h2bg1hotvh ; tmp_p = ( ssIsMajorTimeStep ( rtS ) && ssGetRTWSolverInfo ( rtS
) -> foundContZcEvents ) ; simulationData -> mData -> mFoundZcEvents = tmp_p
; simulationData -> mData -> mIsMajorTimeStep = ssIsMajorTimeStep ( rtS ) ;
tmp_p = ( ssGetMdlInfoPtr ( rtS ) -> mdlFlags . solverAssertCheck == 1U ) ;
simulationData -> mData -> mIsSolverAssertCheck = tmp_p ; tmp_p =
ssIsSolverCheckingCIC ( rtS ) ; simulationData -> mData ->
mIsSolverCheckingCIC = tmp_p ; tmp_p = ssIsSolverComputingJacobian ( rtS ) ;
simulationData -> mData -> mIsComputingJacobian = tmp_p ; simulationData ->
mData -> mIsEvaluatingF0 = ( ssGetEvaluatingF0ForJacobian ( rtS ) != 0 ) ;
tmp_p = ssIsSolverRequestingReset ( rtS ) ; simulationData -> mData ->
mIsSolverRequestingReset = tmp_p ; simulationData -> mData ->
mIsModeUpdateTimeStep = ssIsModeUpdateTimeStep ( rtS ) ; tmp_g [ 0 ] = 0 ;
tmp_m [ 0 ] = rtB . c513o5pg5r [ 0 ] ; tmp_m [ 1 ] = rtB . c513o5pg5r [ 1 ] ;
tmp_m [ 2 ] = rtB . c513o5pg5r [ 2 ] ; tmp_m [ 3 ] = rtB . c513o5pg5r [ 3 ] ;
tmp_g [ 1 ] = 4 ; tmp_m [ 4 ] = rtB . pfmi5qvzhf [ 0 ] ; tmp_m [ 5 ] = rtB .
pfmi5qvzhf [ 1 ] ; tmp_m [ 6 ] = rtB . pfmi5qvzhf [ 2 ] ; tmp_m [ 7 ] = rtB .
pfmi5qvzhf [ 3 ] ; tmp_g [ 2 ] = 8 ; tmp_m [ 8 ] = rtB . gi35bo1nw4 [ 0 ] ;
tmp_m [ 9 ] = rtB . gi35bo1nw4 [ 1 ] ; tmp_m [ 10 ] = rtB . gi35bo1nw4 [ 2 ]
; tmp_m [ 11 ] = rtB . gi35bo1nw4 [ 3 ] ; tmp_g [ 3 ] = 12 ; tmp_m [ 12 ] =
rtB . iii5xyl5wn [ 0 ] ; tmp_m [ 13 ] = rtB . iii5xyl5wn [ 1 ] ; tmp_m [ 14 ]
= rtB . iii5xyl5wn [ 2 ] ; tmp_m [ 15 ] = rtB . iii5xyl5wn [ 3 ] ; tmp_g [ 4
] = 16 ; tmp_m [ 16 ] = rtB . mjqwg1vruy [ 0 ] ; tmp_m [ 17 ] = rtB .
mjqwg1vruy [ 1 ] ; tmp_m [ 18 ] = rtB . mjqwg1vruy [ 2 ] ; tmp_m [ 19 ] = rtB
. mjqwg1vruy [ 3 ] ; tmp_g [ 5 ] = 20 ; tmp_m [ 20 ] = rtB . eo0yhyzcnn [ 0 ]
; tmp_m [ 21 ] = rtB . eo0yhyzcnn [ 1 ] ; tmp_m [ 22 ] = rtB . eo0yhyzcnn [ 2
] ; tmp_m [ 23 ] = rtB . eo0yhyzcnn [ 3 ] ; tmp_g [ 6 ] = 24 ; tmp_m [ 24 ] =
rtB . c5dvbg4axm [ 0 ] ; tmp_m [ 25 ] = rtB . c5dvbg4axm [ 1 ] ; tmp_m [ 26 ]
= rtB . c5dvbg4axm [ 2 ] ; tmp_m [ 27 ] = rtB . c5dvbg4axm [ 3 ] ; tmp_g [ 7
] = 28 ; tmp_m [ 28 ] = rtB . bxbu5xamaf [ 0 ] ; tmp_m [ 29 ] = rtB .
bxbu5xamaf [ 1 ] ; tmp_m [ 30 ] = rtB . bxbu5xamaf [ 2 ] ; tmp_m [ 31 ] = rtB
. bxbu5xamaf [ 3 ] ; tmp_g [ 8 ] = 32 ; tmp_m [ 32 ] = rtB . foqtyz1zpt [ 0 ]
; tmp_m [ 33 ] = rtB . foqtyz1zpt [ 1 ] ; tmp_m [ 34 ] = rtB . foqtyz1zpt [ 2
] ; tmp_m [ 35 ] = rtB . foqtyz1zpt [ 3 ] ; tmp_g [ 9 ] = 36 ; tmp_m [ 36 ] =
rtB . gy3u0julyl [ 0 ] ; tmp_m [ 37 ] = rtB . gy3u0julyl [ 1 ] ; tmp_m [ 38 ]
= rtB . gy3u0julyl [ 2 ] ; tmp_m [ 39 ] = rtB . gy3u0julyl [ 3 ] ; tmp_g [ 10
] = 40 ; tmp_m [ 40 ] = rtB . hytxiurw0s [ 0 ] ; tmp_m [ 41 ] = rtB .
hytxiurw0s [ 1 ] ; tmp_m [ 42 ] = rtB . hytxiurw0s [ 2 ] ; tmp_m [ 43 ] = rtB
. hytxiurw0s [ 3 ] ; tmp_g [ 11 ] = 44 ; tmp_m [ 44 ] = rtB . nubcraskrw [ 0
] ; tmp_m [ 45 ] = rtB . nubcraskrw [ 1 ] ; tmp_m [ 46 ] = rtB . nubcraskrw [
2 ] ; tmp_m [ 47 ] = rtB . nubcraskrw [ 3 ] ; tmp_g [ 12 ] = 48 ; tmp_m [ 48
] = rtB . lp0cxari03 [ 0 ] ; tmp_m [ 49 ] = rtB . lp0cxari03 [ 1 ] ; tmp_m [
50 ] = rtB . lp0cxari03 [ 2 ] ; tmp_m [ 51 ] = rtB . lp0cxari03 [ 3 ] ; tmp_g
[ 13 ] = 52 ; tmp_m [ 52 ] = rtB . phx5lsta5w [ 0 ] ; tmp_m [ 53 ] = rtB .
phx5lsta5w [ 1 ] ; tmp_m [ 54 ] = rtB . phx5lsta5w [ 2 ] ; tmp_m [ 55 ] = rtB
. phx5lsta5w [ 3 ] ; tmp_g [ 14 ] = 56 ; tmp_m [ 56 ] = rtB . hoap0zr3eu [ 0
] ; tmp_m [ 57 ] = rtB . hoap0zr3eu [ 1 ] ; tmp_m [ 58 ] = rtB . hoap0zr3eu [
2 ] ; tmp_m [ 59 ] = rtB . hoap0zr3eu [ 3 ] ; tmp_g [ 15 ] = 60 ; tmp_m [ 60
] = rtB . enp5f3s002 [ 0 ] ; tmp_m [ 61 ] = rtB . enp5f3s002 [ 1 ] ; tmp_m [
62 ] = rtB . enp5f3s002 [ 2 ] ; tmp_m [ 63 ] = rtB . enp5f3s002 [ 3 ] ; tmp_g
[ 16 ] = 64 ; tmp_m [ 64 ] = rtB . nsl31ekm0f [ 0 ] ; tmp_m [ 65 ] = rtB .
nsl31ekm0f [ 1 ] ; tmp_m [ 66 ] = rtB . nsl31ekm0f [ 2 ] ; tmp_m [ 67 ] = rtB
. nsl31ekm0f [ 3 ] ; tmp_g [ 17 ] = 68 ; tmp_m [ 68 ] = rtB . hrqvwki2jl [ 0
] ; tmp_m [ 69 ] = rtB . hrqvwki2jl [ 1 ] ; tmp_m [ 70 ] = rtB . hrqvwki2jl [
2 ] ; tmp_m [ 71 ] = rtB . hrqvwki2jl [ 3 ] ; tmp_g [ 18 ] = 72 ; tmp_m [ 72
] = rtB . dypcby02yn [ 0 ] ; tmp_m [ 73 ] = rtB . dypcby02yn [ 1 ] ; tmp_m [
74 ] = rtB . dypcby02yn [ 2 ] ; tmp_m [ 75 ] = rtB . dypcby02yn [ 3 ] ; tmp_g
[ 19 ] = 76 ; tmp_m [ 76 ] = rtB . a100glwqkz [ 0 ] ; tmp_m [ 77 ] = rtB .
a100glwqkz [ 1 ] ; tmp_m [ 78 ] = rtB . a100glwqkz [ 2 ] ; tmp_m [ 79 ] = rtB
. a100glwqkz [ 3 ] ; tmp_g [ 20 ] = 80 ; tmp_m [ 80 ] = rtB . c3twv1nzdj [ 0
] ; tmp_m [ 81 ] = rtB . c3twv1nzdj [ 1 ] ; tmp_m [ 82 ] = rtB . c3twv1nzdj [
2 ] ; tmp_m [ 83 ] = rtB . c3twv1nzdj [ 3 ] ; tmp_g [ 21 ] = 84 ; tmp_m [ 84
] = rtB . iurh4h1cla [ 0 ] ; tmp_m [ 85 ] = rtB . iurh4h1cla [ 1 ] ; tmp_m [
86 ] = rtB . iurh4h1cla [ 2 ] ; tmp_m [ 87 ] = rtB . iurh4h1cla [ 3 ] ; tmp_g
[ 22 ] = 88 ; tmp_m [ 88 ] = rtB . har0a2wrmi [ 0 ] ; tmp_m [ 89 ] = rtB .
har0a2wrmi [ 1 ] ; tmp_m [ 90 ] = rtB . har0a2wrmi [ 2 ] ; tmp_m [ 91 ] = rtB
. har0a2wrmi [ 3 ] ; tmp_g [ 23 ] = 92 ; tmp_m [ 92 ] = rtB . ib2ieuf2an [ 0
] ; tmp_m [ 93 ] = rtB . ib2ieuf2an [ 1 ] ; tmp_m [ 94 ] = rtB . ib2ieuf2an [
2 ] ; tmp_m [ 95 ] = rtB . ib2ieuf2an [ 3 ] ; tmp_g [ 24 ] = 96 ; tmp_m [ 96
] = rtB . ojpovrodvr [ 0 ] ; tmp_m [ 97 ] = rtB . ojpovrodvr [ 1 ] ; tmp_m [
98 ] = rtB . ojpovrodvr [ 2 ] ; tmp_m [ 99 ] = rtB . ojpovrodvr [ 3 ] ; tmp_g
[ 25 ] = 100 ; tmp_m [ 100 ] = rtB . fcgmpow5uz [ 0 ] ; tmp_m [ 101 ] = rtB .
fcgmpow5uz [ 1 ] ; tmp_m [ 102 ] = rtB . fcgmpow5uz [ 2 ] ; tmp_m [ 103 ] =
rtB . fcgmpow5uz [ 3 ] ; tmp_g [ 26 ] = 104 ; tmp_m [ 104 ] = rtB .
jylntbxfpj [ 0 ] ; tmp_m [ 105 ] = rtB . jylntbxfpj [ 1 ] ; tmp_m [ 106 ] =
rtB . jylntbxfpj [ 2 ] ; tmp_m [ 107 ] = rtB . jylntbxfpj [ 3 ] ; tmp_g [ 27
] = 108 ; tmp_m [ 108 ] = rtB . g1yo2p3pxt [ 0 ] ; tmp_m [ 109 ] = rtB .
g1yo2p3pxt [ 1 ] ; tmp_m [ 110 ] = rtB . g1yo2p3pxt [ 2 ] ; tmp_m [ 111 ] =
rtB . g1yo2p3pxt [ 3 ] ; tmp_g [ 28 ] = 112 ; tmp_m [ 112 ] = rtB .
az4eo3zh5o [ 0 ] ; tmp_m [ 113 ] = rtB . az4eo3zh5o [ 1 ] ; tmp_m [ 114 ] =
rtB . az4eo3zh5o [ 2 ] ; tmp_m [ 115 ] = rtB . az4eo3zh5o [ 3 ] ; tmp_g [ 29
] = 116 ; tmp_m [ 116 ] = rtB . ch1tbnhyon [ 0 ] ; tmp_m [ 117 ] = rtB .
ch1tbnhyon [ 1 ] ; tmp_m [ 118 ] = rtB . ch1tbnhyon [ 2 ] ; tmp_m [ 119 ] =
rtB . ch1tbnhyon [ 3 ] ; tmp_g [ 30 ] = 120 ; tmp_m [ 120 ] = rtB .
eoseoj2f3b [ 0 ] ; tmp_m [ 121 ] = rtB . eoseoj2f3b [ 1 ] ; tmp_m [ 122 ] =
rtB . eoseoj2f3b [ 2 ] ; tmp_m [ 123 ] = rtB . eoseoj2f3b [ 3 ] ; tmp_g [ 31
] = 124 ; tmp_m [ 124 ] = rtB . aigpe1li2p [ 0 ] ; tmp_m [ 125 ] = rtB .
aigpe1li2p [ 1 ] ; tmp_m [ 126 ] = rtB . aigpe1li2p [ 2 ] ; tmp_m [ 127 ] =
rtB . aigpe1li2p [ 3 ] ; tmp_g [ 32 ] = 128 ; tmp_m [ 128 ] = rtB .
pbn2nzk2k3 [ 0 ] ; tmp_m [ 129 ] = rtB . pbn2nzk2k3 [ 1 ] ; tmp_m [ 130 ] =
rtB . pbn2nzk2k3 [ 2 ] ; tmp_m [ 131 ] = rtB . pbn2nzk2k3 [ 3 ] ; tmp_g [ 33
] = 132 ; tmp_m [ 132 ] = rtB . hffhwq1imf [ 0 ] ; tmp_m [ 133 ] = rtB .
hffhwq1imf [ 1 ] ; tmp_m [ 134 ] = rtB . hffhwq1imf [ 2 ] ; tmp_m [ 135 ] =
rtB . hffhwq1imf [ 3 ] ; tmp_g [ 34 ] = 136 ; tmp_m [ 136 ] = rtB .
gbxuu4wtte [ 0 ] ; tmp_m [ 137 ] = rtB . gbxuu4wtte [ 1 ] ; tmp_m [ 138 ] =
rtB . gbxuu4wtte [ 2 ] ; tmp_m [ 139 ] = rtB . gbxuu4wtte [ 3 ] ; tmp_g [ 35
] = 140 ; tmp_m [ 140 ] = rtB . k5nlme3kns [ 0 ] ; tmp_m [ 141 ] = rtB .
k5nlme3kns [ 1 ] ; tmp_m [ 142 ] = rtB . k5nlme3kns [ 2 ] ; tmp_m [ 143 ] =
rtB . k5nlme3kns [ 3 ] ; tmp_g [ 36 ] = 144 ; tmp_m [ 144 ] = rtB .
j52cfj2d1w [ 0 ] ; tmp_m [ 145 ] = rtB . j52cfj2d1w [ 1 ] ; tmp_m [ 146 ] =
rtB . j52cfj2d1w [ 2 ] ; tmp_m [ 147 ] = rtB . j52cfj2d1w [ 3 ] ; tmp_g [ 37
] = 148 ; tmp_m [ 148 ] = rtB . ftvz2ekipg [ 0 ] ; tmp_m [ 149 ] = rtB .
ftvz2ekipg [ 1 ] ; tmp_m [ 150 ] = rtB . ftvz2ekipg [ 2 ] ; tmp_m [ 151 ] =
rtB . ftvz2ekipg [ 3 ] ; tmp_g [ 38 ] = 152 ; tmp_m [ 152 ] = rtB .
awgskjjmsb [ 0 ] ; tmp_m [ 153 ] = rtB . awgskjjmsb [ 1 ] ; tmp_m [ 154 ] =
rtB . awgskjjmsb [ 2 ] ; tmp_m [ 155 ] = rtB . awgskjjmsb [ 3 ] ; tmp_g [ 39
] = 156 ; tmp_m [ 156 ] = rtB . fl4rvtkarn [ 0 ] ; tmp_m [ 157 ] = rtB .
fl4rvtkarn [ 1 ] ; tmp_m [ 158 ] = rtB . fl4rvtkarn [ 2 ] ; tmp_m [ 159 ] =
rtB . fl4rvtkarn [ 3 ] ; tmp_g [ 40 ] = 160 ; tmp_m [ 160 ] = rtB .
jp22vgjwwm [ 0 ] ; tmp_m [ 161 ] = rtB . jp22vgjwwm [ 1 ] ; tmp_m [ 162 ] =
rtB . jp22vgjwwm [ 2 ] ; tmp_m [ 163 ] = rtB . jp22vgjwwm [ 3 ] ; tmp_g [ 41
] = 164 ; tmp_m [ 164 ] = rtB . pzpjv5s5n2 [ 0 ] ; tmp_m [ 165 ] = rtB .
pzpjv5s5n2 [ 1 ] ; tmp_m [ 166 ] = rtB . pzpjv5s5n2 [ 2 ] ; tmp_m [ 167 ] =
rtB . pzpjv5s5n2 [ 3 ] ; tmp_g [ 42 ] = 168 ; simulationData -> mData ->
mInputValues . mN = 168 ; simulationData -> mData -> mInputValues . mX = &
tmp_m [ 0 ] ; simulationData -> mData -> mInputOffsets . mN = 43 ;
simulationData -> mData -> mInputOffsets . mX = & tmp_g [ 0 ] ;
simulationData -> mData -> mNumjacDxLo . mN = 18 ; simulationData -> mData ->
mNumjacDxLo . mX = & _rtXPerturbMin -> eouclruk0d [ 0 ] ; simulationData ->
mData -> mNumjacDxHi . mN = 18 ; simulationData -> mData -> mNumjacDxHi . mX
= & _rtXPerturbMax -> eouclruk0d [ 0 ] ; diagnosticManager = (
NeuDiagnosticManager * ) rtDW . oqiguspxfn ; diagnosticTree_p =
neu_diagnostic_manager_get_initial_tree ( diagnosticManager ) ; tmp_i =
ne_simulator_method ( ( NeslSimulator * ) rtDW . ctycytv1sy ,
NESL_SIM_NUMJAC_DX_BOUNDS , simulationData , diagnosticManager ) ; if ( tmp_i
!= 0 ) { tmp_p = error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if (
tmp_p ) { msg_p = rtw_diagnostics_msg ( diagnosticTree_p ) ; ssSetErrorStatus
( rtS , msg_p ) ; } } tmp = nesl_lease_simulator (
"MagneticBasket_Simscape_Optimizer/Solver Configuration_1" , 1 , 0 ) ; rtDW .
dikp21y4fz = ( void * ) tmp ; tmp_p = pointer_is_null ( rtDW . dikp21y4fz ) ;
if ( tmp_p ) { MagneticBasket_Simscape_Optimizer_dda62cd9_1_gateway ( ) ; tmp
= nesl_lease_simulator (
"MagneticBasket_Simscape_Optimizer/Solver Configuration_1" , 1 , 0 ) ; rtDW .
dikp21y4fz = ( void * ) tmp ; } slsaSaveRawMemoryForSimTargetOP ( rtS ,
"MagneticBasket_Simscape_Optimizer/Solver Configuration_110" , ( void * * ) (
& rtDW . dikp21y4fz ) , 0U * sizeof ( real_T ) , nesl_save_simdata ,
nesl_restore_simdata ) ; simulationData = nesl_create_simulation_data ( ) ;
rtDW . fqoysoyi51 = ( void * ) simulationData ; diagnosticManager =
rtw_create_diagnostics ( ) ; rtDW . dcyq5ws3qs = ( void * ) diagnosticManager
; modelParameters_p . mSolverType = NE_SOLVER_TYPE_DAE ; modelParameters_p .
mSolverAbsTol = 0.001 ; modelParameters_p . mSolverRelTol = 0.001 ;
modelParameters_p . mSolverModifyAbsTol = NE_MODIFY_ABS_TOL_MAYBE ;
modelParameters_p . mStartTime = 0.0 ; modelParameters_p . mLoadInitialState
= false ; modelParameters_p . mUseSimState = false ; modelParameters_p .
mLinTrimCompile = false ; modelParameters_p . mLoggingMode = SSC_LOGGING_OFF
; modelParameters_p . mRTWModifiedTimeStamp = 6.61247773E+8 ;
modelParameters_p . mZcDisabled = false ; modelParameters_p .
mUseModelRefSolver = false ; modelParameters_p . mTargetFPGAHIL = false ;
tmp_e = 0.001 ; modelParameters_p . mSolverTolerance = tmp_e ; tmp_e = 0.0 ;
modelParameters_p . mFixedStepSize = tmp_e ; tmp_p = true ; modelParameters_p
. mVariableStepSolver = tmp_p ; tmp_p = false ; modelParameters_p .
mIsUsingODEN = tmp_p ; tmp_p = slIsRapidAcceleratorSimulating ( ) ; val =
ssGetGlobalInitialStatesAvailable ( rtS ) ; if ( tmp_p ) { val = ( val &&
ssIsFirstInitCond ( rtS ) ) ; } modelParameters_p . mLoadInitialState = val ;
modelParameters_p . mZcDisabled = false ; diagnosticManager = (
NeuDiagnosticManager * ) rtDW . dcyq5ws3qs ; diagnosticTree_e =
neu_diagnostic_manager_get_initial_tree ( diagnosticManager ) ; tmp_i =
nesl_initialize_simulator ( ( NeslSimulator * ) rtDW . dikp21y4fz , &
modelParameters_p , diagnosticManager ) ; if ( tmp_i != 0 ) { tmp_p =
error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if ( tmp_p ) { msg_e =
rtw_diagnostics_msg ( diagnosticTree_e ) ; ssSetErrorStatus ( rtS , msg_e ) ;
} } rtDW . o3hz2jqijw = 0 ; rtDW . a1ctm05tly = 0 ; rtDW . eiqm0ndfyg = 0 ;
rtDW . pibyff0eql = 0 ; rtDW . l4r1nqhb1g = 0 ; rtDW . emuswyjmhs = 0 ; rtDW
. ogqlv3h0ke = 0 ; rtDW . hgiqtwygy1 = 0 ; rtDW . a1er4qtcuo = 0 ; rtDW .
jtymq0r3wm = 0 ; rtDW . bj4ejpvtiv = 0 ; rtDW . eyhr14ycpb = 0 ; rtDW .
gmlynrr0kl = 0 ; rtDW . dfembtpulp = 0 ; } void MdlOutputs ( int_T tid ) {
NeslSimulationData * simulationData ; NeuDiagnosticManager *
diagnosticManager ; NeuDiagnosticTree * diagnosticTree ; NeuDiagnosticTree *
diagnosticTree_p ; char * msg ; char * msg_p ; real_T tmp_i [ 186 ] ; real_T
tmp_p [ 168 ] ; real_T cazgyvwuzb [ 3 ] ; real_T czut5pi4ll ; real_T
jj1qbnwy23_idx_1 ; real_T jj1qbnwy23_idx_2 ; real_T jj1qbnwy23_idx_3 ; real_T
nsmqj1yv3h ; real_T time ; real_T time_e ; real_T time_i ; real_T time_p ;
real_T * tmp_g ; int32_T yIdx ; int_T tmp_m [ 44 ] ; int_T tmp_e [ 43 ] ;
boolean_T tmp ; simulationData = ( NeslSimulationData * ) rtDW . m2wgpmvjvy ;
time = ssGetT ( rtS ) ; simulationData -> mData -> mTime . mN = 1 ;
simulationData -> mData -> mTime . mX = & time ; simulationData -> mData ->
mContStates . mN = 18 ; simulationData -> mData -> mContStates . mX = & rtX .
eouclruk0d [ 0 ] ; simulationData -> mData -> mDiscStates . mN = 0 ;
simulationData -> mData -> mDiscStates . mX = & rtDW . exk4rrasy3 ;
simulationData -> mData -> mModeVector . mN = 0 ; simulationData -> mData ->
mModeVector . mX = & rtDW . h2bg1hotvh ; tmp = ( ssIsMajorTimeStep ( rtS ) &&
ssGetRTWSolverInfo ( rtS ) -> foundContZcEvents ) ; simulationData -> mData
-> mFoundZcEvents = tmp ; simulationData -> mData -> mIsMajorTimeStep =
ssIsMajorTimeStep ( rtS ) ; tmp = ( ssGetMdlInfoPtr ( rtS ) -> mdlFlags .
solverAssertCheck == 1U ) ; simulationData -> mData -> mIsSolverAssertCheck =
tmp ; tmp = ssIsSolverCheckingCIC ( rtS ) ; simulationData -> mData ->
mIsSolverCheckingCIC = tmp ; tmp = ssIsSolverComputingJacobian ( rtS ) ;
simulationData -> mData -> mIsComputingJacobian = tmp ; simulationData ->
mData -> mIsEvaluatingF0 = ( ssGetEvaluatingF0ForJacobian ( rtS ) != 0 ) ;
tmp = ssIsSolverRequestingReset ( rtS ) ; simulationData -> mData ->
mIsSolverRequestingReset = tmp ; simulationData -> mData ->
mIsModeUpdateTimeStep = ssIsModeUpdateTimeStep ( rtS ) ; tmp_e [ 0 ] = 0 ;
tmp_p [ 0 ] = rtB . c513o5pg5r [ 0 ] ; tmp_p [ 1 ] = rtB . c513o5pg5r [ 1 ] ;
tmp_p [ 2 ] = rtB . c513o5pg5r [ 2 ] ; tmp_p [ 3 ] = rtB . c513o5pg5r [ 3 ] ;
tmp_e [ 1 ] = 4 ; tmp_p [ 4 ] = rtB . pfmi5qvzhf [ 0 ] ; tmp_p [ 5 ] = rtB .
pfmi5qvzhf [ 1 ] ; tmp_p [ 6 ] = rtB . pfmi5qvzhf [ 2 ] ; tmp_p [ 7 ] = rtB .
pfmi5qvzhf [ 3 ] ; tmp_e [ 2 ] = 8 ; tmp_p [ 8 ] = rtB . gi35bo1nw4 [ 0 ] ;
tmp_p [ 9 ] = rtB . gi35bo1nw4 [ 1 ] ; tmp_p [ 10 ] = rtB . gi35bo1nw4 [ 2 ]
; tmp_p [ 11 ] = rtB . gi35bo1nw4 [ 3 ] ; tmp_e [ 3 ] = 12 ; tmp_p [ 12 ] =
rtB . iii5xyl5wn [ 0 ] ; tmp_p [ 13 ] = rtB . iii5xyl5wn [ 1 ] ; tmp_p [ 14 ]
= rtB . iii5xyl5wn [ 2 ] ; tmp_p [ 15 ] = rtB . iii5xyl5wn [ 3 ] ; tmp_e [ 4
] = 16 ; tmp_p [ 16 ] = rtB . mjqwg1vruy [ 0 ] ; tmp_p [ 17 ] = rtB .
mjqwg1vruy [ 1 ] ; tmp_p [ 18 ] = rtB . mjqwg1vruy [ 2 ] ; tmp_p [ 19 ] = rtB
. mjqwg1vruy [ 3 ] ; tmp_e [ 5 ] = 20 ; tmp_p [ 20 ] = rtB . eo0yhyzcnn [ 0 ]
; tmp_p [ 21 ] = rtB . eo0yhyzcnn [ 1 ] ; tmp_p [ 22 ] = rtB . eo0yhyzcnn [ 2
] ; tmp_p [ 23 ] = rtB . eo0yhyzcnn [ 3 ] ; tmp_e [ 6 ] = 24 ; tmp_p [ 24 ] =
rtB . c5dvbg4axm [ 0 ] ; tmp_p [ 25 ] = rtB . c5dvbg4axm [ 1 ] ; tmp_p [ 26 ]
= rtB . c5dvbg4axm [ 2 ] ; tmp_p [ 27 ] = rtB . c5dvbg4axm [ 3 ] ; tmp_e [ 7
] = 28 ; tmp_p [ 28 ] = rtB . bxbu5xamaf [ 0 ] ; tmp_p [ 29 ] = rtB .
bxbu5xamaf [ 1 ] ; tmp_p [ 30 ] = rtB . bxbu5xamaf [ 2 ] ; tmp_p [ 31 ] = rtB
. bxbu5xamaf [ 3 ] ; tmp_e [ 8 ] = 32 ; tmp_p [ 32 ] = rtB . foqtyz1zpt [ 0 ]
; tmp_p [ 33 ] = rtB . foqtyz1zpt [ 1 ] ; tmp_p [ 34 ] = rtB . foqtyz1zpt [ 2
] ; tmp_p [ 35 ] = rtB . foqtyz1zpt [ 3 ] ; tmp_e [ 9 ] = 36 ; tmp_p [ 36 ] =
rtB . gy3u0julyl [ 0 ] ; tmp_p [ 37 ] = rtB . gy3u0julyl [ 1 ] ; tmp_p [ 38 ]
= rtB . gy3u0julyl [ 2 ] ; tmp_p [ 39 ] = rtB . gy3u0julyl [ 3 ] ; tmp_e [ 10
] = 40 ; tmp_p [ 40 ] = rtB . hytxiurw0s [ 0 ] ; tmp_p [ 41 ] = rtB .
hytxiurw0s [ 1 ] ; tmp_p [ 42 ] = rtB . hytxiurw0s [ 2 ] ; tmp_p [ 43 ] = rtB
. hytxiurw0s [ 3 ] ; tmp_e [ 11 ] = 44 ; tmp_p [ 44 ] = rtB . nubcraskrw [ 0
] ; tmp_p [ 45 ] = rtB . nubcraskrw [ 1 ] ; tmp_p [ 46 ] = rtB . nubcraskrw [
2 ] ; tmp_p [ 47 ] = rtB . nubcraskrw [ 3 ] ; tmp_e [ 12 ] = 48 ; tmp_p [ 48
] = rtB . lp0cxari03 [ 0 ] ; tmp_p [ 49 ] = rtB . lp0cxari03 [ 1 ] ; tmp_p [
50 ] = rtB . lp0cxari03 [ 2 ] ; tmp_p [ 51 ] = rtB . lp0cxari03 [ 3 ] ; tmp_e
[ 13 ] = 52 ; tmp_p [ 52 ] = rtB . phx5lsta5w [ 0 ] ; tmp_p [ 53 ] = rtB .
phx5lsta5w [ 1 ] ; tmp_p [ 54 ] = rtB . phx5lsta5w [ 2 ] ; tmp_p [ 55 ] = rtB
. phx5lsta5w [ 3 ] ; tmp_e [ 14 ] = 56 ; tmp_p [ 56 ] = rtB . hoap0zr3eu [ 0
] ; tmp_p [ 57 ] = rtB . hoap0zr3eu [ 1 ] ; tmp_p [ 58 ] = rtB . hoap0zr3eu [
2 ] ; tmp_p [ 59 ] = rtB . hoap0zr3eu [ 3 ] ; tmp_e [ 15 ] = 60 ; tmp_p [ 60
] = rtB . enp5f3s002 [ 0 ] ; tmp_p [ 61 ] = rtB . enp5f3s002 [ 1 ] ; tmp_p [
62 ] = rtB . enp5f3s002 [ 2 ] ; tmp_p [ 63 ] = rtB . enp5f3s002 [ 3 ] ; tmp_e
[ 16 ] = 64 ; tmp_p [ 64 ] = rtB . nsl31ekm0f [ 0 ] ; tmp_p [ 65 ] = rtB .
nsl31ekm0f [ 1 ] ; tmp_p [ 66 ] = rtB . nsl31ekm0f [ 2 ] ; tmp_p [ 67 ] = rtB
. nsl31ekm0f [ 3 ] ; tmp_e [ 17 ] = 68 ; tmp_p [ 68 ] = rtB . hrqvwki2jl [ 0
] ; tmp_p [ 69 ] = rtB . hrqvwki2jl [ 1 ] ; tmp_p [ 70 ] = rtB . hrqvwki2jl [
2 ] ; tmp_p [ 71 ] = rtB . hrqvwki2jl [ 3 ] ; tmp_e [ 18 ] = 72 ; tmp_p [ 72
] = rtB . dypcby02yn [ 0 ] ; tmp_p [ 73 ] = rtB . dypcby02yn [ 1 ] ; tmp_p [
74 ] = rtB . dypcby02yn [ 2 ] ; tmp_p [ 75 ] = rtB . dypcby02yn [ 3 ] ; tmp_e
[ 19 ] = 76 ; tmp_p [ 76 ] = rtB . a100glwqkz [ 0 ] ; tmp_p [ 77 ] = rtB .
a100glwqkz [ 1 ] ; tmp_p [ 78 ] = rtB . a100glwqkz [ 2 ] ; tmp_p [ 79 ] = rtB
. a100glwqkz [ 3 ] ; tmp_e [ 20 ] = 80 ; tmp_p [ 80 ] = rtB . c3twv1nzdj [ 0
] ; tmp_p [ 81 ] = rtB . c3twv1nzdj [ 1 ] ; tmp_p [ 82 ] = rtB . c3twv1nzdj [
2 ] ; tmp_p [ 83 ] = rtB . c3twv1nzdj [ 3 ] ; tmp_e [ 21 ] = 84 ; tmp_p [ 84
] = rtB . iurh4h1cla [ 0 ] ; tmp_p [ 85 ] = rtB . iurh4h1cla [ 1 ] ; tmp_p [
86 ] = rtB . iurh4h1cla [ 2 ] ; tmp_p [ 87 ] = rtB . iurh4h1cla [ 3 ] ; tmp_e
[ 22 ] = 88 ; tmp_p [ 88 ] = rtB . har0a2wrmi [ 0 ] ; tmp_p [ 89 ] = rtB .
har0a2wrmi [ 1 ] ; tmp_p [ 90 ] = rtB . har0a2wrmi [ 2 ] ; tmp_p [ 91 ] = rtB
. har0a2wrmi [ 3 ] ; tmp_e [ 23 ] = 92 ; tmp_p [ 92 ] = rtB . ib2ieuf2an [ 0
] ; tmp_p [ 93 ] = rtB . ib2ieuf2an [ 1 ] ; tmp_p [ 94 ] = rtB . ib2ieuf2an [
2 ] ; tmp_p [ 95 ] = rtB . ib2ieuf2an [ 3 ] ; tmp_e [ 24 ] = 96 ; tmp_p [ 96
] = rtB . ojpovrodvr [ 0 ] ; tmp_p [ 97 ] = rtB . ojpovrodvr [ 1 ] ; tmp_p [
98 ] = rtB . ojpovrodvr [ 2 ] ; tmp_p [ 99 ] = rtB . ojpovrodvr [ 3 ] ; tmp_e
[ 25 ] = 100 ; tmp_p [ 100 ] = rtB . fcgmpow5uz [ 0 ] ; tmp_p [ 101 ] = rtB .
fcgmpow5uz [ 1 ] ; tmp_p [ 102 ] = rtB . fcgmpow5uz [ 2 ] ; tmp_p [ 103 ] =
rtB . fcgmpow5uz [ 3 ] ; tmp_e [ 26 ] = 104 ; tmp_p [ 104 ] = rtB .
jylntbxfpj [ 0 ] ; tmp_p [ 105 ] = rtB . jylntbxfpj [ 1 ] ; tmp_p [ 106 ] =
rtB . jylntbxfpj [ 2 ] ; tmp_p [ 107 ] = rtB . jylntbxfpj [ 3 ] ; tmp_e [ 27
] = 108 ; tmp_p [ 108 ] = rtB . g1yo2p3pxt [ 0 ] ; tmp_p [ 109 ] = rtB .
g1yo2p3pxt [ 1 ] ; tmp_p [ 110 ] = rtB . g1yo2p3pxt [ 2 ] ; tmp_p [ 111 ] =
rtB . g1yo2p3pxt [ 3 ] ; tmp_e [ 28 ] = 112 ; tmp_p [ 112 ] = rtB .
az4eo3zh5o [ 0 ] ; tmp_p [ 113 ] = rtB . az4eo3zh5o [ 1 ] ; tmp_p [ 114 ] =
rtB . az4eo3zh5o [ 2 ] ; tmp_p [ 115 ] = rtB . az4eo3zh5o [ 3 ] ; tmp_e [ 29
] = 116 ; tmp_p [ 116 ] = rtB . ch1tbnhyon [ 0 ] ; tmp_p [ 117 ] = rtB .
ch1tbnhyon [ 1 ] ; tmp_p [ 118 ] = rtB . ch1tbnhyon [ 2 ] ; tmp_p [ 119 ] =
rtB . ch1tbnhyon [ 3 ] ; tmp_e [ 30 ] = 120 ; tmp_p [ 120 ] = rtB .
eoseoj2f3b [ 0 ] ; tmp_p [ 121 ] = rtB . eoseoj2f3b [ 1 ] ; tmp_p [ 122 ] =
rtB . eoseoj2f3b [ 2 ] ; tmp_p [ 123 ] = rtB . eoseoj2f3b [ 3 ] ; tmp_e [ 31
] = 124 ; tmp_p [ 124 ] = rtB . aigpe1li2p [ 0 ] ; tmp_p [ 125 ] = rtB .
aigpe1li2p [ 1 ] ; tmp_p [ 126 ] = rtB . aigpe1li2p [ 2 ] ; tmp_p [ 127 ] =
rtB . aigpe1li2p [ 3 ] ; tmp_e [ 32 ] = 128 ; tmp_p [ 128 ] = rtB .
pbn2nzk2k3 [ 0 ] ; tmp_p [ 129 ] = rtB . pbn2nzk2k3 [ 1 ] ; tmp_p [ 130 ] =
rtB . pbn2nzk2k3 [ 2 ] ; tmp_p [ 131 ] = rtB . pbn2nzk2k3 [ 3 ] ; tmp_e [ 33
] = 132 ; tmp_p [ 132 ] = rtB . hffhwq1imf [ 0 ] ; tmp_p [ 133 ] = rtB .
hffhwq1imf [ 1 ] ; tmp_p [ 134 ] = rtB . hffhwq1imf [ 2 ] ; tmp_p [ 135 ] =
rtB . hffhwq1imf [ 3 ] ; tmp_e [ 34 ] = 136 ; tmp_p [ 136 ] = rtB .
gbxuu4wtte [ 0 ] ; tmp_p [ 137 ] = rtB . gbxuu4wtte [ 1 ] ; tmp_p [ 138 ] =
rtB . gbxuu4wtte [ 2 ] ; tmp_p [ 139 ] = rtB . gbxuu4wtte [ 3 ] ; tmp_e [ 35
] = 140 ; tmp_p [ 140 ] = rtB . k5nlme3kns [ 0 ] ; tmp_p [ 141 ] = rtB .
k5nlme3kns [ 1 ] ; tmp_p [ 142 ] = rtB . k5nlme3kns [ 2 ] ; tmp_p [ 143 ] =
rtB . k5nlme3kns [ 3 ] ; tmp_e [ 36 ] = 144 ; tmp_p [ 144 ] = rtB .
j52cfj2d1w [ 0 ] ; tmp_p [ 145 ] = rtB . j52cfj2d1w [ 1 ] ; tmp_p [ 146 ] =
rtB . j52cfj2d1w [ 2 ] ; tmp_p [ 147 ] = rtB . j52cfj2d1w [ 3 ] ; tmp_e [ 37
] = 148 ; tmp_p [ 148 ] = rtB . ftvz2ekipg [ 0 ] ; tmp_p [ 149 ] = rtB .
ftvz2ekipg [ 1 ] ; tmp_p [ 150 ] = rtB . ftvz2ekipg [ 2 ] ; tmp_p [ 151 ] =
rtB . ftvz2ekipg [ 3 ] ; tmp_e [ 38 ] = 152 ; tmp_p [ 152 ] = rtB .
awgskjjmsb [ 0 ] ; tmp_p [ 153 ] = rtB . awgskjjmsb [ 1 ] ; tmp_p [ 154 ] =
rtB . awgskjjmsb [ 2 ] ; tmp_p [ 155 ] = rtB . awgskjjmsb [ 3 ] ; tmp_e [ 39
] = 156 ; tmp_p [ 156 ] = rtB . fl4rvtkarn [ 0 ] ; tmp_p [ 157 ] = rtB .
fl4rvtkarn [ 1 ] ; tmp_p [ 158 ] = rtB . fl4rvtkarn [ 2 ] ; tmp_p [ 159 ] =
rtB . fl4rvtkarn [ 3 ] ; tmp_e [ 40 ] = 160 ; tmp_p [ 160 ] = rtB .
jp22vgjwwm [ 0 ] ; tmp_p [ 161 ] = rtB . jp22vgjwwm [ 1 ] ; tmp_p [ 162 ] =
rtB . jp22vgjwwm [ 2 ] ; tmp_p [ 163 ] = rtB . jp22vgjwwm [ 3 ] ; tmp_e [ 41
] = 164 ; tmp_p [ 164 ] = rtB . pzpjv5s5n2 [ 0 ] ; tmp_p [ 165 ] = rtB .
pzpjv5s5n2 [ 1 ] ; tmp_p [ 166 ] = rtB . pzpjv5s5n2 [ 2 ] ; tmp_p [ 167 ] =
rtB . pzpjv5s5n2 [ 3 ] ; tmp_e [ 42 ] = 168 ; simulationData -> mData ->
mInputValues . mN = 168 ; simulationData -> mData -> mInputValues . mX = &
tmp_p [ 0 ] ; simulationData -> mData -> mInputOffsets . mN = 43 ;
simulationData -> mData -> mInputOffsets . mX = & tmp_e [ 0 ] ;
simulationData -> mData -> mOutputs . mN = 18 ; simulationData -> mData ->
mOutputs . mX = & rtB . aay54cmgia [ 0 ] ; simulationData -> mData ->
mTolerances . mN = 0 ; simulationData -> mData -> mTolerances . mX = NULL ;
simulationData -> mData -> mCstateHasChanged = false ; time_p = ssGetTaskTime
( rtS , 0 ) ; simulationData -> mData -> mTime . mN = 1 ; simulationData ->
mData -> mTime . mX = & time_p ; simulationData -> mData -> mSampleHits . mN
= 0 ; simulationData -> mData -> mSampleHits . mX = NULL ; simulationData ->
mData -> mIsFundamentalSampleHit = false ; diagnosticManager = (
NeuDiagnosticManager * ) rtDW . oqiguspxfn ; diagnosticTree =
neu_diagnostic_manager_get_initial_tree ( diagnosticManager ) ; yIdx =
ne_simulator_method ( ( NeslSimulator * ) rtDW . ctycytv1sy ,
NESL_SIM_OUTPUTS , simulationData , diagnosticManager ) ; if ( yIdx != 0 ) {
tmp = error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if ( tmp ) { msg =
rtw_diagnostics_msg ( diagnosticTree ) ; ssSetErrorStatus ( rtS , msg ) ; } }
if ( ssIsMajorTimeStep ( rtS ) && simulationData -> mData ->
mCstateHasChanged ) { ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ; }
simulationData = ( NeslSimulationData * ) rtDW . fqoysoyi51 ; time_e = ssGetT
( rtS ) ; simulationData -> mData -> mTime . mN = 1 ; simulationData -> mData
-> mTime . mX = & time_e ; simulationData -> mData -> mContStates . mN = 0 ;
simulationData -> mData -> mContStates . mX = NULL ; simulationData -> mData
-> mDiscStates . mN = 0 ; simulationData -> mData -> mDiscStates . mX = &
rtDW . oz0f0tgc3u ; simulationData -> mData -> mModeVector . mN = 0 ;
simulationData -> mData -> mModeVector . mX = & rtDW . glrxfbh4lh ; tmp = (
ssIsMajorTimeStep ( rtS ) && ssGetRTWSolverInfo ( rtS ) -> foundContZcEvents
) ; simulationData -> mData -> mFoundZcEvents = tmp ; simulationData -> mData
-> mIsMajorTimeStep = ssIsMajorTimeStep ( rtS ) ; tmp = ( ssGetMdlInfoPtr (
rtS ) -> mdlFlags . solverAssertCheck == 1U ) ; simulationData -> mData ->
mIsSolverAssertCheck = tmp ; tmp = ssIsSolverCheckingCIC ( rtS ) ;
simulationData -> mData -> mIsSolverCheckingCIC = tmp ; simulationData ->
mData -> mIsComputingJacobian = false ; simulationData -> mData ->
mIsEvaluatingF0 = false ; tmp = ssIsSolverRequestingReset ( rtS ) ;
simulationData -> mData -> mIsSolverRequestingReset = tmp ; simulationData ->
mData -> mIsModeUpdateTimeStep = ssIsModeUpdateTimeStep ( rtS ) ; tmp_m [ 0 ]
= 0 ; tmp_i [ 0 ] = rtB . c513o5pg5r [ 0 ] ; tmp_i [ 1 ] = rtB . c513o5pg5r [
1 ] ; tmp_i [ 2 ] = rtB . c513o5pg5r [ 2 ] ; tmp_i [ 3 ] = rtB . c513o5pg5r [
3 ] ; tmp_m [ 1 ] = 4 ; tmp_i [ 4 ] = rtB . pfmi5qvzhf [ 0 ] ; tmp_i [ 5 ] =
rtB . pfmi5qvzhf [ 1 ] ; tmp_i [ 6 ] = rtB . pfmi5qvzhf [ 2 ] ; tmp_i [ 7 ] =
rtB . pfmi5qvzhf [ 3 ] ; tmp_m [ 2 ] = 8 ; tmp_i [ 8 ] = rtB . gi35bo1nw4 [ 0
] ; tmp_i [ 9 ] = rtB . gi35bo1nw4 [ 1 ] ; tmp_i [ 10 ] = rtB . gi35bo1nw4 [
2 ] ; tmp_i [ 11 ] = rtB . gi35bo1nw4 [ 3 ] ; tmp_m [ 3 ] = 12 ; tmp_i [ 12 ]
= rtB . iii5xyl5wn [ 0 ] ; tmp_i [ 13 ] = rtB . iii5xyl5wn [ 1 ] ; tmp_i [ 14
] = rtB . iii5xyl5wn [ 2 ] ; tmp_i [ 15 ] = rtB . iii5xyl5wn [ 3 ] ; tmp_m [
4 ] = 16 ; tmp_i [ 16 ] = rtB . mjqwg1vruy [ 0 ] ; tmp_i [ 17 ] = rtB .
mjqwg1vruy [ 1 ] ; tmp_i [ 18 ] = rtB . mjqwg1vruy [ 2 ] ; tmp_i [ 19 ] = rtB
. mjqwg1vruy [ 3 ] ; tmp_m [ 5 ] = 20 ; tmp_i [ 20 ] = rtB . eo0yhyzcnn [ 0 ]
; tmp_i [ 21 ] = rtB . eo0yhyzcnn [ 1 ] ; tmp_i [ 22 ] = rtB . eo0yhyzcnn [ 2
] ; tmp_i [ 23 ] = rtB . eo0yhyzcnn [ 3 ] ; tmp_m [ 6 ] = 24 ; tmp_i [ 24 ] =
rtB . c5dvbg4axm [ 0 ] ; tmp_i [ 25 ] = rtB . c5dvbg4axm [ 1 ] ; tmp_i [ 26 ]
= rtB . c5dvbg4axm [ 2 ] ; tmp_i [ 27 ] = rtB . c5dvbg4axm [ 3 ] ; tmp_m [ 7
] = 28 ; tmp_i [ 28 ] = rtB . bxbu5xamaf [ 0 ] ; tmp_i [ 29 ] = rtB .
bxbu5xamaf [ 1 ] ; tmp_i [ 30 ] = rtB . bxbu5xamaf [ 2 ] ; tmp_i [ 31 ] = rtB
. bxbu5xamaf [ 3 ] ; tmp_m [ 8 ] = 32 ; tmp_i [ 32 ] = rtB . foqtyz1zpt [ 0 ]
; tmp_i [ 33 ] = rtB . foqtyz1zpt [ 1 ] ; tmp_i [ 34 ] = rtB . foqtyz1zpt [ 2
] ; tmp_i [ 35 ] = rtB . foqtyz1zpt [ 3 ] ; tmp_m [ 9 ] = 36 ; tmp_i [ 36 ] =
rtB . gy3u0julyl [ 0 ] ; tmp_i [ 37 ] = rtB . gy3u0julyl [ 1 ] ; tmp_i [ 38 ]
= rtB . gy3u0julyl [ 2 ] ; tmp_i [ 39 ] = rtB . gy3u0julyl [ 3 ] ; tmp_m [ 10
] = 40 ; tmp_i [ 40 ] = rtB . hytxiurw0s [ 0 ] ; tmp_i [ 41 ] = rtB .
hytxiurw0s [ 1 ] ; tmp_i [ 42 ] = rtB . hytxiurw0s [ 2 ] ; tmp_i [ 43 ] = rtB
. hytxiurw0s [ 3 ] ; tmp_m [ 11 ] = 44 ; tmp_i [ 44 ] = rtB . nubcraskrw [ 0
] ; tmp_i [ 45 ] = rtB . nubcraskrw [ 1 ] ; tmp_i [ 46 ] = rtB . nubcraskrw [
2 ] ; tmp_i [ 47 ] = rtB . nubcraskrw [ 3 ] ; tmp_m [ 12 ] = 48 ; tmp_i [ 48
] = rtB . lp0cxari03 [ 0 ] ; tmp_i [ 49 ] = rtB . lp0cxari03 [ 1 ] ; tmp_i [
50 ] = rtB . lp0cxari03 [ 2 ] ; tmp_i [ 51 ] = rtB . lp0cxari03 [ 3 ] ; tmp_m
[ 13 ] = 52 ; tmp_i [ 52 ] = rtB . phx5lsta5w [ 0 ] ; tmp_i [ 53 ] = rtB .
phx5lsta5w [ 1 ] ; tmp_i [ 54 ] = rtB . phx5lsta5w [ 2 ] ; tmp_i [ 55 ] = rtB
. phx5lsta5w [ 3 ] ; tmp_m [ 14 ] = 56 ; tmp_i [ 56 ] = rtB . hoap0zr3eu [ 0
] ; tmp_i [ 57 ] = rtB . hoap0zr3eu [ 1 ] ; tmp_i [ 58 ] = rtB . hoap0zr3eu [
2 ] ; tmp_i [ 59 ] = rtB . hoap0zr3eu [ 3 ] ; tmp_m [ 15 ] = 60 ; tmp_i [ 60
] = rtB . enp5f3s002 [ 0 ] ; tmp_i [ 61 ] = rtB . enp5f3s002 [ 1 ] ; tmp_i [
62 ] = rtB . enp5f3s002 [ 2 ] ; tmp_i [ 63 ] = rtB . enp5f3s002 [ 3 ] ; tmp_m
[ 16 ] = 64 ; tmp_i [ 64 ] = rtB . nsl31ekm0f [ 0 ] ; tmp_i [ 65 ] = rtB .
nsl31ekm0f [ 1 ] ; tmp_i [ 66 ] = rtB . nsl31ekm0f [ 2 ] ; tmp_i [ 67 ] = rtB
. nsl31ekm0f [ 3 ] ; tmp_m [ 17 ] = 68 ; tmp_i [ 68 ] = rtB . hrqvwki2jl [ 0
] ; tmp_i [ 69 ] = rtB . hrqvwki2jl [ 1 ] ; tmp_i [ 70 ] = rtB . hrqvwki2jl [
2 ] ; tmp_i [ 71 ] = rtB . hrqvwki2jl [ 3 ] ; tmp_m [ 18 ] = 72 ; tmp_i [ 72
] = rtB . dypcby02yn [ 0 ] ; tmp_i [ 73 ] = rtB . dypcby02yn [ 1 ] ; tmp_i [
74 ] = rtB . dypcby02yn [ 2 ] ; tmp_i [ 75 ] = rtB . dypcby02yn [ 3 ] ; tmp_m
[ 19 ] = 76 ; tmp_i [ 76 ] = rtB . a100glwqkz [ 0 ] ; tmp_i [ 77 ] = rtB .
a100glwqkz [ 1 ] ; tmp_i [ 78 ] = rtB . a100glwqkz [ 2 ] ; tmp_i [ 79 ] = rtB
. a100glwqkz [ 3 ] ; tmp_m [ 20 ] = 80 ; tmp_i [ 80 ] = rtB . c3twv1nzdj [ 0
] ; tmp_i [ 81 ] = rtB . c3twv1nzdj [ 1 ] ; tmp_i [ 82 ] = rtB . c3twv1nzdj [
2 ] ; tmp_i [ 83 ] = rtB . c3twv1nzdj [ 3 ] ; tmp_m [ 21 ] = 84 ; tmp_i [ 84
] = rtB . iurh4h1cla [ 0 ] ; tmp_i [ 85 ] = rtB . iurh4h1cla [ 1 ] ; tmp_i [
86 ] = rtB . iurh4h1cla [ 2 ] ; tmp_i [ 87 ] = rtB . iurh4h1cla [ 3 ] ; tmp_m
[ 22 ] = 88 ; tmp_i [ 88 ] = rtB . har0a2wrmi [ 0 ] ; tmp_i [ 89 ] = rtB .
har0a2wrmi [ 1 ] ; tmp_i [ 90 ] = rtB . har0a2wrmi [ 2 ] ; tmp_i [ 91 ] = rtB
. har0a2wrmi [ 3 ] ; tmp_m [ 23 ] = 92 ; tmp_i [ 92 ] = rtB . ib2ieuf2an [ 0
] ; tmp_i [ 93 ] = rtB . ib2ieuf2an [ 1 ] ; tmp_i [ 94 ] = rtB . ib2ieuf2an [
2 ] ; tmp_i [ 95 ] = rtB . ib2ieuf2an [ 3 ] ; tmp_m [ 24 ] = 96 ; tmp_i [ 96
] = rtB . ojpovrodvr [ 0 ] ; tmp_i [ 97 ] = rtB . ojpovrodvr [ 1 ] ; tmp_i [
98 ] = rtB . ojpovrodvr [ 2 ] ; tmp_i [ 99 ] = rtB . ojpovrodvr [ 3 ] ; tmp_m
[ 25 ] = 100 ; tmp_i [ 100 ] = rtB . fcgmpow5uz [ 0 ] ; tmp_i [ 101 ] = rtB .
fcgmpow5uz [ 1 ] ; tmp_i [ 102 ] = rtB . fcgmpow5uz [ 2 ] ; tmp_i [ 103 ] =
rtB . fcgmpow5uz [ 3 ] ; tmp_m [ 26 ] = 104 ; tmp_i [ 104 ] = rtB .
jylntbxfpj [ 0 ] ; tmp_i [ 105 ] = rtB . jylntbxfpj [ 1 ] ; tmp_i [ 106 ] =
rtB . jylntbxfpj [ 2 ] ; tmp_i [ 107 ] = rtB . jylntbxfpj [ 3 ] ; tmp_m [ 27
] = 108 ; tmp_i [ 108 ] = rtB . g1yo2p3pxt [ 0 ] ; tmp_i [ 109 ] = rtB .
g1yo2p3pxt [ 1 ] ; tmp_i [ 110 ] = rtB . g1yo2p3pxt [ 2 ] ; tmp_i [ 111 ] =
rtB . g1yo2p3pxt [ 3 ] ; tmp_m [ 28 ] = 112 ; tmp_i [ 112 ] = rtB .
az4eo3zh5o [ 0 ] ; tmp_i [ 113 ] = rtB . az4eo3zh5o [ 1 ] ; tmp_i [ 114 ] =
rtB . az4eo3zh5o [ 2 ] ; tmp_i [ 115 ] = rtB . az4eo3zh5o [ 3 ] ; tmp_m [ 29
] = 116 ; tmp_i [ 116 ] = rtB . ch1tbnhyon [ 0 ] ; tmp_i [ 117 ] = rtB .
ch1tbnhyon [ 1 ] ; tmp_i [ 118 ] = rtB . ch1tbnhyon [ 2 ] ; tmp_i [ 119 ] =
rtB . ch1tbnhyon [ 3 ] ; tmp_m [ 30 ] = 120 ; tmp_i [ 120 ] = rtB .
eoseoj2f3b [ 0 ] ; tmp_i [ 121 ] = rtB . eoseoj2f3b [ 1 ] ; tmp_i [ 122 ] =
rtB . eoseoj2f3b [ 2 ] ; tmp_i [ 123 ] = rtB . eoseoj2f3b [ 3 ] ; tmp_m [ 31
] = 124 ; tmp_i [ 124 ] = rtB . aigpe1li2p [ 0 ] ; tmp_i [ 125 ] = rtB .
aigpe1li2p [ 1 ] ; tmp_i [ 126 ] = rtB . aigpe1li2p [ 2 ] ; tmp_i [ 127 ] =
rtB . aigpe1li2p [ 3 ] ; tmp_m [ 32 ] = 128 ; tmp_i [ 128 ] = rtB .
pbn2nzk2k3 [ 0 ] ; tmp_i [ 129 ] = rtB . pbn2nzk2k3 [ 1 ] ; tmp_i [ 130 ] =
rtB . pbn2nzk2k3 [ 2 ] ; tmp_i [ 131 ] = rtB . pbn2nzk2k3 [ 3 ] ; tmp_m [ 33
] = 132 ; tmp_i [ 132 ] = rtB . hffhwq1imf [ 0 ] ; tmp_i [ 133 ] = rtB .
hffhwq1imf [ 1 ] ; tmp_i [ 134 ] = rtB . hffhwq1imf [ 2 ] ; tmp_i [ 135 ] =
rtB . hffhwq1imf [ 3 ] ; tmp_m [ 34 ] = 136 ; tmp_i [ 136 ] = rtB .
gbxuu4wtte [ 0 ] ; tmp_i [ 137 ] = rtB . gbxuu4wtte [ 1 ] ; tmp_i [ 138 ] =
rtB . gbxuu4wtte [ 2 ] ; tmp_i [ 139 ] = rtB . gbxuu4wtte [ 3 ] ; tmp_m [ 35
] = 140 ; tmp_i [ 140 ] = rtB . k5nlme3kns [ 0 ] ; tmp_i [ 141 ] = rtB .
k5nlme3kns [ 1 ] ; tmp_i [ 142 ] = rtB . k5nlme3kns [ 2 ] ; tmp_i [ 143 ] =
rtB . k5nlme3kns [ 3 ] ; tmp_m [ 36 ] = 144 ; tmp_i [ 144 ] = rtB .
j52cfj2d1w [ 0 ] ; tmp_i [ 145 ] = rtB . j52cfj2d1w [ 1 ] ; tmp_i [ 146 ] =
rtB . j52cfj2d1w [ 2 ] ; tmp_i [ 147 ] = rtB . j52cfj2d1w [ 3 ] ; tmp_m [ 37
] = 148 ; tmp_i [ 148 ] = rtB . ftvz2ekipg [ 0 ] ; tmp_i [ 149 ] = rtB .
ftvz2ekipg [ 1 ] ; tmp_i [ 150 ] = rtB . ftvz2ekipg [ 2 ] ; tmp_i [ 151 ] =
rtB . ftvz2ekipg [ 3 ] ; tmp_m [ 38 ] = 152 ; tmp_i [ 152 ] = rtB .
awgskjjmsb [ 0 ] ; tmp_i [ 153 ] = rtB . awgskjjmsb [ 1 ] ; tmp_i [ 154 ] =
rtB . awgskjjmsb [ 2 ] ; tmp_i [ 155 ] = rtB . awgskjjmsb [ 3 ] ; tmp_m [ 39
] = 156 ; tmp_i [ 156 ] = rtB . fl4rvtkarn [ 0 ] ; tmp_i [ 157 ] = rtB .
fl4rvtkarn [ 1 ] ; tmp_i [ 158 ] = rtB . fl4rvtkarn [ 2 ] ; tmp_i [ 159 ] =
rtB . fl4rvtkarn [ 3 ] ; tmp_m [ 40 ] = 160 ; tmp_i [ 160 ] = rtB .
jp22vgjwwm [ 0 ] ; tmp_i [ 161 ] = rtB . jp22vgjwwm [ 1 ] ; tmp_i [ 162 ] =
rtB . jp22vgjwwm [ 2 ] ; tmp_i [ 163 ] = rtB . jp22vgjwwm [ 3 ] ; tmp_m [ 41
] = 164 ; tmp_i [ 164 ] = rtB . pzpjv5s5n2 [ 0 ] ; tmp_i [ 165 ] = rtB .
pzpjv5s5n2 [ 1 ] ; tmp_i [ 166 ] = rtB . pzpjv5s5n2 [ 2 ] ; tmp_i [ 167 ] =
rtB . pzpjv5s5n2 [ 3 ] ; tmp_m [ 42 ] = 168 ; memcpy ( & tmp_i [ 168 ] , &
rtB . aay54cmgia [ 0 ] , 18U * sizeof ( real_T ) ) ; tmp_m [ 43 ] = 186 ;
simulationData -> mData -> mInputValues . mN = 186 ; simulationData -> mData
-> mInputValues . mX = & tmp_i [ 0 ] ; simulationData -> mData ->
mInputOffsets . mN = 44 ; simulationData -> mData -> mInputOffsets . mX = &
tmp_m [ 0 ] ; simulationData -> mData -> mOutputs . mN = 154 ; simulationData
-> mData -> mOutputs . mX = & rtB . nfyrk0gene [ 0 ] ; simulationData ->
mData -> mTolerances . mN = 0 ; simulationData -> mData -> mTolerances . mX =
NULL ; simulationData -> mData -> mCstateHasChanged = false ; time_i =
ssGetTaskTime ( rtS , 0 ) ; simulationData -> mData -> mTime . mN = 1 ;
simulationData -> mData -> mTime . mX = & time_i ; simulationData -> mData ->
mSampleHits . mN = 0 ; simulationData -> mData -> mSampleHits . mX = NULL ;
simulationData -> mData -> mIsFundamentalSampleHit = false ;
diagnosticManager = ( NeuDiagnosticManager * ) rtDW . dcyq5ws3qs ;
diagnosticTree_p = neu_diagnostic_manager_get_initial_tree (
diagnosticManager ) ; yIdx = ne_simulator_method ( ( NeslSimulator * ) rtDW .
dikp21y4fz , NESL_SIM_OUTPUTS , simulationData , diagnosticManager ) ; if (
yIdx != 0 ) { tmp = error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if (
tmp ) { msg_p = rtw_diagnostics_msg ( diagnosticTree_p ) ; ssSetErrorStatus (
rtS , msg_p ) ; } } if ( ssIsMajorTimeStep ( rtS ) && simulationData -> mData
-> mCstateHasChanged ) { ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
} for ( yIdx = 0 ; yIdx < 7 ; yIdx ++ ) { czut5pi4ll = rtP . Gain_Gain * rtB
. nfyrk0gene [ yIdx ] ; rtB . gdleddk2bu [ yIdx ] = czut5pi4ll ; rtB .
o3cnonqf2y [ yIdx ] = rtP . Gain_Gain_iby0aikld1 * czut5pi4ll ; } { if ( rtDW
. l4r0r31ubh . AQHandles && ssGetLogOutput ( rtS ) ) { sdiWriteSignal ( rtDW
. l4r0r31ubh . AQHandles , ssGetTaskTime ( rtS , 0 ) , ( char * ) & rtB .
gdleddk2bu [ 0 ] + 0 ) ; } } rtB . bjyrxtyjst [ 0 ] = rtB . nfyrk0gene [ 16 ]
; rtB . bjyrxtyjst [ 1 ] = rtB . nfyrk0gene [ 17 ] ; rtB . bjyrxtyjst [ 2 ] =
rtB . nfyrk0gene [ 18 ] ; czut5pi4ll = rtB . bjyrxtyjst [ 0 ] * rtP . x [ 7 ]
; rtB . ecnsyysan1 [ 0 ] = czut5pi4ll ; rtB . d11srgkxv1 [ 0 ] = rtB .
nfyrk0gene [ 7 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll =
rtB . bjyrxtyjst [ 1 ] * rtP . x [ 7 ] ; rtB . ecnsyysan1 [ 1 ] = czut5pi4ll
; rtB . d11srgkxv1 [ 1 ] = rtB . nfyrk0gene [ 8 ] ; cazgyvwuzb [ 1 ] =
czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB . bjyrxtyjst [ 2 ] * rtP . x [ 7 ]
; rtB . ecnsyysan1 [ 2 ] = czut5pi4ll ; rtB . d11srgkxv1 [ 2 ] = rtB .
nfyrk0gene [ 9 ] ; czut5pi4ll = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) +
czut5pi4ll * czut5pi4ll ; if ( czut5pi4ll < 0.0 ) { rtB . oywwaak55z = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) ) ; } else { rtB .
oywwaak55z = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . oxke0tvvta = ( rtB . oywwaak55z >
rtP . NormalizeVector_maxzero ) ; } if ( rtDW . oxke0tvvta ) { czut5pi4ll =
rtB . ecnsyysan1 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . ecnsyysan1 [ 1 ] ;
jj1qbnwy23_idx_2 = rtB . ecnsyysan1 [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
oywwaak55z ; } else { jj1qbnwy23_idx_2 = rtB . ecnsyysan1 [ 0 ] * 0.0 ; rtB .
itxt4ytwah [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . ecnsyysan1 [ 1 ] * 0.0 ; rtB . itxt4ytwah [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . ecnsyysan1 [ 2 ] * 0.0 ; rtB . itxt4ytwah [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_lu1rr1djht ; } rtB . n3uyejqyiv [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . n3uyejqyiv [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . n3uyejqyiv [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 19 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . dg2tz20qh3 [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
dg2tz20qh3 [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . dg2tz20qh3 [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . d11srgkxv1 [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
jd0epngxde = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . jd0epngxde < 0.0 ) { rtB . i1wwep0quu = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . jd0epngxde ) ) ; } else { rtB . i1wwep0quu =
muDoubleScalarSqrt ( rtB . jd0epngxde ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . kr4a3cr0mf = ( rtB . i1wwep0quu > rtP .
NormalizeVector_maxzero_ab1bjfgz5c ) ; } if ( rtDW . kr4a3cr0mf ) { rtB .
jytdis5y4s [ 0 ] = rtB . d11srgkxv1 [ 0 ] ; rtB . jytdis5y4s [ 1 ] = rtB .
d11srgkxv1 [ 1 ] ; rtB . jytdis5y4s [ 2 ] = rtB . d11srgkxv1 [ 2 ] ; rtB .
jytdis5y4s [ 3 ] = rtB . i1wwep0quu ; } else { czut5pi4ll = rtB . d11srgkxv1
[ 0 ] * 0.0 ; rtB . kdrism1ogn [ 0 ] = czut5pi4ll ; rtB . jytdis5y4s [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . d11srgkxv1 [ 1 ] * 0.0 ; rtB . kdrism1ogn [ 1
] = czut5pi4ll ; rtB . jytdis5y4s [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
d11srgkxv1 [ 2 ] * 0.0 ; rtB . kdrism1ogn [ 2 ] = czut5pi4ll ; rtB .
jytdis5y4s [ 2 ] = czut5pi4ll ; rtB . jytdis5y4s [ 3 ] = rtP .
Constant_Value_aldt0rz0ce ; } czut5pi4ll = rtB . jytdis5y4s [ 0 ] / rtB .
jytdis5y4s [ 3 ] ; rtB . dvlkqx1wrk [ 0 ] = czut5pi4ll ; rtB . axlew0gacz [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . jytdis5y4s [ 1 ] / rtB . jytdis5y4s [ 3 ]
; rtB . dvlkqx1wrk [ 1 ] = czut5pi4ll ; rtB . axlew0gacz [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . jytdis5y4s [ 2 ] / rtB . jytdis5y4s [ 3 ] ; rtB .
dvlkqx1wrk [ 2 ] = czut5pi4ll ; rtB . axlew0gacz [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . dvlkqx1wrk [ 0 ] ; jj1qbnwy23_idx_1 = rtB . dvlkqx1wrk [ 1
] ; jj1qbnwy23_idx_2 = rtB . dvlkqx1wrk [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . axlew0gacz [ yIdx ] ; rtB . ffk0uso1jv [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . ffk0uso1jv [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . ffk0uso1jv [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain * rtB . ffk0uso1jv [ yIdx ] ; rtB .
j2caeelflf [ yIdx ] = czut5pi4ll ; rtB . gzdxak0rgb [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value [ yIdx ] ; } czut5pi4ll = 0.0 ; jj1qbnwy23_idx_2 = rtB
. nsmqj1yv3h [ 1 ] ; jj1qbnwy23_idx_3 = rtB . nsmqj1yv3h [ 0 ] ; nsmqj1yv3h =
rtB . nsmqj1yv3h [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ; yIdx ++ ) { rtB .
nhpvoyt4rs [ yIdx ] = ( rtB . gzdxak0rgb [ yIdx + 3 ] * jj1qbnwy23_idx_2 +
rtB . gzdxak0rgb [ yIdx ] * jj1qbnwy23_idx_3 ) + rtB . gzdxak0rgb [ yIdx + 6
] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB . d11srgkxv1 [ yIdx ] ; czut5pi4ll +=
jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } if ( ssIsMajorTimeStep ( rtS ) ) { if
( rtDW . o3hz2jqijw != 0 ) { ssSetBlockStateForSolverChangedAtMajorStep ( rtS
) ; ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
o3hz2jqijw = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . o3hz2jqijw = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . mnl5wynqkw =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_ffpsqzlbdy ) ; rtB . n5zmerun4k [ 0 ]
= rtB . nhpvoyt4rs [ 0 ] * rtB . mnl5wynqkw ; rtB . n5zmerun4k [ 1 ] = rtB .
nhpvoyt4rs [ 1 ] * rtB . mnl5wynqkw ; rtB . n5zmerun4k [ 2 ] = rtB .
nhpvoyt4rs [ 2 ] * rtB . mnl5wynqkw ; rtB . i210rcby1e [ 0 ] = rtB .
ecnsyysan1 [ 1 ] * rtB . n5zmerun4k [ 2 ] ; rtB . i210rcby1e [ 1 ] = rtB .
n5zmerun4k [ 0 ] * rtB . ecnsyysan1 [ 2 ] ; rtB . i210rcby1e [ 2 ] = rtB .
ecnsyysan1 [ 0 ] * rtB . n5zmerun4k [ 1 ] ; rtB . i210rcby1e [ 3 ] = rtB .
n5zmerun4k [ 1 ] * rtB . ecnsyysan1 [ 2 ] ; rtB . i210rcby1e [ 4 ] = rtB .
ecnsyysan1 [ 0 ] * rtB . n5zmerun4k [ 2 ] ; rtB . i210rcby1e [ 5 ] = rtB .
n5zmerun4k [ 0 ] * rtB . ecnsyysan1 [ 1 ] ; rtB . njlotg3eyk [ 0 ] = rtB .
i210rcby1e [ 0 ] - rtB . i210rcby1e [ 3 ] ; rtB . njlotg3eyk [ 1 ] = rtB .
i210rcby1e [ 1 ] - rtB . i210rcby1e [ 4 ] ; rtB . njlotg3eyk [ 2 ] = rtB .
i210rcby1e [ 2 ] - rtB . i210rcby1e [ 5 ] ; czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . njlotg3eyk [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
njlotg3eyk [ 0 ] ; nsmqj1yv3h = rtB . njlotg3eyk [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 = ( rtB . dg2tz20qh3 [ yIdx + 3 ] *
jj1qbnwy23_idx_2 + rtB . dg2tz20qh3 [ yIdx ] * jj1qbnwy23_idx_3 ) + rtB .
dg2tz20qh3 [ yIdx + 6 ] * nsmqj1yv3h ; rtB . ciemy5f0t3 [ yIdx ] =
jj1qbnwy23_idx_1 ; rtB . i3fsmpf0l4 [ yIdx ] = rtP . SimParams . torque_onoff
* jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB . d11srgkxv1 [ yIdx ] ;
czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } if ( ssIsMajorTimeStep
( rtS ) ) { if ( rtDW . a1ctm05tly != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
a1ctm05tly = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . a1ctm05tly = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . p1t00nt2c4 =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value ) ; rtB . lkk1p3hkce = ( rtB .
d11srgkxv1 [ 0 ] * rtB . d11srgkxv1 [ 0 ] + rtB . d11srgkxv1 [ 1 ] * rtB .
d11srgkxv1 [ 1 ] ) + rtB . d11srgkxv1 [ 2 ] * rtB . d11srgkxv1 [ 2 ] ; if (
rtB . lkk1p3hkce < 0.0 ) { rtB . j0nzmiwajs = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . lkk1p3hkce ) ) ; } else { rtB . j0nzmiwajs =
muDoubleScalarSqrt ( rtB . lkk1p3hkce ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . gdblfyh2ib = ( rtB . j0nzmiwajs > rtP .
NormalizeVector_maxzero_jpxndxsk21 ) ; } if ( rtDW . gdblfyh2ib ) { rtB .
p4b33uxmxl [ 0 ] = rtB . d11srgkxv1 [ 0 ] ; rtB . p4b33uxmxl [ 1 ] = rtB .
d11srgkxv1 [ 1 ] ; rtB . p4b33uxmxl [ 2 ] = rtB . d11srgkxv1 [ 2 ] ; rtB .
p4b33uxmxl [ 3 ] = rtB . j0nzmiwajs ; } else { czut5pi4ll = rtB . d11srgkxv1
[ 0 ] * 0.0 ; rtB . nuca2ykly0 [ 0 ] = czut5pi4ll ; rtB . p4b33uxmxl [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . d11srgkxv1 [ 1 ] * 0.0 ; rtB . nuca2ykly0 [ 1
] = czut5pi4ll ; rtB . p4b33uxmxl [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
d11srgkxv1 [ 2 ] * 0.0 ; rtB . nuca2ykly0 [ 2 ] = czut5pi4ll ; rtB .
p4b33uxmxl [ 2 ] = czut5pi4ll ; rtB . p4b33uxmxl [ 3 ] = rtP .
Constant_Value_foidw3bu40 ; } jj1qbnwy23_idx_1 = rtB . p4b33uxmxl [ 0 ] / rtB
. p4b33uxmxl [ 3 ] ; rtB . kerf1yzlz0 [ 0 ] = jj1qbnwy23_idx_1 ; czut5pi4ll =
jj1qbnwy23_idx_1 * rtB . nsmqj1yv3h [ 0 ] ; jj1qbnwy23_idx_1 = rtB .
p4b33uxmxl [ 1 ] / rtB . p4b33uxmxl [ 3 ] ; rtB . kerf1yzlz0 [ 1 ] =
jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB . nsmqj1yv3h [ 1 ] ;
jj1qbnwy23_idx_1 = rtB . p4b33uxmxl [ 2 ] / rtB . p4b33uxmxl [ 3 ] ; rtB .
kerf1yzlz0 [ 2 ] = jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB .
nsmqj1yv3h [ 2 ] ; rtB . dxk5dr2guu [ 0 ] = czut5pi4ll * rtB . ecnsyysan1 [ 0
] ; rtB . dxk5dr2guu [ 1 ] = czut5pi4ll * rtB . ecnsyysan1 [ 1 ] ; rtB .
dxk5dr2guu [ 2 ] = czut5pi4ll * rtB . ecnsyysan1 [ 2 ] ; jj1qbnwy23_idx_1 = (
rtB . kerf1yzlz0 [ 0 ] * rtB . ecnsyysan1 [ 0 ] + rtB . kerf1yzlz0 [ 1 ] *
rtB . ecnsyysan1 [ 1 ] ) + rtB . kerf1yzlz0 [ 2 ] * rtB . ecnsyysan1 [ 2 ] ;
rtB . jmgra3bzio = czut5pi4ll * jj1qbnwy23_idx_1 * rtP . Gain_Gain_o3ussx2pey
; rtB . gfwkvsxkxh [ 0 ] = jj1qbnwy23_idx_1 * rtB . nsmqj1yv3h [ 0 ] ; rtB .
gfwkvsxkxh [ 1 ] = jj1qbnwy23_idx_1 * rtB . nsmqj1yv3h [ 1 ] ; rtB .
gfwkvsxkxh [ 2 ] = jj1qbnwy23_idx_1 * rtB . nsmqj1yv3h [ 2 ] ; rtB .
n43m55vtpu = ( ( rtB . nsmqj1yv3h [ 0 ] * rtB . ecnsyysan1 [ 0 ] + rtB .
nsmqj1yv3h [ 1 ] * rtB . ecnsyysan1 [ 1 ] ) + rtB . nsmqj1yv3h [ 2 ] * rtB .
ecnsyysan1 [ 2 ] ) - rtB . jmgra3bzio ; czut5pi4ll = rtB . kerf1yzlz0 [ 0 ] *
rtB . n43m55vtpu ; rtB . dm0nsvocoi [ 0 ] = czut5pi4ll ; czut5pi4ll += rtB .
dxk5dr2guu [ 0 ] + rtB . gfwkvsxkxh [ 0 ] ; rtB . pzoqxgrplf [ 0 ] =
czut5pi4ll ; czut5pi4ll *= rtB . p1t00nt2c4 ; rtB . pi3p3kites [ 0 ] =
czut5pi4ll ; rtB . j4tq1cqvap [ 0 ] = rtP . SimParams . force_onoff *
czut5pi4ll ; czut5pi4ll = rtB . kerf1yzlz0 [ 1 ] * rtB . n43m55vtpu ; rtB .
dm0nsvocoi [ 1 ] = czut5pi4ll ; czut5pi4ll += rtB . dxk5dr2guu [ 1 ] + rtB .
gfwkvsxkxh [ 1 ] ; rtB . pzoqxgrplf [ 1 ] = czut5pi4ll ; czut5pi4ll *= rtB .
p1t00nt2c4 ; rtB . pi3p3kites [ 1 ] = czut5pi4ll ; rtB . j4tq1cqvap [ 1 ] =
rtP . SimParams . force_onoff * czut5pi4ll ; czut5pi4ll = rtB . kerf1yzlz0 [
2 ] * rtB . n43m55vtpu ; rtB . dm0nsvocoi [ 2 ] = czut5pi4ll ; czut5pi4ll +=
rtB . dxk5dr2guu [ 2 ] + rtB . gfwkvsxkxh [ 2 ] ; rtB . pzoqxgrplf [ 2 ] =
czut5pi4ll ; czut5pi4ll *= rtB . p1t00nt2c4 ; rtB . pi3p3kites [ 2 ] =
czut5pi4ll ; rtB . j4tq1cqvap [ 2 ] = rtP . SimParams . force_onoff *
czut5pi4ll ; czut5pi4ll = ( rtB . ecnsyysan1 [ 0 ] * rtB . ecnsyysan1 [ 0 ] +
rtB . ecnsyysan1 [ 1 ] * rtB . ecnsyysan1 [ 1 ] ) + rtB . ecnsyysan1 [ 2 ] *
rtB . ecnsyysan1 [ 2 ] ; if ( czut5pi4ll < 0.0 ) { rtB . d4sl1kvckp = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) ) ; } else { rtB .
d4sl1kvckp = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . mmq5ttry5y = ( rtB . d4sl1kvckp >
rtP . NormalizeVector1_maxzero ) ; } if ( rtDW . mmq5ttry5y ) { czut5pi4ll =
rtB . ecnsyysan1 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . ecnsyysan1 [ 1 ] ;
jj1qbnwy23_idx_2 = rtB . ecnsyysan1 [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
d4sl1kvckp ; } else { jj1qbnwy23_idx_2 = rtB . ecnsyysan1 [ 0 ] * 0.0 ; rtB .
nl3crkyc5i [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . ecnsyysan1 [ 1 ] * 0.0 ; rtB . nl3crkyc5i [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . ecnsyysan1 [ 2 ] * 0.0 ; rtB . nl3crkyc5i [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_hsgo25r4sk ; } rtB . k0e5tb5t4r [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . mrse3gkxz4 [ 0 ] = rtB . nfyrk0gene [
37 ] ; rtB . k0e5tb5t4r [ 1 ] = jj1qbnwy23_idx_1 / jj1qbnwy23_idx_3 ; rtB .
mrse3gkxz4 [ 1 ] = rtB . nfyrk0gene [ 38 ] ; rtB . k0e5tb5t4r [ 2 ] =
jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB . mrse3gkxz4 [ 2 ] = rtB .
nfyrk0gene [ 39 ] ; czut5pi4ll = rtB . mrse3gkxz4 [ 0 ] * rtP . x [ 8 ] ; rtB
. hjkncigkim [ 0 ] = czut5pi4ll ; rtB . bv20qr1sq1 [ 0 ] = rtB . nfyrk0gene [
28 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB .
mrse3gkxz4 [ 1 ] * rtP . x [ 8 ] ; rtB . hjkncigkim [ 1 ] = czut5pi4ll ; rtB
. bv20qr1sq1 [ 1 ] = rtB . nfyrk0gene [ 29 ] ; cazgyvwuzb [ 1 ] = czut5pi4ll
* czut5pi4ll ; czut5pi4ll = rtB . mrse3gkxz4 [ 2 ] * rtP . x [ 8 ] ; rtB .
hjkncigkim [ 2 ] = czut5pi4ll ; rtB . bv20qr1sq1 [ 2 ] = rtB . nfyrk0gene [
30 ] ; jj1qbnwy23_idx_1 = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) +
czut5pi4ll * czut5pi4ll ; if ( jj1qbnwy23_idx_1 < 0.0 ) { rtB . jknqpcyts5 =
- muDoubleScalarSqrt ( muDoubleScalarAbs ( jj1qbnwy23_idx_1 ) ) ; } else {
rtB . jknqpcyts5 = muDoubleScalarSqrt ( jj1qbnwy23_idx_1 ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . agwigalrld = ( rtB . jknqpcyts5 >
rtP . NormalizeVector_maxzero_hach1vuaw5 ) ; } if ( rtDW . agwigalrld ) {
czut5pi4ll = rtB . hjkncigkim [ 0 ] ; jj1qbnwy23_idx_1 = rtB . hjkncigkim [ 1
] ; jj1qbnwy23_idx_2 = rtB . hjkncigkim [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
jknqpcyts5 ; } else { jj1qbnwy23_idx_2 = rtB . hjkncigkim [ 0 ] * 0.0 ; rtB .
p5hpxiz5u3 [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . hjkncigkim [ 1 ] * 0.0 ; rtB . p5hpxiz5u3 [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . hjkncigkim [ 2 ] * 0.0 ; rtB . p5hpxiz5u3 [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_kebaex1duz ; } rtB . fo44k43jik [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . fo44k43jik [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . fo44k43jik [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 40 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . leoubbnbua [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
leoubbnbua [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . leoubbnbua [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . bv20qr1sq1 [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
eb5u3a2x55 = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . eb5u3a2x55 < 0.0 ) { rtB . eded4se0e2 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . eb5u3a2x55 ) ) ; } else { rtB . eded4se0e2 =
muDoubleScalarSqrt ( rtB . eb5u3a2x55 ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . fudyiowmus = ( rtB . eded4se0e2 > rtP .
NormalizeVector_maxzero_nzh4nf0nvr ) ; } if ( rtDW . fudyiowmus ) { rtB .
epqt4glrc2 [ 0 ] = rtB . bv20qr1sq1 [ 0 ] ; rtB . epqt4glrc2 [ 1 ] = rtB .
bv20qr1sq1 [ 1 ] ; rtB . epqt4glrc2 [ 2 ] = rtB . bv20qr1sq1 [ 2 ] ; rtB .
epqt4glrc2 [ 3 ] = rtB . eded4se0e2 ; } else { czut5pi4ll = rtB . bv20qr1sq1
[ 0 ] * 0.0 ; rtB . jqlckttaam [ 0 ] = czut5pi4ll ; rtB . epqt4glrc2 [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . bv20qr1sq1 [ 1 ] * 0.0 ; rtB . jqlckttaam [ 1
] = czut5pi4ll ; rtB . epqt4glrc2 [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
bv20qr1sq1 [ 2 ] * 0.0 ; rtB . jqlckttaam [ 2 ] = czut5pi4ll ; rtB .
epqt4glrc2 [ 2 ] = czut5pi4ll ; rtB . epqt4glrc2 [ 3 ] = rtP .
Constant_Value_ag31g0q0wg ; } czut5pi4ll = rtB . epqt4glrc2 [ 0 ] / rtB .
epqt4glrc2 [ 3 ] ; rtB . huyyr23tvs [ 0 ] = czut5pi4ll ; rtB . a5iwna5ofw [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . epqt4glrc2 [ 1 ] / rtB . epqt4glrc2 [ 3 ]
; rtB . huyyr23tvs [ 1 ] = czut5pi4ll ; rtB . a5iwna5ofw [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . epqt4glrc2 [ 2 ] / rtB . epqt4glrc2 [ 3 ] ; rtB .
huyyr23tvs [ 2 ] = czut5pi4ll ; rtB . a5iwna5ofw [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . huyyr23tvs [ 0 ] ; jj1qbnwy23_idx_1 = rtB . huyyr23tvs [ 1
] ; jj1qbnwy23_idx_2 = rtB . huyyr23tvs [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . a5iwna5ofw [ yIdx ] ; rtB . n5qkyzyc0p [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . n5qkyzyc0p [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . n5qkyzyc0p [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain_inzuw1ftmc * rtB . n5qkyzyc0p [ yIdx ] ; rtB
. cviiadazje [ yIdx ] = czut5pi4ll ; rtB . hrdkh50xl2 [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value_mat3dcwilq [ yIdx ] ; } czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . mdrxplc5qx [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
mdrxplc5qx [ 0 ] ; nsmqj1yv3h = rtB . mdrxplc5qx [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { rtB . p4bi154apu [ yIdx ] = ( rtB . hrdkh50xl2 [ yIdx
+ 3 ] * jj1qbnwy23_idx_2 + rtB . hrdkh50xl2 [ yIdx ] * jj1qbnwy23_idx_3 ) +
rtB . hrdkh50xl2 [ yIdx + 6 ] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB .
bv20qr1sq1 [ yIdx ] ; czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; }
if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . eiqm0ndfyg != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
eiqm0ndfyg = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . eiqm0ndfyg = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . ollqza00rh =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_ccdbgnkv2m ) ; rtB . kypr3ufomq [ 0 ]
= rtB . p4bi154apu [ 0 ] * rtB . ollqza00rh ; rtB . kypr3ufomq [ 1 ] = rtB .
p4bi154apu [ 1 ] * rtB . ollqza00rh ; rtB . kypr3ufomq [ 2 ] = rtB .
p4bi154apu [ 2 ] * rtB . ollqza00rh ; rtB . ljklcjruel [ 0 ] = rtB .
hjkncigkim [ 1 ] * rtB . kypr3ufomq [ 2 ] ; rtB . ljklcjruel [ 1 ] = rtB .
kypr3ufomq [ 0 ] * rtB . hjkncigkim [ 2 ] ; rtB . ljklcjruel [ 2 ] = rtB .
hjkncigkim [ 0 ] * rtB . kypr3ufomq [ 1 ] ; rtB . ljklcjruel [ 3 ] = rtB .
kypr3ufomq [ 1 ] * rtB . hjkncigkim [ 2 ] ; rtB . ljklcjruel [ 4 ] = rtB .
hjkncigkim [ 0 ] * rtB . kypr3ufomq [ 2 ] ; rtB . ljklcjruel [ 5 ] = rtB .
kypr3ufomq [ 0 ] * rtB . hjkncigkim [ 1 ] ; rtB . giyavk0qxh [ 0 ] = rtB .
ljklcjruel [ 0 ] - rtB . ljklcjruel [ 3 ] ; rtB . giyavk0qxh [ 1 ] = rtB .
ljklcjruel [ 1 ] - rtB . ljklcjruel [ 4 ] ; rtB . giyavk0qxh [ 2 ] = rtB .
ljklcjruel [ 2 ] - rtB . ljklcjruel [ 5 ] ; czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . giyavk0qxh [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
giyavk0qxh [ 0 ] ; nsmqj1yv3h = rtB . giyavk0qxh [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 = ( rtB . leoubbnbua [ yIdx + 3 ] *
jj1qbnwy23_idx_2 + rtB . leoubbnbua [ yIdx ] * jj1qbnwy23_idx_3 ) + rtB .
leoubbnbua [ yIdx + 6 ] * nsmqj1yv3h ; rtB . iwp0nxhxfg [ yIdx ] =
jj1qbnwy23_idx_1 ; rtB . cwwkrda1qf [ yIdx ] = rtP . SimParams . torque_onoff
* jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB . bv20qr1sq1 [ yIdx ] ;
czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } if ( ssIsMajorTimeStep
( rtS ) ) { if ( rtDW . pibyff0eql != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
pibyff0eql = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . pibyff0eql = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . ojzjawyp2p =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_hrhqyti5xa ) ; rtB . j22ln0sz4g = (
rtB . bv20qr1sq1 [ 0 ] * rtB . bv20qr1sq1 [ 0 ] + rtB . bv20qr1sq1 [ 1 ] *
rtB . bv20qr1sq1 [ 1 ] ) + rtB . bv20qr1sq1 [ 2 ] * rtB . bv20qr1sq1 [ 2 ] ;
if ( rtB . j22ln0sz4g < 0.0 ) { rtB . ltesgv5ex5 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . j22ln0sz4g ) ) ; } else { rtB . ltesgv5ex5 =
muDoubleScalarSqrt ( rtB . j22ln0sz4g ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . ftrjzsptcb = ( rtB . ltesgv5ex5 > rtP .
NormalizeVector_maxzero_pjti1ch2qi ) ; } if ( rtDW . ftrjzsptcb ) { rtB .
e4k3bfcn1z [ 0 ] = rtB . bv20qr1sq1 [ 0 ] ; rtB . e4k3bfcn1z [ 1 ] = rtB .
bv20qr1sq1 [ 1 ] ; rtB . e4k3bfcn1z [ 2 ] = rtB . bv20qr1sq1 [ 2 ] ; rtB .
e4k3bfcn1z [ 3 ] = rtB . ltesgv5ex5 ; } else { czut5pi4ll = rtB . bv20qr1sq1
[ 0 ] * 0.0 ; rtB . npob1uuqas [ 0 ] = czut5pi4ll ; rtB . e4k3bfcn1z [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . bv20qr1sq1 [ 1 ] * 0.0 ; rtB . npob1uuqas [ 1
] = czut5pi4ll ; rtB . e4k3bfcn1z [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
bv20qr1sq1 [ 2 ] * 0.0 ; rtB . npob1uuqas [ 2 ] = czut5pi4ll ; rtB .
e4k3bfcn1z [ 2 ] = czut5pi4ll ; rtB . e4k3bfcn1z [ 3 ] = rtP .
Constant_Value_bgbfcrykss ; } jj1qbnwy23_idx_1 = rtB . e4k3bfcn1z [ 0 ] / rtB
. e4k3bfcn1z [ 3 ] ; rtB . fi20a5zfv5 [ 0 ] = jj1qbnwy23_idx_1 ; czut5pi4ll =
jj1qbnwy23_idx_1 * rtB . mdrxplc5qx [ 0 ] ; jj1qbnwy23_idx_1 = rtB .
e4k3bfcn1z [ 1 ] / rtB . e4k3bfcn1z [ 3 ] ; rtB . fi20a5zfv5 [ 1 ] =
jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB . mdrxplc5qx [ 1 ] ;
jj1qbnwy23_idx_1 = rtB . e4k3bfcn1z [ 2 ] / rtB . e4k3bfcn1z [ 3 ] ; rtB .
fi20a5zfv5 [ 2 ] = jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB .
mdrxplc5qx [ 2 ] ; rtB . oexwmjcbqo [ 0 ] = czut5pi4ll * rtB . hjkncigkim [ 0
] ; rtB . oexwmjcbqo [ 1 ] = czut5pi4ll * rtB . hjkncigkim [ 1 ] ; rtB .
oexwmjcbqo [ 2 ] = czut5pi4ll * rtB . hjkncigkim [ 2 ] ; jj1qbnwy23_idx_1 = (
rtB . fi20a5zfv5 [ 0 ] * rtB . hjkncigkim [ 0 ] + rtB . fi20a5zfv5 [ 1 ] *
rtB . hjkncigkim [ 1 ] ) + rtB . fi20a5zfv5 [ 2 ] * rtB . hjkncigkim [ 2 ] ;
rtB . nj4zxwlnaz = czut5pi4ll * jj1qbnwy23_idx_1 * rtP . Gain_Gain_ch2juginai
; rtB . ljdiomnr0e [ 0 ] = jj1qbnwy23_idx_1 * rtB . mdrxplc5qx [ 0 ] ; rtB .
ljdiomnr0e [ 1 ] = jj1qbnwy23_idx_1 * rtB . mdrxplc5qx [ 1 ] ; rtB .
ljdiomnr0e [ 2 ] = jj1qbnwy23_idx_1 * rtB . mdrxplc5qx [ 2 ] ; rtB .
b4eypdlrbg = ( ( rtB . mdrxplc5qx [ 0 ] * rtB . hjkncigkim [ 0 ] + rtB .
mdrxplc5qx [ 1 ] * rtB . hjkncigkim [ 1 ] ) + rtB . mdrxplc5qx [ 2 ] * rtB .
hjkncigkim [ 2 ] ) - rtB . nj4zxwlnaz ; czut5pi4ll = rtB . fi20a5zfv5 [ 0 ] *
rtB . b4eypdlrbg ; rtB . mswvoafvyk [ 0 ] = czut5pi4ll ; czut5pi4ll += rtB .
oexwmjcbqo [ 0 ] + rtB . ljdiomnr0e [ 0 ] ; rtB . ga3a4xdbbc [ 0 ] =
czut5pi4ll ; rtB . eecvbm2tyj [ 0 ] = rtB . ojzjawyp2p * czut5pi4ll ;
czut5pi4ll = rtB . fi20a5zfv5 [ 1 ] * rtB . b4eypdlrbg ; rtB . mswvoafvyk [ 1
] = czut5pi4ll ; czut5pi4ll += rtB . oexwmjcbqo [ 1 ] + rtB . ljdiomnr0e [ 1
] ; rtB . ga3a4xdbbc [ 1 ] = czut5pi4ll ; rtB . eecvbm2tyj [ 1 ] = rtB .
ojzjawyp2p * czut5pi4ll ; czut5pi4ll = rtB . fi20a5zfv5 [ 2 ] * rtB .
b4eypdlrbg ; rtB . mswvoafvyk [ 2 ] = czut5pi4ll ; czut5pi4ll += rtB .
oexwmjcbqo [ 2 ] + rtB . ljdiomnr0e [ 2 ] ; rtB . ga3a4xdbbc [ 2 ] =
czut5pi4ll ; rtB . eecvbm2tyj [ 2 ] = rtB . ojzjawyp2p * czut5pi4ll ;
czut5pi4ll = ( rtB . hjkncigkim [ 0 ] * rtB . hjkncigkim [ 0 ] + rtB .
hjkncigkim [ 1 ] * rtB . hjkncigkim [ 1 ] ) + rtB . hjkncigkim [ 2 ] * rtB .
hjkncigkim [ 2 ] ; if ( czut5pi4ll < 0.0 ) { rtB . bkxjvl2wwl = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) ) ; } else { rtB .
bkxjvl2wwl = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . ciu1ngvzgz = ( rtB . bkxjvl2wwl >
rtP . NormalizeVector1_maxzero_hfyr12og40 ) ; } if ( rtDW . ciu1ngvzgz ) {
czut5pi4ll = rtB . hjkncigkim [ 0 ] ; jj1qbnwy23_idx_1 = rtB . hjkncigkim [ 1
] ; jj1qbnwy23_idx_2 = rtB . hjkncigkim [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
bkxjvl2wwl ; } else { jj1qbnwy23_idx_2 = rtB . hjkncigkim [ 0 ] * 0.0 ; rtB .
h2ica3jpzg [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . hjkncigkim [ 1 ] * 0.0 ; rtB . h2ica3jpzg [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . hjkncigkim [ 2 ] * 0.0 ; rtB . h2ica3jpzg [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_lpo3xrcqjh ; } rtB . jcdegmljub [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . dakqbtmvna [ 0 ] = rtP . SimParams .
force_onoff * rtB . eecvbm2tyj [ 0 ] ; rtB . j54ujtxdmj [ 0 ] = rtB .
nfyrk0gene [ 58 ] ; rtB . jcdegmljub [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . dakqbtmvna [ 1 ] = rtP . SimParams . force_onoff *
rtB . eecvbm2tyj [ 1 ] ; rtB . j54ujtxdmj [ 1 ] = rtB . nfyrk0gene [ 59 ] ;
rtB . jcdegmljub [ 2 ] = jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB .
dakqbtmvna [ 2 ] = rtP . SimParams . force_onoff * rtB . eecvbm2tyj [ 2 ] ;
rtB . j54ujtxdmj [ 2 ] = rtB . nfyrk0gene [ 60 ] ; czut5pi4ll = rtB .
j54ujtxdmj [ 0 ] * rtP . x [ 9 ] ; rtB . ol0mge5uht [ 0 ] = czut5pi4ll ; rtB
. fusdzrgc4y [ 0 ] = rtB . nfyrk0gene [ 49 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll
* czut5pi4ll ; czut5pi4ll = rtB . j54ujtxdmj [ 1 ] * rtP . x [ 9 ] ; rtB .
ol0mge5uht [ 1 ] = czut5pi4ll ; rtB . fusdzrgc4y [ 1 ] = rtB . nfyrk0gene [
50 ] ; cazgyvwuzb [ 1 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB .
j54ujtxdmj [ 2 ] * rtP . x [ 9 ] ; rtB . ol0mge5uht [ 2 ] = czut5pi4ll ; rtB
. fusdzrgc4y [ 2 ] = rtB . nfyrk0gene [ 51 ] ; jj1qbnwy23_idx_1 = (
cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + czut5pi4ll * czut5pi4ll ; if (
jj1qbnwy23_idx_1 < 0.0 ) { rtB . n2cw123fsj = - muDoubleScalarSqrt (
muDoubleScalarAbs ( jj1qbnwy23_idx_1 ) ) ; } else { rtB . n2cw123fsj =
muDoubleScalarSqrt ( jj1qbnwy23_idx_1 ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . csdsboi20g = ( rtB . n2cw123fsj > rtP .
NormalizeVector_maxzero_jvlgwkux2b ) ; } if ( rtDW . csdsboi20g ) {
czut5pi4ll = rtB . ol0mge5uht [ 0 ] ; jj1qbnwy23_idx_1 = rtB . ol0mge5uht [ 1
] ; jj1qbnwy23_idx_2 = rtB . ol0mge5uht [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
n2cw123fsj ; } else { jj1qbnwy23_idx_2 = rtB . ol0mge5uht [ 0 ] * 0.0 ; rtB .
inhs5zftab [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . ol0mge5uht [ 1 ] * 0.0 ; rtB . inhs5zftab [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . ol0mge5uht [ 2 ] * 0.0 ; rtB . inhs5zftab [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_mqmlcee2pu ; } rtB . gtc4x0vmrj [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . gtc4x0vmrj [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . gtc4x0vmrj [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 61 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . igdgu1b2hk [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
igdgu1b2hk [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . igdgu1b2hk [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . fusdzrgc4y [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
lntlwwhs4g = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . lntlwwhs4g < 0.0 ) { rtB . lvldddhyjb = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . lntlwwhs4g ) ) ; } else { rtB . lvldddhyjb =
muDoubleScalarSqrt ( rtB . lntlwwhs4g ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . pxsktqyoas = ( rtB . lvldddhyjb > rtP .
NormalizeVector_maxzero_k3zfgywmre ) ; } if ( rtDW . pxsktqyoas ) { rtB .
j1izfy2apg [ 0 ] = rtB . fusdzrgc4y [ 0 ] ; rtB . j1izfy2apg [ 1 ] = rtB .
fusdzrgc4y [ 1 ] ; rtB . j1izfy2apg [ 2 ] = rtB . fusdzrgc4y [ 2 ] ; rtB .
j1izfy2apg [ 3 ] = rtB . lvldddhyjb ; } else { czut5pi4ll = rtB . fusdzrgc4y
[ 0 ] * 0.0 ; rtB . cschdiniut [ 0 ] = czut5pi4ll ; rtB . j1izfy2apg [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . fusdzrgc4y [ 1 ] * 0.0 ; rtB . cschdiniut [ 1
] = czut5pi4ll ; rtB . j1izfy2apg [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
fusdzrgc4y [ 2 ] * 0.0 ; rtB . cschdiniut [ 2 ] = czut5pi4ll ; rtB .
j1izfy2apg [ 2 ] = czut5pi4ll ; rtB . j1izfy2apg [ 3 ] = rtP .
Constant_Value_a4fqp5ywdt ; } czut5pi4ll = rtB . j1izfy2apg [ 0 ] / rtB .
j1izfy2apg [ 3 ] ; rtB . f2e0m2twoo [ 0 ] = czut5pi4ll ; rtB . pctih4a2gv [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . j1izfy2apg [ 1 ] / rtB . j1izfy2apg [ 3 ]
; rtB . f2e0m2twoo [ 1 ] = czut5pi4ll ; rtB . pctih4a2gv [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . j1izfy2apg [ 2 ] / rtB . j1izfy2apg [ 3 ] ; rtB .
f2e0m2twoo [ 2 ] = czut5pi4ll ; rtB . pctih4a2gv [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . f2e0m2twoo [ 0 ] ; jj1qbnwy23_idx_1 = rtB . f2e0m2twoo [ 1
] ; jj1qbnwy23_idx_2 = rtB . f2e0m2twoo [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . pctih4a2gv [ yIdx ] ; rtB . d15mjte2mm [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . d15mjte2mm [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . d15mjte2mm [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain_iijf52nd0p * rtB . d15mjte2mm [ yIdx ] ; rtB
. o0sl3javac [ yIdx ] = czut5pi4ll ; rtB . lk2g50sct1 [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value_obh4dra0ft [ yIdx ] ; } czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . idcqc0uibn [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
idcqc0uibn [ 0 ] ; nsmqj1yv3h = rtB . idcqc0uibn [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { rtB . kucavif4r4 [ yIdx ] = ( rtB . lk2g50sct1 [ yIdx
+ 3 ] * jj1qbnwy23_idx_2 + rtB . lk2g50sct1 [ yIdx ] * jj1qbnwy23_idx_3 ) +
rtB . lk2g50sct1 [ yIdx + 6 ] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB .
fusdzrgc4y [ yIdx ] ; czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; }
if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . l4r1nqhb1g != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
l4r1nqhb1g = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . l4r1nqhb1g = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . fd02dku0jj =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_pxe2jyk15f ) ; rtB . j3daqzj5si [ 0 ]
= rtB . kucavif4r4 [ 0 ] * rtB . fd02dku0jj ; rtB . j3daqzj5si [ 1 ] = rtB .
kucavif4r4 [ 1 ] * rtB . fd02dku0jj ; rtB . j3daqzj5si [ 2 ] = rtB .
kucavif4r4 [ 2 ] * rtB . fd02dku0jj ; rtB . e30x0ooogi [ 0 ] = rtB .
ol0mge5uht [ 1 ] * rtB . j3daqzj5si [ 2 ] ; rtB . e30x0ooogi [ 1 ] = rtB .
j3daqzj5si [ 0 ] * rtB . ol0mge5uht [ 2 ] ; rtB . e30x0ooogi [ 2 ] = rtB .
ol0mge5uht [ 0 ] * rtB . j3daqzj5si [ 1 ] ; rtB . e30x0ooogi [ 3 ] = rtB .
j3daqzj5si [ 1 ] * rtB . ol0mge5uht [ 2 ] ; rtB . e30x0ooogi [ 4 ] = rtB .
ol0mge5uht [ 0 ] * rtB . j3daqzj5si [ 2 ] ; rtB . e30x0ooogi [ 5 ] = rtB .
j3daqzj5si [ 0 ] * rtB . ol0mge5uht [ 1 ] ; rtB . dgv155wtql [ 0 ] = rtB .
e30x0ooogi [ 0 ] - rtB . e30x0ooogi [ 3 ] ; rtB . dgv155wtql [ 1 ] = rtB .
e30x0ooogi [ 1 ] - rtB . e30x0ooogi [ 4 ] ; rtB . dgv155wtql [ 2 ] = rtB .
e30x0ooogi [ 2 ] - rtB . e30x0ooogi [ 5 ] ; czut5pi4ll = rtB . dgv155wtql [ 1
] ; jj1qbnwy23_idx_2 = rtB . dgv155wtql [ 0 ] ; jj1qbnwy23_idx_3 = rtB .
dgv155wtql [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 =
( rtB . igdgu1b2hk [ yIdx + 3 ] * czut5pi4ll + rtB . igdgu1b2hk [ yIdx ] *
jj1qbnwy23_idx_2 ) + rtB . igdgu1b2hk [ yIdx + 6 ] * jj1qbnwy23_idx_3 ; rtB .
bcueeadxdi [ yIdx ] = jj1qbnwy23_idx_1 ; rtB . nalsmrc4um [ yIdx ] = rtP .
SimParams . torque_onoff * jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB .
fusdzrgc4y [ yIdx ] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 *
jj1qbnwy23_idx_1 ; } rtB . gmj04jp0xd = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ]
) + cazgyvwuzb [ 2 ] ; if ( rtB . gmj04jp0xd < 0.0 ) { rtB . ij2rjc05nz = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( rtB . gmj04jp0xd ) ) ; } else { rtB
. ij2rjc05nz = muDoubleScalarSqrt ( rtB . gmj04jp0xd ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . bya2xiv5j3 = ( rtB . ij2rjc05nz >
rtP . NormalizeVector_maxzero_jceokergvn ) ; } if ( rtDW . bya2xiv5j3 ) { rtB
. pyxwiqcuay [ 0 ] = rtB . fusdzrgc4y [ 0 ] ; rtB . pyxwiqcuay [ 1 ] = rtB .
fusdzrgc4y [ 1 ] ; rtB . pyxwiqcuay [ 2 ] = rtB . fusdzrgc4y [ 2 ] ; rtB .
pyxwiqcuay [ 3 ] = rtB . ij2rjc05nz ; } else { czut5pi4ll = rtB . fusdzrgc4y
[ 0 ] * 0.0 ; rtB . o4wni3515w [ 0 ] = czut5pi4ll ; rtB . pyxwiqcuay [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . fusdzrgc4y [ 1 ] * 0.0 ; rtB . o4wni3515w [ 1
] = czut5pi4ll ; rtB . pyxwiqcuay [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
fusdzrgc4y [ 2 ] * 0.0 ; rtB . o4wni3515w [ 2 ] = czut5pi4ll ; rtB .
pyxwiqcuay [ 2 ] = czut5pi4ll ; rtB . pyxwiqcuay [ 3 ] = rtP .
Constant_Value_p4qexwk3wd ; } jj1qbnwy23_idx_2 = rtB . pyxwiqcuay [ 0 ] / rtB
. pyxwiqcuay [ 3 ] ; rtB . ndcqoyp0vx [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll =
jj1qbnwy23_idx_2 * rtB . idcqc0uibn [ 0 ] ; jj1qbnwy23_idx_1 =
jj1qbnwy23_idx_2 * rtB . ol0mge5uht [ 0 ] ; jj1qbnwy23_idx_2 = rtB .
pyxwiqcuay [ 1 ] / rtB . pyxwiqcuay [ 3 ] ; rtB . ndcqoyp0vx [ 1 ] =
jj1qbnwy23_idx_2 ; czut5pi4ll += jj1qbnwy23_idx_2 * rtB . idcqc0uibn [ 1 ] ;
jj1qbnwy23_idx_1 += jj1qbnwy23_idx_2 * rtB . ol0mge5uht [ 1 ] ;
jj1qbnwy23_idx_2 = rtB . pyxwiqcuay [ 2 ] / rtB . pyxwiqcuay [ 3 ] ; rtB .
ndcqoyp0vx [ 2 ] = jj1qbnwy23_idx_2 ; czut5pi4ll += jj1qbnwy23_idx_2 * rtB .
idcqc0uibn [ 2 ] ; jj1qbnwy23_idx_1 += jj1qbnwy23_idx_2 * rtB . ol0mge5uht [
2 ] ; rtB . hv30h5tsxp = czut5pi4ll * jj1qbnwy23_idx_1 * rtP .
Gain_Gain_cnz04octod ; rtB . iqlb144ibi = ( ( rtB . idcqc0uibn [ 0 ] * rtB .
ol0mge5uht [ 0 ] + rtB . idcqc0uibn [ 1 ] * rtB . ol0mge5uht [ 1 ] ) + rtB .
idcqc0uibn [ 2 ] * rtB . ol0mge5uht [ 2 ] ) - rtB . hv30h5tsxp ;
jj1qbnwy23_idx_2 = rtB . ndcqoyp0vx [ 0 ] * rtB . iqlb144ibi ; rtB .
pl32irduyc [ 0 ] = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 = czut5pi4ll * rtB .
ol0mge5uht [ 0 ] ; rtB . lgifpnih0e [ 0 ] = jj1qbnwy23_idx_3 ; nsmqj1yv3h =
jj1qbnwy23_idx_1 * rtB . idcqc0uibn [ 0 ] ; rtB . l0eiwonrmi [ 0 ] =
nsmqj1yv3h ; rtB . hrtkouuc3i [ 0 ] = ( jj1qbnwy23_idx_3 + nsmqj1yv3h ) +
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 = rtB . ndcqoyp0vx [ 1 ] * rtB .
iqlb144ibi ; rtB . pl32irduyc [ 1 ] = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 =
czut5pi4ll * rtB . ol0mge5uht [ 1 ] ; rtB . lgifpnih0e [ 1 ] =
jj1qbnwy23_idx_3 ; nsmqj1yv3h = jj1qbnwy23_idx_1 * rtB . idcqc0uibn [ 1 ] ;
rtB . l0eiwonrmi [ 1 ] = nsmqj1yv3h ; rtB . hrtkouuc3i [ 1 ] = (
jj1qbnwy23_idx_3 + nsmqj1yv3h ) + jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 = rtB .
ndcqoyp0vx [ 2 ] * rtB . iqlb144ibi ; rtB . pl32irduyc [ 2 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 = czut5pi4ll * rtB . ol0mge5uht [ 2 ] ;
rtB . lgifpnih0e [ 2 ] = jj1qbnwy23_idx_3 ; nsmqj1yv3h = jj1qbnwy23_idx_1 *
rtB . idcqc0uibn [ 2 ] ; rtB . l0eiwonrmi [ 2 ] = nsmqj1yv3h ; rtB .
hrtkouuc3i [ 2 ] = ( jj1qbnwy23_idx_3 + nsmqj1yv3h ) + jj1qbnwy23_idx_2 ;
czut5pi4ll = ( rtB . fusdzrgc4y [ 0 ] * rtB . fusdzrgc4y [ 0 ] + rtB .
fusdzrgc4y [ 1 ] * rtB . fusdzrgc4y [ 1 ] ) + rtB . fusdzrgc4y [ 2 ] * rtB .
fusdzrgc4y [ 2 ] ; if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . emuswyjmhs
!= 0 ) { ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
emuswyjmhs = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . emuswyjmhs = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . gmo4itlisv =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_crcwx2inuv ) ; rtB . nexxpkogi1 [ 0 ]
= rtB . gmo4itlisv * rtB . hrtkouuc3i [ 0 ] ; rtB . nexxpkogi1 [ 1 ] = rtB .
gmo4itlisv * rtB . hrtkouuc3i [ 1 ] ; rtB . nexxpkogi1 [ 2 ] = rtB .
gmo4itlisv * rtB . hrtkouuc3i [ 2 ] ; czut5pi4ll = ( rtB . ol0mge5uht [ 0 ] *
rtB . ol0mge5uht [ 0 ] + rtB . ol0mge5uht [ 1 ] * rtB . ol0mge5uht [ 1 ] ) +
rtB . ol0mge5uht [ 2 ] * rtB . ol0mge5uht [ 2 ] ; if ( czut5pi4ll < 0.0 ) {
rtB . idwieepnv4 = - muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) )
; } else { rtB . idwieepnv4 = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . hnl0z3aet1 = ( rtB . idwieepnv4 >
rtP . NormalizeVector1_maxzero_edxcmuydlr ) ; } if ( rtDW . hnl0z3aet1 ) {
czut5pi4ll = rtB . ol0mge5uht [ 0 ] ; jj1qbnwy23_idx_1 = rtB . ol0mge5uht [ 1
] ; jj1qbnwy23_idx_2 = rtB . ol0mge5uht [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
idwieepnv4 ; } else { jj1qbnwy23_idx_2 = rtB . ol0mge5uht [ 0 ] * 0.0 ; rtB .
hrudfjbz3i [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . ol0mge5uht [ 1 ] * 0.0 ; rtB . hrudfjbz3i [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . ol0mge5uht [ 2 ] * 0.0 ; rtB . hrudfjbz3i [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_lwhq5q2fck ; } rtB . h2vbgw30tf [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . dndiwk2b1s [ 0 ] = rtP . SimParams .
force_onoff * rtB . nexxpkogi1 [ 0 ] ; rtB . gvavv3jbs0 [ 0 ] = rtB .
nfyrk0gene [ 79 ] ; rtB . h2vbgw30tf [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . dndiwk2b1s [ 1 ] = rtP . SimParams . force_onoff *
rtB . nexxpkogi1 [ 1 ] ; rtB . gvavv3jbs0 [ 1 ] = rtB . nfyrk0gene [ 80 ] ;
rtB . h2vbgw30tf [ 2 ] = jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB .
dndiwk2b1s [ 2 ] = rtP . SimParams . force_onoff * rtB . nexxpkogi1 [ 2 ] ;
rtB . gvavv3jbs0 [ 2 ] = rtB . nfyrk0gene [ 81 ] ; czut5pi4ll = rtB .
gvavv3jbs0 [ 0 ] * rtP . x [ 10 ] ; rtB . dw12uizmd2 [ 0 ] = czut5pi4ll ; rtB
. bdswxd3b20 [ 0 ] = rtB . nfyrk0gene [ 70 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll
* czut5pi4ll ; czut5pi4ll = rtB . gvavv3jbs0 [ 1 ] * rtP . x [ 10 ] ; rtB .
dw12uizmd2 [ 1 ] = czut5pi4ll ; rtB . bdswxd3b20 [ 1 ] = rtB . nfyrk0gene [
71 ] ; cazgyvwuzb [ 1 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB .
gvavv3jbs0 [ 2 ] * rtP . x [ 10 ] ; rtB . dw12uizmd2 [ 2 ] = czut5pi4ll ; rtB
. bdswxd3b20 [ 2 ] = rtB . nfyrk0gene [ 72 ] ; jj1qbnwy23_idx_1 = (
cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + czut5pi4ll * czut5pi4ll ; if (
jj1qbnwy23_idx_1 < 0.0 ) { rtB . fh04srz5g0 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( jj1qbnwy23_idx_1 ) ) ; } else { rtB . fh04srz5g0 =
muDoubleScalarSqrt ( jj1qbnwy23_idx_1 ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . ecvugwr44v = ( rtB . fh04srz5g0 > rtP .
NormalizeVector_maxzero_levzkkbbsj ) ; } if ( rtDW . ecvugwr44v ) {
czut5pi4ll = rtB . dw12uizmd2 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . dw12uizmd2 [ 1
] ; jj1qbnwy23_idx_2 = rtB . dw12uizmd2 [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
fh04srz5g0 ; } else { jj1qbnwy23_idx_2 = rtB . dw12uizmd2 [ 0 ] * 0.0 ; rtB .
itx5qkpva4 [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . dw12uizmd2 [ 1 ] * 0.0 ; rtB . itx5qkpva4 [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . dw12uizmd2 [ 2 ] * 0.0 ; rtB . itx5qkpva4 [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_iinhm22lbs ; } rtB . mkzi43umox [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . mkzi43umox [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . mkzi43umox [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 82 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . cpzvyphms1 [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
cpzvyphms1 [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . cpzvyphms1 [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . bdswxd3b20 [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
oc1jrzdnhe = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . oc1jrzdnhe < 0.0 ) { rtB . hlnjhscgik = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . oc1jrzdnhe ) ) ; } else { rtB . hlnjhscgik =
muDoubleScalarSqrt ( rtB . oc1jrzdnhe ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . dfdi1srria = ( rtB . hlnjhscgik > rtP .
NormalizeVector_maxzero_kt51ehrppa ) ; } if ( rtDW . dfdi1srria ) { rtB .
fkzq14ohjb [ 0 ] = rtB . bdswxd3b20 [ 0 ] ; rtB . fkzq14ohjb [ 1 ] = rtB .
bdswxd3b20 [ 1 ] ; rtB . fkzq14ohjb [ 2 ] = rtB . bdswxd3b20 [ 2 ] ; rtB .
fkzq14ohjb [ 3 ] = rtB . hlnjhscgik ; } else { czut5pi4ll = rtB . bdswxd3b20
[ 0 ] * 0.0 ; rtB . mdxnlbsfvm [ 0 ] = czut5pi4ll ; rtB . fkzq14ohjb [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . bdswxd3b20 [ 1 ] * 0.0 ; rtB . mdxnlbsfvm [ 1
] = czut5pi4ll ; rtB . fkzq14ohjb [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
bdswxd3b20 [ 2 ] * 0.0 ; rtB . mdxnlbsfvm [ 2 ] = czut5pi4ll ; rtB .
fkzq14ohjb [ 2 ] = czut5pi4ll ; rtB . fkzq14ohjb [ 3 ] = rtP .
Constant_Value_hooldt5mbh ; } czut5pi4ll = rtB . fkzq14ohjb [ 0 ] / rtB .
fkzq14ohjb [ 3 ] ; rtB . cafqksoz00 [ 0 ] = czut5pi4ll ; rtB . c402o0md5n [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . fkzq14ohjb [ 1 ] / rtB . fkzq14ohjb [ 3 ]
; rtB . cafqksoz00 [ 1 ] = czut5pi4ll ; rtB . c402o0md5n [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . fkzq14ohjb [ 2 ] / rtB . fkzq14ohjb [ 3 ] ; rtB .
cafqksoz00 [ 2 ] = czut5pi4ll ; rtB . c402o0md5n [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . cafqksoz00 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . cafqksoz00 [ 1
] ; jj1qbnwy23_idx_2 = rtB . cafqksoz00 [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . c402o0md5n [ yIdx ] ; rtB . hwppgfjqwy [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . hwppgfjqwy [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . hwppgfjqwy [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain_kgqmwgbxmo * rtB . hwppgfjqwy [ yIdx ] ; rtB
. i1lnvldtbv [ yIdx ] = czut5pi4ll ; rtB . ahw1dydhgo [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value_htaifltrf1 [ yIdx ] ; } czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . lkwtwvkutl [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
lkwtwvkutl [ 0 ] ; nsmqj1yv3h = rtB . lkwtwvkutl [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { rtB . ort4y5zlaf [ yIdx ] = ( rtB . ahw1dydhgo [ yIdx
+ 3 ] * jj1qbnwy23_idx_2 + rtB . ahw1dydhgo [ yIdx ] * jj1qbnwy23_idx_3 ) +
rtB . ahw1dydhgo [ yIdx + 6 ] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB .
bdswxd3b20 [ yIdx ] ; czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; }
if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . ogqlv3h0ke != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
ogqlv3h0ke = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . ogqlv3h0ke = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . f1arg5nwa3 =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_kyzvzubrmq ) ; rtB . eckq3h2tzb [ 0 ]
= rtB . ort4y5zlaf [ 0 ] * rtB . f1arg5nwa3 ; rtB . eckq3h2tzb [ 1 ] = rtB .
ort4y5zlaf [ 1 ] * rtB . f1arg5nwa3 ; rtB . eckq3h2tzb [ 2 ] = rtB .
ort4y5zlaf [ 2 ] * rtB . f1arg5nwa3 ; rtB . bsaqv0u04c [ 0 ] = rtB .
dw12uizmd2 [ 1 ] * rtB . eckq3h2tzb [ 2 ] ; rtB . bsaqv0u04c [ 1 ] = rtB .
eckq3h2tzb [ 0 ] * rtB . dw12uizmd2 [ 2 ] ; rtB . bsaqv0u04c [ 2 ] = rtB .
dw12uizmd2 [ 0 ] * rtB . eckq3h2tzb [ 1 ] ; rtB . bsaqv0u04c [ 3 ] = rtB .
eckq3h2tzb [ 1 ] * rtB . dw12uizmd2 [ 2 ] ; rtB . bsaqv0u04c [ 4 ] = rtB .
dw12uizmd2 [ 0 ] * rtB . eckq3h2tzb [ 2 ] ; rtB . bsaqv0u04c [ 5 ] = rtB .
eckq3h2tzb [ 0 ] * rtB . dw12uizmd2 [ 1 ] ; rtB . alzuotumkg [ 0 ] = rtB .
bsaqv0u04c [ 0 ] - rtB . bsaqv0u04c [ 3 ] ; rtB . alzuotumkg [ 1 ] = rtB .
bsaqv0u04c [ 1 ] - rtB . bsaqv0u04c [ 4 ] ; rtB . alzuotumkg [ 2 ] = rtB .
bsaqv0u04c [ 2 ] - rtB . bsaqv0u04c [ 5 ] ; czut5pi4ll = rtB . alzuotumkg [ 1
] ; jj1qbnwy23_idx_2 = rtB . alzuotumkg [ 0 ] ; jj1qbnwy23_idx_3 = rtB .
alzuotumkg [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 =
( rtB . cpzvyphms1 [ yIdx + 3 ] * czut5pi4ll + rtB . cpzvyphms1 [ yIdx ] *
jj1qbnwy23_idx_2 ) + rtB . cpzvyphms1 [ yIdx + 6 ] * jj1qbnwy23_idx_3 ; rtB .
ley50ghvfj [ yIdx ] = jj1qbnwy23_idx_1 ; rtB . b0dwpzqyrm [ yIdx ] = rtP .
SimParams . torque_onoff * jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB .
bdswxd3b20 [ yIdx ] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 *
jj1qbnwy23_idx_1 ; } rtB . fbdf2tllgu = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ]
) + cazgyvwuzb [ 2 ] ; if ( rtB . fbdf2tllgu < 0.0 ) { rtB . my44vuapi0 = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( rtB . fbdf2tllgu ) ) ; } else { rtB
. my44vuapi0 = muDoubleScalarSqrt ( rtB . fbdf2tllgu ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . ibwargsoog = ( rtB . my44vuapi0 >
rtP . NormalizeVector_maxzero_hbfzo5qsh2 ) ; } if ( rtDW . ibwargsoog ) { rtB
. k3faf5atq0 [ 0 ] = rtB . bdswxd3b20 [ 0 ] ; rtB . k3faf5atq0 [ 1 ] = rtB .
bdswxd3b20 [ 1 ] ; rtB . k3faf5atq0 [ 2 ] = rtB . bdswxd3b20 [ 2 ] ; rtB .
k3faf5atq0 [ 3 ] = rtB . my44vuapi0 ; } else { czut5pi4ll = rtB . bdswxd3b20
[ 0 ] * 0.0 ; rtB . fvtfdqea2r [ 0 ] = czut5pi4ll ; rtB . k3faf5atq0 [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . bdswxd3b20 [ 1 ] * 0.0 ; rtB . fvtfdqea2r [ 1
] = czut5pi4ll ; rtB . k3faf5atq0 [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
bdswxd3b20 [ 2 ] * 0.0 ; rtB . fvtfdqea2r [ 2 ] = czut5pi4ll ; rtB .
k3faf5atq0 [ 2 ] = czut5pi4ll ; rtB . k3faf5atq0 [ 3 ] = rtP .
Constant_Value_f4bgccfy2v ; } jj1qbnwy23_idx_2 = rtB . k3faf5atq0 [ 0 ] / rtB
. k3faf5atq0 [ 3 ] ; rtB . cg3eiq4i4s [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll =
jj1qbnwy23_idx_2 * rtB . lkwtwvkutl [ 0 ] ; jj1qbnwy23_idx_1 =
jj1qbnwy23_idx_2 * rtB . dw12uizmd2 [ 0 ] ; jj1qbnwy23_idx_2 = rtB .
k3faf5atq0 [ 1 ] / rtB . k3faf5atq0 [ 3 ] ; rtB . cg3eiq4i4s [ 1 ] =
jj1qbnwy23_idx_2 ; czut5pi4ll += jj1qbnwy23_idx_2 * rtB . lkwtwvkutl [ 1 ] ;
jj1qbnwy23_idx_1 += jj1qbnwy23_idx_2 * rtB . dw12uizmd2 [ 1 ] ;
jj1qbnwy23_idx_2 = rtB . k3faf5atq0 [ 2 ] / rtB . k3faf5atq0 [ 3 ] ; rtB .
cg3eiq4i4s [ 2 ] = jj1qbnwy23_idx_2 ; czut5pi4ll += jj1qbnwy23_idx_2 * rtB .
lkwtwvkutl [ 2 ] ; jj1qbnwy23_idx_1 += jj1qbnwy23_idx_2 * rtB . dw12uizmd2 [
2 ] ; rtB . encqrxuu3y = czut5pi4ll * jj1qbnwy23_idx_1 * rtP .
Gain_Gain_mxfyehzrbo ; rtB . bdozk1fpny = ( ( rtB . lkwtwvkutl [ 0 ] * rtB .
dw12uizmd2 [ 0 ] + rtB . lkwtwvkutl [ 1 ] * rtB . dw12uizmd2 [ 1 ] ) + rtB .
lkwtwvkutl [ 2 ] * rtB . dw12uizmd2 [ 2 ] ) - rtB . encqrxuu3y ;
jj1qbnwy23_idx_2 = rtB . cg3eiq4i4s [ 0 ] * rtB . bdozk1fpny ; rtB .
nbrrv2d3zg [ 0 ] = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 = czut5pi4ll * rtB .
dw12uizmd2 [ 0 ] ; rtB . cuxk4j3qys [ 0 ] = jj1qbnwy23_idx_3 ; nsmqj1yv3h =
jj1qbnwy23_idx_1 * rtB . lkwtwvkutl [ 0 ] ; rtB . eegu2ofo43 [ 0 ] =
nsmqj1yv3h ; rtB . hacw1qj3aj [ 0 ] = ( jj1qbnwy23_idx_3 + nsmqj1yv3h ) +
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 = rtB . cg3eiq4i4s [ 1 ] * rtB .
bdozk1fpny ; rtB . nbrrv2d3zg [ 1 ] = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 =
czut5pi4ll * rtB . dw12uizmd2 [ 1 ] ; rtB . cuxk4j3qys [ 1 ] =
jj1qbnwy23_idx_3 ; nsmqj1yv3h = jj1qbnwy23_idx_1 * rtB . lkwtwvkutl [ 1 ] ;
rtB . eegu2ofo43 [ 1 ] = nsmqj1yv3h ; rtB . hacw1qj3aj [ 1 ] = (
jj1qbnwy23_idx_3 + nsmqj1yv3h ) + jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 = rtB .
cg3eiq4i4s [ 2 ] * rtB . bdozk1fpny ; rtB . nbrrv2d3zg [ 2 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 = czut5pi4ll * rtB . dw12uizmd2 [ 2 ] ;
rtB . cuxk4j3qys [ 2 ] = jj1qbnwy23_idx_3 ; nsmqj1yv3h = jj1qbnwy23_idx_1 *
rtB . lkwtwvkutl [ 2 ] ; rtB . eegu2ofo43 [ 2 ] = nsmqj1yv3h ; rtB .
hacw1qj3aj [ 2 ] = ( jj1qbnwy23_idx_3 + nsmqj1yv3h ) + jj1qbnwy23_idx_2 ;
czut5pi4ll = ( rtB . bdswxd3b20 [ 0 ] * rtB . bdswxd3b20 [ 0 ] + rtB .
bdswxd3b20 [ 1 ] * rtB . bdswxd3b20 [ 1 ] ) + rtB . bdswxd3b20 [ 2 ] * rtB .
bdswxd3b20 [ 2 ] ; if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . hgiqtwygy1
!= 0 ) { ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
hgiqtwygy1 = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . hgiqtwygy1 = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . pxg5b2slaw =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_mulpe1a5eh ) ; rtB . h4h1ivsmcr [ 0 ]
= rtB . pxg5b2slaw * rtB . hacw1qj3aj [ 0 ] ; rtB . h4h1ivsmcr [ 1 ] = rtB .
pxg5b2slaw * rtB . hacw1qj3aj [ 1 ] ; rtB . h4h1ivsmcr [ 2 ] = rtB .
pxg5b2slaw * rtB . hacw1qj3aj [ 2 ] ; czut5pi4ll = ( rtB . dw12uizmd2 [ 0 ] *
rtB . dw12uizmd2 [ 0 ] + rtB . dw12uizmd2 [ 1 ] * rtB . dw12uizmd2 [ 1 ] ) +
rtB . dw12uizmd2 [ 2 ] * rtB . dw12uizmd2 [ 2 ] ; if ( czut5pi4ll < 0.0 ) {
rtB . jqzltcnumt = - muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) )
; } else { rtB . jqzltcnumt = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . cb3tpkjx5q = ( rtB . jqzltcnumt >
rtP . NormalizeVector1_maxzero_os1qknjcxt ) ; } if ( rtDW . cb3tpkjx5q ) {
czut5pi4ll = rtB . dw12uizmd2 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . dw12uizmd2 [ 1
] ; jj1qbnwy23_idx_2 = rtB . dw12uizmd2 [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
jqzltcnumt ; } else { jj1qbnwy23_idx_2 = rtB . dw12uizmd2 [ 0 ] * 0.0 ; rtB .
e4tk4xaa2e [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . dw12uizmd2 [ 1 ] * 0.0 ; rtB . e4tk4xaa2e [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . dw12uizmd2 [ 2 ] * 0.0 ; rtB . e4tk4xaa2e [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_fyrymlyvwp ; } rtB . etawd0whn4 [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . cbmxpbd2iu [ 0 ] = rtP . SimParams .
force_onoff * rtB . h4h1ivsmcr [ 0 ] ; rtB . l1wqaovs4e [ 0 ] = rtB .
nfyrk0gene [ 100 ] ; rtB . etawd0whn4 [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . cbmxpbd2iu [ 1 ] = rtP . SimParams . force_onoff *
rtB . h4h1ivsmcr [ 1 ] ; rtB . l1wqaovs4e [ 1 ] = rtB . nfyrk0gene [ 101 ] ;
rtB . etawd0whn4 [ 2 ] = jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB .
cbmxpbd2iu [ 2 ] = rtP . SimParams . force_onoff * rtB . h4h1ivsmcr [ 2 ] ;
rtB . l1wqaovs4e [ 2 ] = rtB . nfyrk0gene [ 102 ] ; czut5pi4ll = rtB .
l1wqaovs4e [ 0 ] * rtP . x [ 11 ] ; rtB . f34tk41pk1 [ 0 ] = czut5pi4ll ; rtB
. kup2makpml [ 0 ] = rtB . nfyrk0gene [ 91 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll
* czut5pi4ll ; czut5pi4ll = rtB . l1wqaovs4e [ 1 ] * rtP . x [ 11 ] ; rtB .
f34tk41pk1 [ 1 ] = czut5pi4ll ; rtB . kup2makpml [ 1 ] = rtB . nfyrk0gene [
92 ] ; cazgyvwuzb [ 1 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB .
l1wqaovs4e [ 2 ] * rtP . x [ 11 ] ; rtB . f34tk41pk1 [ 2 ] = czut5pi4ll ; rtB
. kup2makpml [ 2 ] = rtB . nfyrk0gene [ 93 ] ; jj1qbnwy23_idx_1 = (
cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + czut5pi4ll * czut5pi4ll ; if (
jj1qbnwy23_idx_1 < 0.0 ) { rtB . d4it13dar5 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( jj1qbnwy23_idx_1 ) ) ; } else { rtB . d4it13dar5 =
muDoubleScalarSqrt ( jj1qbnwy23_idx_1 ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . fg4y55crvv = ( rtB . d4it13dar5 > rtP .
NormalizeVector_maxzero_anlse3xi4h ) ; } if ( rtDW . fg4y55crvv ) {
czut5pi4ll = rtB . f34tk41pk1 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . f34tk41pk1 [ 1
] ; jj1qbnwy23_idx_2 = rtB . f34tk41pk1 [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
d4it13dar5 ; } else { jj1qbnwy23_idx_2 = rtB . f34tk41pk1 [ 0 ] * 0.0 ; rtB .
jltrdwyzmw [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . f34tk41pk1 [ 1 ] * 0.0 ; rtB . jltrdwyzmw [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . f34tk41pk1 [ 2 ] * 0.0 ; rtB . jltrdwyzmw [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_e34xedeirm ; } rtB . dl4tl2eodk [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . dl4tl2eodk [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . dl4tl2eodk [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 103 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . obk5xjqd1g [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
obk5xjqd1g [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . obk5xjqd1g [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . kup2makpml [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
p00fanezjt = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . p00fanezjt < 0.0 ) { rtB . b4b10et3s3 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . p00fanezjt ) ) ; } else { rtB . b4b10et3s3 =
muDoubleScalarSqrt ( rtB . p00fanezjt ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . h4wxpouwz2 = ( rtB . b4b10et3s3 > rtP .
NormalizeVector_maxzero_dfin15pekf ) ; } if ( rtDW . h4wxpouwz2 ) { rtB .
flzta24bjr [ 0 ] = rtB . kup2makpml [ 0 ] ; rtB . flzta24bjr [ 1 ] = rtB .
kup2makpml [ 1 ] ; rtB . flzta24bjr [ 2 ] = rtB . kup2makpml [ 2 ] ; rtB .
flzta24bjr [ 3 ] = rtB . b4b10et3s3 ; } else { czut5pi4ll = rtB . kup2makpml
[ 0 ] * 0.0 ; rtB . lebemiabnp [ 0 ] = czut5pi4ll ; rtB . flzta24bjr [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . kup2makpml [ 1 ] * 0.0 ; rtB . lebemiabnp [ 1
] = czut5pi4ll ; rtB . flzta24bjr [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
kup2makpml [ 2 ] * 0.0 ; rtB . lebemiabnp [ 2 ] = czut5pi4ll ; rtB .
flzta24bjr [ 2 ] = czut5pi4ll ; rtB . flzta24bjr [ 3 ] = rtP .
Constant_Value_gxy2qdgtgj ; } czut5pi4ll = rtB . flzta24bjr [ 0 ] / rtB .
flzta24bjr [ 3 ] ; rtB . dv2svmenmy [ 0 ] = czut5pi4ll ; rtB . odlamrju54 [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . flzta24bjr [ 1 ] / rtB . flzta24bjr [ 3 ]
; rtB . dv2svmenmy [ 1 ] = czut5pi4ll ; rtB . odlamrju54 [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . flzta24bjr [ 2 ] / rtB . flzta24bjr [ 3 ] ; rtB .
dv2svmenmy [ 2 ] = czut5pi4ll ; rtB . odlamrju54 [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . dv2svmenmy [ 0 ] ; jj1qbnwy23_idx_1 = rtB . dv2svmenmy [ 1
] ; jj1qbnwy23_idx_2 = rtB . dv2svmenmy [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . odlamrju54 [ yIdx ] ; rtB . p3ru4cddgi [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . p3ru4cddgi [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . p3ru4cddgi [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain_pklggewqk2 * rtB . p3ru4cddgi [ yIdx ] ; rtB
. bqjv3qidt0 [ yIdx ] = czut5pi4ll ; rtB . lkoo4fzueb [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value_ajigyq4ohn [ yIdx ] ; } czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . dunjuwexfg [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
dunjuwexfg [ 0 ] ; nsmqj1yv3h = rtB . dunjuwexfg [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { rtB . gxha2ko2vz [ yIdx ] = ( rtB . lkoo4fzueb [ yIdx
+ 3 ] * jj1qbnwy23_idx_2 + rtB . lkoo4fzueb [ yIdx ] * jj1qbnwy23_idx_3 ) +
rtB . lkoo4fzueb [ yIdx + 6 ] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB .
kup2makpml [ yIdx ] ; czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; }
if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . a1er4qtcuo != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
a1er4qtcuo = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . a1er4qtcuo = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . lffwcdbemm =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_c50rewkplm ) ; rtB . chig3jpwbn [ 0 ]
= rtB . gxha2ko2vz [ 0 ] * rtB . lffwcdbemm ; rtB . chig3jpwbn [ 1 ] = rtB .
gxha2ko2vz [ 1 ] * rtB . lffwcdbemm ; rtB . chig3jpwbn [ 2 ] = rtB .
gxha2ko2vz [ 2 ] * rtB . lffwcdbemm ; rtB . fcsl1vvaio [ 0 ] = rtB .
f34tk41pk1 [ 1 ] * rtB . chig3jpwbn [ 2 ] ; rtB . fcsl1vvaio [ 1 ] = rtB .
chig3jpwbn [ 0 ] * rtB . f34tk41pk1 [ 2 ] ; rtB . fcsl1vvaio [ 2 ] = rtB .
f34tk41pk1 [ 0 ] * rtB . chig3jpwbn [ 1 ] ; rtB . fcsl1vvaio [ 3 ] = rtB .
chig3jpwbn [ 1 ] * rtB . f34tk41pk1 [ 2 ] ; rtB . fcsl1vvaio [ 4 ] = rtB .
f34tk41pk1 [ 0 ] * rtB . chig3jpwbn [ 2 ] ; rtB . fcsl1vvaio [ 5 ] = rtB .
chig3jpwbn [ 0 ] * rtB . f34tk41pk1 [ 1 ] ; rtB . pu15hoqemb [ 0 ] = rtB .
fcsl1vvaio [ 0 ] - rtB . fcsl1vvaio [ 3 ] ; rtB . pu15hoqemb [ 1 ] = rtB .
fcsl1vvaio [ 1 ] - rtB . fcsl1vvaio [ 4 ] ; rtB . pu15hoqemb [ 2 ] = rtB .
fcsl1vvaio [ 2 ] - rtB . fcsl1vvaio [ 5 ] ; czut5pi4ll = rtB . pu15hoqemb [ 1
] ; jj1qbnwy23_idx_2 = rtB . pu15hoqemb [ 0 ] ; jj1qbnwy23_idx_3 = rtB .
pu15hoqemb [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 =
( rtB . obk5xjqd1g [ yIdx + 3 ] * czut5pi4ll + rtB . obk5xjqd1g [ yIdx ] *
jj1qbnwy23_idx_2 ) + rtB . obk5xjqd1g [ yIdx + 6 ] * jj1qbnwy23_idx_3 ; rtB .
emdpdwkzg3 [ yIdx ] = jj1qbnwy23_idx_1 ; rtB . lxihdocxga [ yIdx ] = rtP .
SimParams . torque_onoff * jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB .
kup2makpml [ yIdx ] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 *
jj1qbnwy23_idx_1 ; } rtB . nfak5a1vn4 = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ]
) + cazgyvwuzb [ 2 ] ; if ( rtB . nfak5a1vn4 < 0.0 ) { rtB . ntt11gc034 = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( rtB . nfak5a1vn4 ) ) ; } else { rtB
. ntt11gc034 = muDoubleScalarSqrt ( rtB . nfak5a1vn4 ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . loncelqhha = ( rtB . ntt11gc034 >
rtP . NormalizeVector_maxzero_nfiirapz00 ) ; } if ( rtDW . loncelqhha ) { rtB
. hvd2sf2hag [ 0 ] = rtB . kup2makpml [ 0 ] ; rtB . hvd2sf2hag [ 1 ] = rtB .
kup2makpml [ 1 ] ; rtB . hvd2sf2hag [ 2 ] = rtB . kup2makpml [ 2 ] ; rtB .
hvd2sf2hag [ 3 ] = rtB . ntt11gc034 ; } else { czut5pi4ll = rtB . kup2makpml
[ 0 ] * 0.0 ; rtB . p5ihfgrspp [ 0 ] = czut5pi4ll ; rtB . hvd2sf2hag [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . kup2makpml [ 1 ] * 0.0 ; rtB . p5ihfgrspp [ 1
] = czut5pi4ll ; rtB . hvd2sf2hag [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
kup2makpml [ 2 ] * 0.0 ; rtB . p5ihfgrspp [ 2 ] = czut5pi4ll ; rtB .
hvd2sf2hag [ 2 ] = czut5pi4ll ; rtB . hvd2sf2hag [ 3 ] = rtP .
Constant_Value_fkwr40m5ih ; } jj1qbnwy23_idx_2 = rtB . hvd2sf2hag [ 0 ] / rtB
. hvd2sf2hag [ 3 ] ; rtB . bpewjnssvt [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll =
jj1qbnwy23_idx_2 * rtB . dunjuwexfg [ 0 ] ; jj1qbnwy23_idx_1 =
jj1qbnwy23_idx_2 * rtB . f34tk41pk1 [ 0 ] ; jj1qbnwy23_idx_2 = rtB .
hvd2sf2hag [ 1 ] / rtB . hvd2sf2hag [ 3 ] ; rtB . bpewjnssvt [ 1 ] =
jj1qbnwy23_idx_2 ; czut5pi4ll += jj1qbnwy23_idx_2 * rtB . dunjuwexfg [ 1 ] ;
jj1qbnwy23_idx_1 += jj1qbnwy23_idx_2 * rtB . f34tk41pk1 [ 1 ] ;
jj1qbnwy23_idx_2 = rtB . hvd2sf2hag [ 2 ] / rtB . hvd2sf2hag [ 3 ] ; rtB .
bpewjnssvt [ 2 ] = jj1qbnwy23_idx_2 ; czut5pi4ll += jj1qbnwy23_idx_2 * rtB .
dunjuwexfg [ 2 ] ; jj1qbnwy23_idx_1 += jj1qbnwy23_idx_2 * rtB . f34tk41pk1 [
2 ] ; rtB . j1ydxe2usi = czut5pi4ll * jj1qbnwy23_idx_1 * rtP .
Gain_Gain_ch2qj2nra0 ; rtB . ceb1rzrjkf = ( ( rtB . dunjuwexfg [ 0 ] * rtB .
f34tk41pk1 [ 0 ] + rtB . dunjuwexfg [ 1 ] * rtB . f34tk41pk1 [ 1 ] ) + rtB .
dunjuwexfg [ 2 ] * rtB . f34tk41pk1 [ 2 ] ) - rtB . j1ydxe2usi ;
jj1qbnwy23_idx_2 = rtB . bpewjnssvt [ 0 ] * rtB . ceb1rzrjkf ; rtB .
lyhkcvfx5d [ 0 ] = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 = czut5pi4ll * rtB .
f34tk41pk1 [ 0 ] ; rtB . enoycngmy3 [ 0 ] = jj1qbnwy23_idx_3 ; nsmqj1yv3h =
jj1qbnwy23_idx_1 * rtB . dunjuwexfg [ 0 ] ; rtB . mw4hhpohry [ 0 ] =
nsmqj1yv3h ; rtB . bhsuwz2d2l [ 0 ] = ( jj1qbnwy23_idx_3 + nsmqj1yv3h ) +
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 = rtB . bpewjnssvt [ 1 ] * rtB .
ceb1rzrjkf ; rtB . lyhkcvfx5d [ 1 ] = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 =
czut5pi4ll * rtB . f34tk41pk1 [ 1 ] ; rtB . enoycngmy3 [ 1 ] =
jj1qbnwy23_idx_3 ; nsmqj1yv3h = jj1qbnwy23_idx_1 * rtB . dunjuwexfg [ 1 ] ;
rtB . mw4hhpohry [ 1 ] = nsmqj1yv3h ; rtB . bhsuwz2d2l [ 1 ] = (
jj1qbnwy23_idx_3 + nsmqj1yv3h ) + jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 = rtB .
bpewjnssvt [ 2 ] * rtB . ceb1rzrjkf ; rtB . lyhkcvfx5d [ 2 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_3 = czut5pi4ll * rtB . f34tk41pk1 [ 2 ] ;
rtB . enoycngmy3 [ 2 ] = jj1qbnwy23_idx_3 ; nsmqj1yv3h = jj1qbnwy23_idx_1 *
rtB . dunjuwexfg [ 2 ] ; rtB . mw4hhpohry [ 2 ] = nsmqj1yv3h ; rtB .
bhsuwz2d2l [ 2 ] = ( jj1qbnwy23_idx_3 + nsmqj1yv3h ) + jj1qbnwy23_idx_2 ;
czut5pi4ll = ( rtB . kup2makpml [ 0 ] * rtB . kup2makpml [ 0 ] + rtB .
kup2makpml [ 1 ] * rtB . kup2makpml [ 1 ] ) + rtB . kup2makpml [ 2 ] * rtB .
kup2makpml [ 2 ] ; if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . jtymq0r3wm
!= 0 ) { ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
jtymq0r3wm = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . jtymq0r3wm = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . ezvzwoavxd =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_o54qqmia1e ) ; rtB . msqcnwfh2s [ 0 ]
= rtB . ezvzwoavxd * rtB . bhsuwz2d2l [ 0 ] ; rtB . msqcnwfh2s [ 1 ] = rtB .
ezvzwoavxd * rtB . bhsuwz2d2l [ 1 ] ; rtB . msqcnwfh2s [ 2 ] = rtB .
ezvzwoavxd * rtB . bhsuwz2d2l [ 2 ] ; czut5pi4ll = ( rtB . f34tk41pk1 [ 0 ] *
rtB . f34tk41pk1 [ 0 ] + rtB . f34tk41pk1 [ 1 ] * rtB . f34tk41pk1 [ 1 ] ) +
rtB . f34tk41pk1 [ 2 ] * rtB . f34tk41pk1 [ 2 ] ; if ( czut5pi4ll < 0.0 ) {
rtB . mp054zeo2w = - muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) )
; } else { rtB . mp054zeo2w = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . afd0kmksnl = ( rtB . mp054zeo2w >
rtP . NormalizeVector1_maxzero_pbowfqqdx5 ) ; } if ( rtDW . afd0kmksnl ) {
czut5pi4ll = rtB . f34tk41pk1 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . f34tk41pk1 [ 1
] ; jj1qbnwy23_idx_2 = rtB . f34tk41pk1 [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
mp054zeo2w ; } else { jj1qbnwy23_idx_2 = rtB . f34tk41pk1 [ 0 ] * 0.0 ; rtB .
fc15cwanga [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . f34tk41pk1 [ 1 ] * 0.0 ; rtB . fc15cwanga [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . f34tk41pk1 [ 2 ] * 0.0 ; rtB . fc15cwanga [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_oqtflizduh ; } rtB . pfosbwya4i [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . l01qgkxm4q [ 0 ] = rtP . SimParams .
force_onoff * rtB . msqcnwfh2s [ 0 ] ; rtB . babxk0cwue [ 0 ] = rtB .
nfyrk0gene [ 121 ] ; rtB . pfosbwya4i [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . l01qgkxm4q [ 1 ] = rtP . SimParams . force_onoff *
rtB . msqcnwfh2s [ 1 ] ; rtB . babxk0cwue [ 1 ] = rtB . nfyrk0gene [ 122 ] ;
rtB . pfosbwya4i [ 2 ] = jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB .
l01qgkxm4q [ 2 ] = rtP . SimParams . force_onoff * rtB . msqcnwfh2s [ 2 ] ;
rtB . babxk0cwue [ 2 ] = rtB . nfyrk0gene [ 123 ] ; czut5pi4ll = rtB .
babxk0cwue [ 0 ] * rtP . x [ 12 ] ; rtB . fzbgmnan4p [ 0 ] = czut5pi4ll ; rtB
. oflmpiaez1 [ 0 ] = rtB . nfyrk0gene [ 112 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll
* czut5pi4ll ; czut5pi4ll = rtB . babxk0cwue [ 1 ] * rtP . x [ 12 ] ; rtB .
fzbgmnan4p [ 1 ] = czut5pi4ll ; rtB . oflmpiaez1 [ 1 ] = rtB . nfyrk0gene [
113 ] ; cazgyvwuzb [ 1 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB .
babxk0cwue [ 2 ] * rtP . x [ 12 ] ; rtB . fzbgmnan4p [ 2 ] = czut5pi4ll ; rtB
. oflmpiaez1 [ 2 ] = rtB . nfyrk0gene [ 114 ] ; jj1qbnwy23_idx_1 = (
cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + czut5pi4ll * czut5pi4ll ; if (
jj1qbnwy23_idx_1 < 0.0 ) { rtB . efmfmwthqv = - muDoubleScalarSqrt (
muDoubleScalarAbs ( jj1qbnwy23_idx_1 ) ) ; } else { rtB . efmfmwthqv =
muDoubleScalarSqrt ( jj1qbnwy23_idx_1 ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . pv0jmcqzdg = ( rtB . efmfmwthqv > rtP .
NormalizeVector_maxzero_ap3h53mgxx ) ; } if ( rtDW . pv0jmcqzdg ) {
czut5pi4ll = rtB . fzbgmnan4p [ 0 ] ; jj1qbnwy23_idx_1 = rtB . fzbgmnan4p [ 1
] ; jj1qbnwy23_idx_2 = rtB . fzbgmnan4p [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
efmfmwthqv ; } else { jj1qbnwy23_idx_2 = rtB . fzbgmnan4p [ 0 ] * 0.0 ; rtB .
l02fn5ndrj [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . fzbgmnan4p [ 1 ] * 0.0 ; rtB . l02fn5ndrj [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . fzbgmnan4p [ 2 ] * 0.0 ; rtB . l02fn5ndrj [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_gimcsalf5s ; } rtB . ivaqfjshfj [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . ivaqfjshfj [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . ivaqfjshfj [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 124 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . fzqka3pwrj [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
fzqka3pwrj [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . fzqka3pwrj [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . oflmpiaez1 [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
egvr4ozujh = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . egvr4ozujh < 0.0 ) { rtB . cjrpkdpszw = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . egvr4ozujh ) ) ; } else { rtB . cjrpkdpszw =
muDoubleScalarSqrt ( rtB . egvr4ozujh ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . lbzvytb4d4 = ( rtB . cjrpkdpszw > rtP .
NormalizeVector_maxzero_kytveykcqy ) ; } if ( rtDW . lbzvytb4d4 ) { rtB .
logpc2hca5 [ 0 ] = rtB . oflmpiaez1 [ 0 ] ; rtB . logpc2hca5 [ 1 ] = rtB .
oflmpiaez1 [ 1 ] ; rtB . logpc2hca5 [ 2 ] = rtB . oflmpiaez1 [ 2 ] ; rtB .
logpc2hca5 [ 3 ] = rtB . cjrpkdpszw ; } else { czut5pi4ll = rtB . oflmpiaez1
[ 0 ] * 0.0 ; rtB . d1zqfcz4zo [ 0 ] = czut5pi4ll ; rtB . logpc2hca5 [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . oflmpiaez1 [ 1 ] * 0.0 ; rtB . d1zqfcz4zo [ 1
] = czut5pi4ll ; rtB . logpc2hca5 [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
oflmpiaez1 [ 2 ] * 0.0 ; rtB . d1zqfcz4zo [ 2 ] = czut5pi4ll ; rtB .
logpc2hca5 [ 2 ] = czut5pi4ll ; rtB . logpc2hca5 [ 3 ] = rtP .
Constant_Value_ndmgh4xmik ; } czut5pi4ll = rtB . logpc2hca5 [ 0 ] / rtB .
logpc2hca5 [ 3 ] ; rtB . ksrswfmi3x [ 0 ] = czut5pi4ll ; rtB . d2uj0bf0jq [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . logpc2hca5 [ 1 ] / rtB . logpc2hca5 [ 3 ]
; rtB . ksrswfmi3x [ 1 ] = czut5pi4ll ; rtB . d2uj0bf0jq [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . logpc2hca5 [ 2 ] / rtB . logpc2hca5 [ 3 ] ; rtB .
ksrswfmi3x [ 2 ] = czut5pi4ll ; rtB . d2uj0bf0jq [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . ksrswfmi3x [ 0 ] ; jj1qbnwy23_idx_1 = rtB . ksrswfmi3x [ 1
] ; jj1qbnwy23_idx_2 = rtB . ksrswfmi3x [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . d2uj0bf0jq [ yIdx ] ; rtB . joimmmlrrj [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . joimmmlrrj [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . joimmmlrrj [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain_mt5ivn2ky5 * rtB . joimmmlrrj [ yIdx ] ; rtB
. aaw4gtn4yc [ yIdx ] = czut5pi4ll ; rtB . jnz14mmcr3 [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value_dagdlneepy [ yIdx ] ; } czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . oimhqzeh2t [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
oimhqzeh2t [ 0 ] ; nsmqj1yv3h = rtB . oimhqzeh2t [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { rtB . le0lvwzojd [ yIdx ] = ( rtB . jnz14mmcr3 [ yIdx
+ 3 ] * jj1qbnwy23_idx_2 + rtB . jnz14mmcr3 [ yIdx ] * jj1qbnwy23_idx_3 ) +
rtB . jnz14mmcr3 [ yIdx + 6 ] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB .
oflmpiaez1 [ yIdx ] ; czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; }
if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . bj4ejpvtiv != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
bj4ejpvtiv = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . bj4ejpvtiv = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . k3vxq1m344 =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_pik4yifivl ) ; rtB . e3hmodlujb [ 0 ]
= rtB . le0lvwzojd [ 0 ] * rtB . k3vxq1m344 ; rtB . e3hmodlujb [ 1 ] = rtB .
le0lvwzojd [ 1 ] * rtB . k3vxq1m344 ; rtB . e3hmodlujb [ 2 ] = rtB .
le0lvwzojd [ 2 ] * rtB . k3vxq1m344 ; rtB . b3yujzbge2 [ 0 ] = rtB .
fzbgmnan4p [ 1 ] * rtB . e3hmodlujb [ 2 ] ; rtB . b3yujzbge2 [ 1 ] = rtB .
e3hmodlujb [ 0 ] * rtB . fzbgmnan4p [ 2 ] ; rtB . b3yujzbge2 [ 2 ] = rtB .
fzbgmnan4p [ 0 ] * rtB . e3hmodlujb [ 1 ] ; rtB . b3yujzbge2 [ 3 ] = rtB .
e3hmodlujb [ 1 ] * rtB . fzbgmnan4p [ 2 ] ; rtB . b3yujzbge2 [ 4 ] = rtB .
fzbgmnan4p [ 0 ] * rtB . e3hmodlujb [ 2 ] ; rtB . b3yujzbge2 [ 5 ] = rtB .
e3hmodlujb [ 0 ] * rtB . fzbgmnan4p [ 1 ] ; rtB . gf5n1xgexk [ 0 ] = rtB .
b3yujzbge2 [ 0 ] - rtB . b3yujzbge2 [ 3 ] ; rtB . gf5n1xgexk [ 1 ] = rtB .
b3yujzbge2 [ 1 ] - rtB . b3yujzbge2 [ 4 ] ; rtB . gf5n1xgexk [ 2 ] = rtB .
b3yujzbge2 [ 2 ] - rtB . b3yujzbge2 [ 5 ] ; czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . gf5n1xgexk [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
gf5n1xgexk [ 0 ] ; nsmqj1yv3h = rtB . gf5n1xgexk [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 = ( rtB . fzqka3pwrj [ yIdx + 3 ] *
jj1qbnwy23_idx_2 + rtB . fzqka3pwrj [ yIdx ] * jj1qbnwy23_idx_3 ) + rtB .
fzqka3pwrj [ yIdx + 6 ] * nsmqj1yv3h ; rtB . e3buvmif1c [ yIdx ] =
jj1qbnwy23_idx_1 ; rtB . j4fldanieq [ yIdx ] = rtP . SimParams . torque_onoff
* jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB . oflmpiaez1 [ yIdx ] ;
czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } if ( ssIsMajorTimeStep
( rtS ) ) { if ( rtDW . eyhr14ycpb != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
eyhr14ycpb = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . eyhr14ycpb = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . cdgxtvlobx =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_l5p2q1mni4 ) ; rtB . lbmtuwi01h = (
rtB . oflmpiaez1 [ 0 ] * rtB . oflmpiaez1 [ 0 ] + rtB . oflmpiaez1 [ 1 ] *
rtB . oflmpiaez1 [ 1 ] ) + rtB . oflmpiaez1 [ 2 ] * rtB . oflmpiaez1 [ 2 ] ;
if ( rtB . lbmtuwi01h < 0.0 ) { rtB . ieyrepmb05 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . lbmtuwi01h ) ) ; } else { rtB . ieyrepmb05 =
muDoubleScalarSqrt ( rtB . lbmtuwi01h ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . pagreemrym = ( rtB . ieyrepmb05 > rtP .
NormalizeVector_maxzero_ihmtmp1vz1 ) ; } if ( rtDW . pagreemrym ) { rtB .
lyjnnzstpa [ 0 ] = rtB . oflmpiaez1 [ 0 ] ; rtB . lyjnnzstpa [ 1 ] = rtB .
oflmpiaez1 [ 1 ] ; rtB . lyjnnzstpa [ 2 ] = rtB . oflmpiaez1 [ 2 ] ; rtB .
lyjnnzstpa [ 3 ] = rtB . ieyrepmb05 ; } else { czut5pi4ll = rtB . oflmpiaez1
[ 0 ] * 0.0 ; rtB . bjhxaihh05 [ 0 ] = czut5pi4ll ; rtB . lyjnnzstpa [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . oflmpiaez1 [ 1 ] * 0.0 ; rtB . bjhxaihh05 [ 1
] = czut5pi4ll ; rtB . lyjnnzstpa [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
oflmpiaez1 [ 2 ] * 0.0 ; rtB . bjhxaihh05 [ 2 ] = czut5pi4ll ; rtB .
lyjnnzstpa [ 2 ] = czut5pi4ll ; rtB . lyjnnzstpa [ 3 ] = rtP .
Constant_Value_cuvlenn5ze ; } jj1qbnwy23_idx_1 = rtB . lyjnnzstpa [ 0 ] / rtB
. lyjnnzstpa [ 3 ] ; rtB . favu2vq25u [ 0 ] = jj1qbnwy23_idx_1 ; czut5pi4ll =
jj1qbnwy23_idx_1 * rtB . oimhqzeh2t [ 0 ] ; jj1qbnwy23_idx_1 = rtB .
lyjnnzstpa [ 1 ] / rtB . lyjnnzstpa [ 3 ] ; rtB . favu2vq25u [ 1 ] =
jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB . oimhqzeh2t [ 1 ] ;
jj1qbnwy23_idx_1 = rtB . lyjnnzstpa [ 2 ] / rtB . lyjnnzstpa [ 3 ] ; rtB .
favu2vq25u [ 2 ] = jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB .
oimhqzeh2t [ 2 ] ; rtB . h1zdr1aani [ 0 ] = czut5pi4ll * rtB . fzbgmnan4p [ 0
] ; rtB . h1zdr1aani [ 1 ] = czut5pi4ll * rtB . fzbgmnan4p [ 1 ] ; rtB .
h1zdr1aani [ 2 ] = czut5pi4ll * rtB . fzbgmnan4p [ 2 ] ; jj1qbnwy23_idx_1 = (
rtB . favu2vq25u [ 0 ] * rtB . fzbgmnan4p [ 0 ] + rtB . favu2vq25u [ 1 ] *
rtB . fzbgmnan4p [ 1 ] ) + rtB . favu2vq25u [ 2 ] * rtB . fzbgmnan4p [ 2 ] ;
rtB . cab5adcq1r = czut5pi4ll * jj1qbnwy23_idx_1 * rtP . Gain_Gain_kzgbyxe2u1
; rtB . m4bi5zptkh [ 0 ] = jj1qbnwy23_idx_1 * rtB . oimhqzeh2t [ 0 ] ; rtB .
m4bi5zptkh [ 1 ] = jj1qbnwy23_idx_1 * rtB . oimhqzeh2t [ 1 ] ; rtB .
m4bi5zptkh [ 2 ] = jj1qbnwy23_idx_1 * rtB . oimhqzeh2t [ 2 ] ; rtB .
avcgwoxwhd = ( ( rtB . oimhqzeh2t [ 0 ] * rtB . fzbgmnan4p [ 0 ] + rtB .
oimhqzeh2t [ 1 ] * rtB . fzbgmnan4p [ 1 ] ) + rtB . oimhqzeh2t [ 2 ] * rtB .
fzbgmnan4p [ 2 ] ) - rtB . cab5adcq1r ; czut5pi4ll = rtB . favu2vq25u [ 0 ] *
rtB . avcgwoxwhd ; rtB . f4nnctaymu [ 0 ] = czut5pi4ll ; czut5pi4ll += rtB .
h1zdr1aani [ 0 ] + rtB . m4bi5zptkh [ 0 ] ; rtB . mcasp4wuog [ 0 ] =
czut5pi4ll ; czut5pi4ll *= rtB . cdgxtvlobx ; rtB . mrdivqgmvm [ 0 ] =
czut5pi4ll ; rtB . g3vpjafupq [ 0 ] = rtP . SimParams . force_onoff *
czut5pi4ll ; czut5pi4ll = rtB . favu2vq25u [ 1 ] * rtB . avcgwoxwhd ; rtB .
f4nnctaymu [ 1 ] = czut5pi4ll ; czut5pi4ll += rtB . h1zdr1aani [ 1 ] + rtB .
m4bi5zptkh [ 1 ] ; rtB . mcasp4wuog [ 1 ] = czut5pi4ll ; czut5pi4ll *= rtB .
cdgxtvlobx ; rtB . mrdivqgmvm [ 1 ] = czut5pi4ll ; rtB . g3vpjafupq [ 1 ] =
rtP . SimParams . force_onoff * czut5pi4ll ; czut5pi4ll = rtB . favu2vq25u [
2 ] * rtB . avcgwoxwhd ; rtB . f4nnctaymu [ 2 ] = czut5pi4ll ; czut5pi4ll +=
rtB . h1zdr1aani [ 2 ] + rtB . m4bi5zptkh [ 2 ] ; rtB . mcasp4wuog [ 2 ] =
czut5pi4ll ; czut5pi4ll *= rtB . cdgxtvlobx ; rtB . mrdivqgmvm [ 2 ] =
czut5pi4ll ; rtB . g3vpjafupq [ 2 ] = rtP . SimParams . force_onoff *
czut5pi4ll ; czut5pi4ll = ( rtB . fzbgmnan4p [ 0 ] * rtB . fzbgmnan4p [ 0 ] +
rtB . fzbgmnan4p [ 1 ] * rtB . fzbgmnan4p [ 1 ] ) + rtB . fzbgmnan4p [ 2 ] *
rtB . fzbgmnan4p [ 2 ] ; if ( czut5pi4ll < 0.0 ) { rtB . btylsma43j = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) ) ; } else { rtB .
btylsma43j = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . kkiofr3le1 = ( rtB . btylsma43j >
rtP . NormalizeVector1_maxzero_byir5pwee4 ) ; } if ( rtDW . kkiofr3le1 ) {
czut5pi4ll = rtB . fzbgmnan4p [ 0 ] ; jj1qbnwy23_idx_1 = rtB . fzbgmnan4p [ 1
] ; jj1qbnwy23_idx_2 = rtB . fzbgmnan4p [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
btylsma43j ; } else { jj1qbnwy23_idx_2 = rtB . fzbgmnan4p [ 0 ] * 0.0 ; rtB .
emllridil0 [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . fzbgmnan4p [ 1 ] * 0.0 ; rtB . emllridil0 [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . fzbgmnan4p [ 2 ] * 0.0 ; rtB . emllridil0 [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_kqizw1abe3 ; } rtB . my424itvki [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . fmqgelgm0i [ 0 ] = rtB . nfyrk0gene [
142 ] ; rtB . my424itvki [ 1 ] = jj1qbnwy23_idx_1 / jj1qbnwy23_idx_3 ; rtB .
fmqgelgm0i [ 1 ] = rtB . nfyrk0gene [ 143 ] ; rtB . my424itvki [ 2 ] =
jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB . fmqgelgm0i [ 2 ] = rtB .
nfyrk0gene [ 144 ] ; czut5pi4ll = rtB . fmqgelgm0i [ 0 ] * rtP . x [ 13 ] ;
rtB . gmivsxqevf [ 0 ] = czut5pi4ll ; rtB . cxwqk2amda [ 0 ] = rtB .
nfyrk0gene [ 133 ] ; cazgyvwuzb [ 0 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll
= rtB . fmqgelgm0i [ 1 ] * rtP . x [ 13 ] ; rtB . gmivsxqevf [ 1 ] =
czut5pi4ll ; rtB . cxwqk2amda [ 1 ] = rtB . nfyrk0gene [ 134 ] ; cazgyvwuzb [
1 ] = czut5pi4ll * czut5pi4ll ; czut5pi4ll = rtB . fmqgelgm0i [ 2 ] * rtP . x
[ 13 ] ; rtB . gmivsxqevf [ 2 ] = czut5pi4ll ; rtB . cxwqk2amda [ 2 ] = rtB .
nfyrk0gene [ 135 ] ; jj1qbnwy23_idx_1 = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ]
) + czut5pi4ll * czut5pi4ll ; if ( jj1qbnwy23_idx_1 < 0.0 ) { rtB .
jgkbpstejl = - muDoubleScalarSqrt ( muDoubleScalarAbs ( jj1qbnwy23_idx_1 ) )
; } else { rtB . jgkbpstejl = muDoubleScalarSqrt ( jj1qbnwy23_idx_1 ) ; } if
( ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . pqnrnlrnbt = ( rtB . jgkbpstejl >
rtP . NormalizeVector_maxzero_fqwaltgljo ) ; } if ( rtDW . pqnrnlrnbt ) {
czut5pi4ll = rtB . gmivsxqevf [ 0 ] ; jj1qbnwy23_idx_1 = rtB . gmivsxqevf [ 1
] ; jj1qbnwy23_idx_2 = rtB . gmivsxqevf [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
jgkbpstejl ; } else { jj1qbnwy23_idx_2 = rtB . gmivsxqevf [ 0 ] * 0.0 ; rtB .
ibott2adsg [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . gmivsxqevf [ 1 ] * 0.0 ; rtB . ibott2adsg [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . gmivsxqevf [ 2 ] * 0.0 ; rtB . ibott2adsg [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_nfejvnig4t ; } rtB . evftq3sk0q [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . evftq3sk0q [ 1 ] = jj1qbnwy23_idx_1 /
jj1qbnwy23_idx_3 ; rtB . evftq3sk0q [ 2 ] = jj1qbnwy23_idx_2 /
jj1qbnwy23_idx_3 ; tmp_g = & rtB . nfyrk0gene [ 145 ] ; for ( yIdx = 0 ; yIdx
< 3 ; yIdx ++ ) { rtB . hnzvyzmoy4 [ 3 * yIdx ] = tmp_g [ yIdx ] ; rtB .
hnzvyzmoy4 [ 3 * yIdx + 1 ] = tmp_g [ yIdx + 3 ] ; rtB . hnzvyzmoy4 [ 3 *
yIdx + 2 ] = tmp_g [ yIdx + 6 ] ; jj1qbnwy23_idx_1 = rtB . cxwqk2amda [ yIdx
] ; cazgyvwuzb [ yIdx ] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } rtB .
gtu1mufolx = ( cazgyvwuzb [ 0 ] + cazgyvwuzb [ 1 ] ) + cazgyvwuzb [ 2 ] ; if
( rtB . gtu1mufolx < 0.0 ) { rtB . ifl2ongx1d = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . gtu1mufolx ) ) ; } else { rtB . ifl2ongx1d =
muDoubleScalarSqrt ( rtB . gtu1mufolx ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . d1ahpk1qqr = ( rtB . ifl2ongx1d > rtP .
NormalizeVector_maxzero_hkbi1qbiut ) ; } if ( rtDW . d1ahpk1qqr ) { rtB .
a3qlbfbxxo [ 0 ] = rtB . cxwqk2amda [ 0 ] ; rtB . a3qlbfbxxo [ 1 ] = rtB .
cxwqk2amda [ 1 ] ; rtB . a3qlbfbxxo [ 2 ] = rtB . cxwqk2amda [ 2 ] ; rtB .
a3qlbfbxxo [ 3 ] = rtB . ifl2ongx1d ; } else { czut5pi4ll = rtB . cxwqk2amda
[ 0 ] * 0.0 ; rtB . eejq0hbv30 [ 0 ] = czut5pi4ll ; rtB . a3qlbfbxxo [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . cxwqk2amda [ 1 ] * 0.0 ; rtB . eejq0hbv30 [ 1
] = czut5pi4ll ; rtB . a3qlbfbxxo [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
cxwqk2amda [ 2 ] * 0.0 ; rtB . eejq0hbv30 [ 2 ] = czut5pi4ll ; rtB .
a3qlbfbxxo [ 2 ] = czut5pi4ll ; rtB . a3qlbfbxxo [ 3 ] = rtP .
Constant_Value_m0tx3vya3h ; } czut5pi4ll = rtB . a3qlbfbxxo [ 0 ] / rtB .
a3qlbfbxxo [ 3 ] ; rtB . glnk1qsnt1 [ 0 ] = czut5pi4ll ; rtB . li0pbutf3d [ 0
] = czut5pi4ll ; czut5pi4ll = rtB . a3qlbfbxxo [ 1 ] / rtB . a3qlbfbxxo [ 3 ]
; rtB . glnk1qsnt1 [ 1 ] = czut5pi4ll ; rtB . li0pbutf3d [ 1 ] = czut5pi4ll ;
czut5pi4ll = rtB . a3qlbfbxxo [ 2 ] / rtB . a3qlbfbxxo [ 3 ] ; rtB .
glnk1qsnt1 [ 2 ] = czut5pi4ll ; rtB . li0pbutf3d [ 2 ] = czut5pi4ll ;
czut5pi4ll = rtB . glnk1qsnt1 [ 0 ] ; jj1qbnwy23_idx_1 = rtB . glnk1qsnt1 [ 1
] ; jj1qbnwy23_idx_2 = rtB . glnk1qsnt1 [ 2 ] ; for ( yIdx = 0 ; yIdx < 3 ;
yIdx ++ ) { jj1qbnwy23_idx_3 = rtB . li0pbutf3d [ yIdx ] ; rtB . nzexd0tfyn [
3 * yIdx ] = czut5pi4ll * jj1qbnwy23_idx_3 ; rtB . nzexd0tfyn [ 3 * yIdx + 1
] = jj1qbnwy23_idx_1 * jj1qbnwy23_idx_3 ; rtB . nzexd0tfyn [ 3 * yIdx + 2 ] =
jj1qbnwy23_idx_2 * jj1qbnwy23_idx_3 ; } for ( yIdx = 0 ; yIdx < 9 ; yIdx ++ )
{ czut5pi4ll = rtP . Gain1_Gain_fzh3lvxzcm * rtB . nzexd0tfyn [ yIdx ] ; rtB
. lrkvtt341w [ yIdx ] = czut5pi4ll ; rtB . kptzapguuy [ yIdx ] = czut5pi4ll -
rtP . Constant1_Value_k2adofcs3i [ yIdx ] ; } czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . bcudxny0c0 [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
bcudxny0c0 [ 0 ] ; nsmqj1yv3h = rtB . bcudxny0c0 [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { rtB . nwco0rdsmq [ yIdx ] = ( rtB . kptzapguuy [ yIdx
+ 3 ] * jj1qbnwy23_idx_2 + rtB . kptzapguuy [ yIdx ] * jj1qbnwy23_idx_3 ) +
rtB . kptzapguuy [ yIdx + 6 ] * nsmqj1yv3h ; jj1qbnwy23_idx_1 = rtB .
cxwqk2amda [ yIdx ] ; czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; }
if ( ssIsMajorTimeStep ( rtS ) ) { if ( rtDW . gmlynrr0kl != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
gmlynrr0kl = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . gmlynrr0kl = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . j2h4dmgm1o =
rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_nu3mnt2521 ) ; rtB . kdzk5onxbh [ 0 ]
= rtB . nwco0rdsmq [ 0 ] * rtB . j2h4dmgm1o ; rtB . kdzk5onxbh [ 1 ] = rtB .
nwco0rdsmq [ 1 ] * rtB . j2h4dmgm1o ; rtB . kdzk5onxbh [ 2 ] = rtB .
nwco0rdsmq [ 2 ] * rtB . j2h4dmgm1o ; rtB . oiy31a51av [ 0 ] = rtB .
gmivsxqevf [ 1 ] * rtB . kdzk5onxbh [ 2 ] ; rtB . oiy31a51av [ 1 ] = rtB .
kdzk5onxbh [ 0 ] * rtB . gmivsxqevf [ 2 ] ; rtB . oiy31a51av [ 2 ] = rtB .
gmivsxqevf [ 0 ] * rtB . kdzk5onxbh [ 1 ] ; rtB . oiy31a51av [ 3 ] = rtB .
kdzk5onxbh [ 1 ] * rtB . gmivsxqevf [ 2 ] ; rtB . oiy31a51av [ 4 ] = rtB .
gmivsxqevf [ 0 ] * rtB . kdzk5onxbh [ 2 ] ; rtB . oiy31a51av [ 5 ] = rtB .
kdzk5onxbh [ 0 ] * rtB . gmivsxqevf [ 1 ] ; rtB . athb42mz2d [ 0 ] = rtB .
oiy31a51av [ 0 ] - rtB . oiy31a51av [ 3 ] ; rtB . athb42mz2d [ 1 ] = rtB .
oiy31a51av [ 1 ] - rtB . oiy31a51av [ 4 ] ; rtB . athb42mz2d [ 2 ] = rtB .
oiy31a51av [ 2 ] - rtB . oiy31a51av [ 5 ] ; czut5pi4ll = 0.0 ;
jj1qbnwy23_idx_2 = rtB . athb42mz2d [ 1 ] ; jj1qbnwy23_idx_3 = rtB .
athb42mz2d [ 0 ] ; nsmqj1yv3h = rtB . athb42mz2d [ 2 ] ; for ( yIdx = 0 ;
yIdx < 3 ; yIdx ++ ) { jj1qbnwy23_idx_1 = ( rtB . hnzvyzmoy4 [ yIdx + 3 ] *
jj1qbnwy23_idx_2 + rtB . hnzvyzmoy4 [ yIdx ] * jj1qbnwy23_idx_3 ) + rtB .
hnzvyzmoy4 [ yIdx + 6 ] * nsmqj1yv3h ; rtB . exqp2uohe5 [ yIdx ] =
jj1qbnwy23_idx_1 ; rtB . o2znllyvoq [ yIdx ] = rtP . SimParams . torque_onoff
* jj1qbnwy23_idx_1 ; jj1qbnwy23_idx_1 = rtB . cxwqk2amda [ yIdx ] ;
czut5pi4ll += jj1qbnwy23_idx_1 * jj1qbnwy23_idx_1 ; } if ( ssIsMajorTimeStep
( rtS ) ) { if ( rtDW . dfembtpulp != 0 ) {
ssSetBlockStateForSolverChangedAtMajorStep ( rtS ) ;
ssSetContTimeOutputInconsistentWithStateAtMajorStep ( rtS ) ; rtDW .
dfembtpulp = 0 ; } jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; }
else if ( czut5pi4ll < 0.0 ) { jj1qbnwy23_idx_1 = - muDoubleScalarSqrt (
muDoubleScalarAbs ( czut5pi4ll ) ) ; rtDW . dfembtpulp = 1 ; } else {
jj1qbnwy23_idx_1 = muDoubleScalarSqrt ( czut5pi4ll ) ; } rtB . bbtvepo4p4 =
3.0 * rtP . Fixed . myu / 12.566370614359172 * muDoubleScalarPower (
jj1qbnwy23_idx_1 , rtP . Constant_Value_h2tkrf03re ) ; rtB . kq5inrtfnq = (
rtB . cxwqk2amda [ 0 ] * rtB . cxwqk2amda [ 0 ] + rtB . cxwqk2amda [ 1 ] *
rtB . cxwqk2amda [ 1 ] ) + rtB . cxwqk2amda [ 2 ] * rtB . cxwqk2amda [ 2 ] ;
if ( rtB . kq5inrtfnq < 0.0 ) { rtB . c00wpmk13d = - muDoubleScalarSqrt (
muDoubleScalarAbs ( rtB . kq5inrtfnq ) ) ; } else { rtB . c00wpmk13d =
muDoubleScalarSqrt ( rtB . kq5inrtfnq ) ; } if ( ssIsModeUpdateTimeStep ( rtS
) ) { rtDW . hw5n3aggbh = ( rtB . c00wpmk13d > rtP .
NormalizeVector_maxzero_dnlmszn2ri ) ; } if ( rtDW . hw5n3aggbh ) { rtB .
huacygyj0o [ 0 ] = rtB . cxwqk2amda [ 0 ] ; rtB . huacygyj0o [ 1 ] = rtB .
cxwqk2amda [ 1 ] ; rtB . huacygyj0o [ 2 ] = rtB . cxwqk2amda [ 2 ] ; rtB .
huacygyj0o [ 3 ] = rtB . c00wpmk13d ; } else { czut5pi4ll = rtB . cxwqk2amda
[ 0 ] * 0.0 ; rtB . joumqp0zpt [ 0 ] = czut5pi4ll ; rtB . huacygyj0o [ 0 ] =
czut5pi4ll ; czut5pi4ll = rtB . cxwqk2amda [ 1 ] * 0.0 ; rtB . joumqp0zpt [ 1
] = czut5pi4ll ; rtB . huacygyj0o [ 1 ] = czut5pi4ll ; czut5pi4ll = rtB .
cxwqk2amda [ 2 ] * 0.0 ; rtB . joumqp0zpt [ 2 ] = czut5pi4ll ; rtB .
huacygyj0o [ 2 ] = czut5pi4ll ; rtB . huacygyj0o [ 3 ] = rtP .
Constant_Value_oc3ljvgxhv ; } jj1qbnwy23_idx_1 = rtB . huacygyj0o [ 0 ] / rtB
. huacygyj0o [ 3 ] ; rtB . gg2oxfw35t [ 0 ] = jj1qbnwy23_idx_1 ; czut5pi4ll =
jj1qbnwy23_idx_1 * rtB . bcudxny0c0 [ 0 ] ; jj1qbnwy23_idx_1 = rtB .
huacygyj0o [ 1 ] / rtB . huacygyj0o [ 3 ] ; rtB . gg2oxfw35t [ 1 ] =
jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB . bcudxny0c0 [ 1 ] ;
jj1qbnwy23_idx_1 = rtB . huacygyj0o [ 2 ] / rtB . huacygyj0o [ 3 ] ; rtB .
gg2oxfw35t [ 2 ] = jj1qbnwy23_idx_1 ; czut5pi4ll += jj1qbnwy23_idx_1 * rtB .
bcudxny0c0 [ 2 ] ; rtB . lperucgciu [ 0 ] = czut5pi4ll * rtB . gmivsxqevf [ 0
] ; rtB . lperucgciu [ 1 ] = czut5pi4ll * rtB . gmivsxqevf [ 1 ] ; rtB .
lperucgciu [ 2 ] = czut5pi4ll * rtB . gmivsxqevf [ 2 ] ; jj1qbnwy23_idx_1 = (
rtB . gg2oxfw35t [ 0 ] * rtB . gmivsxqevf [ 0 ] + rtB . gg2oxfw35t [ 1 ] *
rtB . gmivsxqevf [ 1 ] ) + rtB . gg2oxfw35t [ 2 ] * rtB . gmivsxqevf [ 2 ] ;
rtB . aizglk3yjb = czut5pi4ll * jj1qbnwy23_idx_1 * rtP . Gain_Gain_gnnygtgrhf
; rtB . jqxy4sgz1r [ 0 ] = jj1qbnwy23_idx_1 * rtB . bcudxny0c0 [ 0 ] ; rtB .
jqxy4sgz1r [ 1 ] = jj1qbnwy23_idx_1 * rtB . bcudxny0c0 [ 1 ] ; rtB .
jqxy4sgz1r [ 2 ] = jj1qbnwy23_idx_1 * rtB . bcudxny0c0 [ 2 ] ; rtB .
avjauevoe4 = ( ( rtB . bcudxny0c0 [ 0 ] * rtB . gmivsxqevf [ 0 ] + rtB .
bcudxny0c0 [ 1 ] * rtB . gmivsxqevf [ 1 ] ) + rtB . bcudxny0c0 [ 2 ] * rtB .
gmivsxqevf [ 2 ] ) - rtB . aizglk3yjb ; czut5pi4ll = rtB . gg2oxfw35t [ 0 ] *
rtB . avjauevoe4 ; rtB . azxdea4e00 [ 0 ] = czut5pi4ll ; czut5pi4ll += rtB .
lperucgciu [ 0 ] + rtB . jqxy4sgz1r [ 0 ] ; rtB . afxh4mt3gu [ 0 ] =
czut5pi4ll ; rtB . gdhtjuaiqp [ 0 ] = rtB . bbtvepo4p4 * czut5pi4ll ;
czut5pi4ll = rtB . gg2oxfw35t [ 1 ] * rtB . avjauevoe4 ; rtB . azxdea4e00 [ 1
] = czut5pi4ll ; czut5pi4ll += rtB . lperucgciu [ 1 ] + rtB . jqxy4sgz1r [ 1
] ; rtB . afxh4mt3gu [ 1 ] = czut5pi4ll ; rtB . gdhtjuaiqp [ 1 ] = rtB .
bbtvepo4p4 * czut5pi4ll ; czut5pi4ll = rtB . gg2oxfw35t [ 2 ] * rtB .
avjauevoe4 ; rtB . azxdea4e00 [ 2 ] = czut5pi4ll ; czut5pi4ll += rtB .
lperucgciu [ 2 ] + rtB . jqxy4sgz1r [ 2 ] ; rtB . afxh4mt3gu [ 2 ] =
czut5pi4ll ; rtB . gdhtjuaiqp [ 2 ] = rtB . bbtvepo4p4 * czut5pi4ll ;
czut5pi4ll = ( rtB . gmivsxqevf [ 0 ] * rtB . gmivsxqevf [ 0 ] + rtB .
gmivsxqevf [ 1 ] * rtB . gmivsxqevf [ 1 ] ) + rtB . gmivsxqevf [ 2 ] * rtB .
gmivsxqevf [ 2 ] ; if ( czut5pi4ll < 0.0 ) { rtB . e0gxuicxzo = -
muDoubleScalarSqrt ( muDoubleScalarAbs ( czut5pi4ll ) ) ; } else { rtB .
e0gxuicxzo = muDoubleScalarSqrt ( czut5pi4ll ) ; } if (
ssIsModeUpdateTimeStep ( rtS ) ) { rtDW . kgduzxmrs5 = ( rtB . e0gxuicxzo >
rtP . NormalizeVector1_maxzero_j5bp1hk5bd ) ; } if ( rtDW . kgduzxmrs5 ) {
czut5pi4ll = rtB . gmivsxqevf [ 0 ] ; jj1qbnwy23_idx_1 = rtB . gmivsxqevf [ 1
] ; jj1qbnwy23_idx_2 = rtB . gmivsxqevf [ 2 ] ; jj1qbnwy23_idx_3 = rtB .
e0gxuicxzo ; } else { jj1qbnwy23_idx_2 = rtB . gmivsxqevf [ 0 ] * 0.0 ; rtB .
ovchpiv2ws [ 0 ] = jj1qbnwy23_idx_2 ; czut5pi4ll = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_2 = rtB . gmivsxqevf [ 1 ] * 0.0 ; rtB . ovchpiv2ws [ 1 ] =
jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_1 = jj1qbnwy23_idx_2 ; jj1qbnwy23_idx_2 =
rtB . gmivsxqevf [ 2 ] * 0.0 ; rtB . ovchpiv2ws [ 2 ] = jj1qbnwy23_idx_2 ;
jj1qbnwy23_idx_3 = rtP . Constant_Value_jgubgdhktv ; } rtB . kbih5351xy [ 0 ]
= czut5pi4ll / jj1qbnwy23_idx_3 ; rtB . j2tnneegmr [ 0 ] = rtP . SimParams .
force_onoff * rtB . gdhtjuaiqp [ 0 ] ; rtB . kbih5351xy [ 1 ] =
jj1qbnwy23_idx_1 / jj1qbnwy23_idx_3 ; rtB . j2tnneegmr [ 1 ] = rtP .
SimParams . force_onoff * rtB . gdhtjuaiqp [ 1 ] ; rtB . kbih5351xy [ 2 ] =
jj1qbnwy23_idx_2 / jj1qbnwy23_idx_3 ; rtB . j2tnneegmr [ 2 ] = rtP .
SimParams . force_onoff * rtB . gdhtjuaiqp [ 2 ] ; rtB . c513o5pg5r [ 0 ] =
rtB . j4tq1cqvap [ 0 ] ; rtB . c513o5pg5r [ 1 ] = 0.0 ; rtB . c513o5pg5r [ 2
] = 0.0 ; rtB . c513o5pg5r [ 3 ] = 0.0 ; rtB . pfmi5qvzhf [ 0 ] = rtB .
j4tq1cqvap [ 1 ] ; rtB . pfmi5qvzhf [ 1 ] = 0.0 ; rtB . pfmi5qvzhf [ 2 ] =
0.0 ; rtB . pfmi5qvzhf [ 3 ] = 0.0 ; rtB . gi35bo1nw4 [ 0 ] = rtB .
j4tq1cqvap [ 2 ] ; rtB . gi35bo1nw4 [ 1 ] = 0.0 ; rtB . gi35bo1nw4 [ 2 ] =
0.0 ; rtB . gi35bo1nw4 [ 3 ] = 0.0 ; rtB . iii5xyl5wn [ 0 ] = rtB .
i3fsmpf0l4 [ 0 ] ; rtB . iii5xyl5wn [ 1 ] = 0.0 ; rtB . iii5xyl5wn [ 2 ] =
0.0 ; rtB . iii5xyl5wn [ 3 ] = 0.0 ; rtB . mjqwg1vruy [ 0 ] = rtB .
i3fsmpf0l4 [ 1 ] ; rtB . mjqwg1vruy [ 1 ] = 0.0 ; rtB . mjqwg1vruy [ 2 ] =
0.0 ; rtB . mjqwg1vruy [ 3 ] = 0.0 ; rtB . eo0yhyzcnn [ 0 ] = rtB .
i3fsmpf0l4 [ 2 ] ; rtB . eo0yhyzcnn [ 1 ] = 0.0 ; rtB . eo0yhyzcnn [ 2 ] =
0.0 ; rtB . eo0yhyzcnn [ 3 ] = 0.0 ; rtB . c5dvbg4axm [ 0 ] = rtB .
dakqbtmvna [ 0 ] ; rtB . c5dvbg4axm [ 1 ] = 0.0 ; rtB . c5dvbg4axm [ 2 ] =
0.0 ; rtB . c5dvbg4axm [ 3 ] = 0.0 ; rtB . bxbu5xamaf [ 0 ] = rtB .
dakqbtmvna [ 1 ] ; rtB . bxbu5xamaf [ 1 ] = 0.0 ; rtB . bxbu5xamaf [ 2 ] =
0.0 ; rtB . bxbu5xamaf [ 3 ] = 0.0 ; rtB . foqtyz1zpt [ 0 ] = rtB .
dakqbtmvna [ 2 ] ; rtB . foqtyz1zpt [ 1 ] = 0.0 ; rtB . foqtyz1zpt [ 2 ] =
0.0 ; rtB . foqtyz1zpt [ 3 ] = 0.0 ; rtB . gy3u0julyl [ 0 ] = rtB .
cwwkrda1qf [ 0 ] ; rtB . gy3u0julyl [ 1 ] = 0.0 ; rtB . gy3u0julyl [ 2 ] =
0.0 ; rtB . gy3u0julyl [ 3 ] = 0.0 ; rtB . hytxiurw0s [ 0 ] = rtB .
cwwkrda1qf [ 1 ] ; rtB . hytxiurw0s [ 1 ] = 0.0 ; rtB . hytxiurw0s [ 2 ] =
0.0 ; rtB . hytxiurw0s [ 3 ] = 0.0 ; rtB . nubcraskrw [ 0 ] = rtB .
cwwkrda1qf [ 2 ] ; rtB . nubcraskrw [ 1 ] = 0.0 ; rtB . nubcraskrw [ 2 ] =
0.0 ; rtB . nubcraskrw [ 3 ] = 0.0 ; rtB . lp0cxari03 [ 0 ] = rtB .
dndiwk2b1s [ 0 ] ; rtB . lp0cxari03 [ 1 ] = 0.0 ; rtB . lp0cxari03 [ 2 ] =
0.0 ; rtB . lp0cxari03 [ 3 ] = 0.0 ; rtB . phx5lsta5w [ 0 ] = rtB .
dndiwk2b1s [ 1 ] ; rtB . phx5lsta5w [ 1 ] = 0.0 ; rtB . phx5lsta5w [ 2 ] =
0.0 ; rtB . phx5lsta5w [ 3 ] = 0.0 ; rtB . hoap0zr3eu [ 0 ] = rtB .
dndiwk2b1s [ 2 ] ; rtB . hoap0zr3eu [ 1 ] = 0.0 ; rtB . hoap0zr3eu [ 2 ] =
0.0 ; rtB . hoap0zr3eu [ 3 ] = 0.0 ; rtB . enp5f3s002 [ 0 ] = rtB .
nalsmrc4um [ 0 ] ; rtB . enp5f3s002 [ 1 ] = 0.0 ; rtB . enp5f3s002 [ 2 ] =
0.0 ; rtB . enp5f3s002 [ 3 ] = 0.0 ; rtB . nsl31ekm0f [ 0 ] = rtB .
nalsmrc4um [ 1 ] ; rtB . nsl31ekm0f [ 1 ] = 0.0 ; rtB . nsl31ekm0f [ 2 ] =
0.0 ; rtB . nsl31ekm0f [ 3 ] = 0.0 ; rtB . hrqvwki2jl [ 0 ] = rtB .
nalsmrc4um [ 2 ] ; rtB . hrqvwki2jl [ 1 ] = 0.0 ; rtB . hrqvwki2jl [ 2 ] =
0.0 ; rtB . hrqvwki2jl [ 3 ] = 0.0 ; rtB . dypcby02yn [ 0 ] = rtB .
cbmxpbd2iu [ 0 ] ; rtB . dypcby02yn [ 1 ] = 0.0 ; rtB . dypcby02yn [ 2 ] =
0.0 ; rtB . dypcby02yn [ 3 ] = 0.0 ; rtB . a100glwqkz [ 0 ] = rtB .
cbmxpbd2iu [ 1 ] ; rtB . a100glwqkz [ 1 ] = 0.0 ; rtB . a100glwqkz [ 2 ] =
0.0 ; rtB . a100glwqkz [ 3 ] = 0.0 ; rtB . c3twv1nzdj [ 0 ] = rtB .
cbmxpbd2iu [ 2 ] ; rtB . c3twv1nzdj [ 1 ] = 0.0 ; rtB . c3twv1nzdj [ 2 ] =
0.0 ; rtB . c3twv1nzdj [ 3 ] = 0.0 ; rtB . iurh4h1cla [ 0 ] = rtB .
b0dwpzqyrm [ 0 ] ; rtB . iurh4h1cla [ 1 ] = 0.0 ; rtB . iurh4h1cla [ 2 ] =
0.0 ; rtB . iurh4h1cla [ 3 ] = 0.0 ; rtB . har0a2wrmi [ 0 ] = rtB .
b0dwpzqyrm [ 1 ] ; rtB . har0a2wrmi [ 1 ] = 0.0 ; rtB . har0a2wrmi [ 2 ] =
0.0 ; rtB . har0a2wrmi [ 3 ] = 0.0 ; rtB . ib2ieuf2an [ 0 ] = rtB .
b0dwpzqyrm [ 2 ] ; rtB . ib2ieuf2an [ 1 ] = 0.0 ; rtB . ib2ieuf2an [ 2 ] =
0.0 ; rtB . ib2ieuf2an [ 3 ] = 0.0 ; rtB . ojpovrodvr [ 0 ] = rtB .
l01qgkxm4q [ 0 ] ; rtB . ojpovrodvr [ 1 ] = 0.0 ; rtB . ojpovrodvr [ 2 ] =
0.0 ; rtB . ojpovrodvr [ 3 ] = 0.0 ; rtB . fcgmpow5uz [ 0 ] = rtB .
l01qgkxm4q [ 1 ] ; rtB . fcgmpow5uz [ 1 ] = 0.0 ; rtB . fcgmpow5uz [ 2 ] =
0.0 ; rtB . fcgmpow5uz [ 3 ] = 0.0 ; rtB . jylntbxfpj [ 0 ] = rtB .
l01qgkxm4q [ 2 ] ; rtB . jylntbxfpj [ 1 ] = 0.0 ; rtB . jylntbxfpj [ 2 ] =
0.0 ; rtB . jylntbxfpj [ 3 ] = 0.0 ; rtB . g1yo2p3pxt [ 0 ] = rtB .
lxihdocxga [ 0 ] ; rtB . g1yo2p3pxt [ 1 ] = 0.0 ; rtB . g1yo2p3pxt [ 2 ] =
0.0 ; rtB . g1yo2p3pxt [ 3 ] = 0.0 ; rtB . az4eo3zh5o [ 0 ] = rtB .
lxihdocxga [ 1 ] ; rtB . az4eo3zh5o [ 1 ] = 0.0 ; rtB . az4eo3zh5o [ 2 ] =
0.0 ; rtB . az4eo3zh5o [ 3 ] = 0.0 ; rtB . ch1tbnhyon [ 0 ] = rtB .
lxihdocxga [ 2 ] ; rtB . ch1tbnhyon [ 1 ] = 0.0 ; rtB . ch1tbnhyon [ 2 ] =
0.0 ; rtB . ch1tbnhyon [ 3 ] = 0.0 ; rtB . eoseoj2f3b [ 0 ] = rtB .
g3vpjafupq [ 0 ] ; rtB . eoseoj2f3b [ 1 ] = 0.0 ; rtB . eoseoj2f3b [ 2 ] =
0.0 ; rtB . eoseoj2f3b [ 3 ] = 0.0 ; rtB . aigpe1li2p [ 0 ] = rtB .
g3vpjafupq [ 1 ] ; rtB . aigpe1li2p [ 1 ] = 0.0 ; rtB . aigpe1li2p [ 2 ] =
0.0 ; rtB . aigpe1li2p [ 3 ] = 0.0 ; rtB . pbn2nzk2k3 [ 0 ] = rtB .
g3vpjafupq [ 2 ] ; rtB . pbn2nzk2k3 [ 1 ] = 0.0 ; rtB . pbn2nzk2k3 [ 2 ] =
0.0 ; rtB . pbn2nzk2k3 [ 3 ] = 0.0 ; rtB . hffhwq1imf [ 0 ] = rtB .
j4fldanieq [ 0 ] ; rtB . hffhwq1imf [ 1 ] = 0.0 ; rtB . hffhwq1imf [ 2 ] =
0.0 ; rtB . hffhwq1imf [ 3 ] = 0.0 ; rtB . gbxuu4wtte [ 0 ] = rtB .
j4fldanieq [ 1 ] ; rtB . gbxuu4wtte [ 1 ] = 0.0 ; rtB . gbxuu4wtte [ 2 ] =
0.0 ; rtB . gbxuu4wtte [ 3 ] = 0.0 ; rtB . k5nlme3kns [ 0 ] = rtB .
j4fldanieq [ 2 ] ; rtB . k5nlme3kns [ 1 ] = 0.0 ; rtB . k5nlme3kns [ 2 ] =
0.0 ; rtB . k5nlme3kns [ 3 ] = 0.0 ; rtB . j52cfj2d1w [ 0 ] = rtB .
j2tnneegmr [ 0 ] ; rtB . j52cfj2d1w [ 1 ] = 0.0 ; rtB . j52cfj2d1w [ 2 ] =
0.0 ; rtB . j52cfj2d1w [ 3 ] = 0.0 ; rtB . ftvz2ekipg [ 0 ] = rtB .
j2tnneegmr [ 1 ] ; rtB . ftvz2ekipg [ 1 ] = 0.0 ; rtB . ftvz2ekipg [ 2 ] =
0.0 ; rtB . ftvz2ekipg [ 3 ] = 0.0 ; rtB . awgskjjmsb [ 0 ] = rtB .
j2tnneegmr [ 2 ] ; rtB . awgskjjmsb [ 1 ] = 0.0 ; rtB . awgskjjmsb [ 2 ] =
0.0 ; rtB . awgskjjmsb [ 3 ] = 0.0 ; rtB . fl4rvtkarn [ 0 ] = rtB .
o2znllyvoq [ 0 ] ; rtB . fl4rvtkarn [ 1 ] = 0.0 ; rtB . fl4rvtkarn [ 2 ] =
0.0 ; rtB . fl4rvtkarn [ 3 ] = 0.0 ; rtB . jp22vgjwwm [ 0 ] = rtB .
o2znllyvoq [ 1 ] ; rtB . jp22vgjwwm [ 1 ] = 0.0 ; rtB . jp22vgjwwm [ 2 ] =
0.0 ; rtB . jp22vgjwwm [ 3 ] = 0.0 ; rtB . pzpjv5s5n2 [ 0 ] = rtB .
o2znllyvoq [ 2 ] ; rtB . pzpjv5s5n2 [ 1 ] = 0.0 ; rtB . pzpjv5s5n2 [ 2 ] =
0.0 ; rtB . pzpjv5s5n2 [ 3 ] = 0.0 ; UNUSED_PARAMETER ( tid ) ; } void
MdlOutputsTID1 ( int_T tid ) { rtB . nsmqj1yv3h [ 0 ] = rtP .
MagnetDipoleMoment_Value [ 0 ] ; rtB . mdrxplc5qx [ 0 ] = rtP .
MagnetDipoleMoment_Value_hus3icj2zx [ 0 ] ; rtB . idcqc0uibn [ 0 ] = rtP .
MagnetDipoleMoment_Value_m4h5ot0d2d [ 0 ] ; rtB . lkwtwvkutl [ 0 ] = rtP .
MagnetDipoleMoment_Value_am3fjqdljc [ 0 ] ; rtB . dunjuwexfg [ 0 ] = rtP .
MagnetDipoleMoment_Value_hqiqmgkprz [ 0 ] ; rtB . oimhqzeh2t [ 0 ] = rtP .
MagnetDipoleMoment_Value_clffg1mnnq [ 0 ] ; rtB . bcudxny0c0 [ 0 ] = rtP .
MagnetDipoleMoment_Value_af4j3ouf2p [ 0 ] ; rtB . nsmqj1yv3h [ 1 ] = rtP .
MagnetDipoleMoment_Value [ 1 ] ; rtB . mdrxplc5qx [ 1 ] = rtP .
MagnetDipoleMoment_Value_hus3icj2zx [ 1 ] ; rtB . idcqc0uibn [ 1 ] = rtP .
MagnetDipoleMoment_Value_m4h5ot0d2d [ 1 ] ; rtB . lkwtwvkutl [ 1 ] = rtP .
MagnetDipoleMoment_Value_am3fjqdljc [ 1 ] ; rtB . dunjuwexfg [ 1 ] = rtP .
MagnetDipoleMoment_Value_hqiqmgkprz [ 1 ] ; rtB . oimhqzeh2t [ 1 ] = rtP .
MagnetDipoleMoment_Value_clffg1mnnq [ 1 ] ; rtB . bcudxny0c0 [ 1 ] = rtP .
MagnetDipoleMoment_Value_af4j3ouf2p [ 1 ] ; rtB . nsmqj1yv3h [ 2 ] = rtP .
MagnetDipoleMoment_Value [ 2 ] ; rtB . mdrxplc5qx [ 2 ] = rtP .
MagnetDipoleMoment_Value_hus3icj2zx [ 2 ] ; rtB . idcqc0uibn [ 2 ] = rtP .
MagnetDipoleMoment_Value_m4h5ot0d2d [ 2 ] ; rtB . lkwtwvkutl [ 2 ] = rtP .
MagnetDipoleMoment_Value_am3fjqdljc [ 2 ] ; rtB . dunjuwexfg [ 2 ] = rtP .
MagnetDipoleMoment_Value_hqiqmgkprz [ 2 ] ; rtB . oimhqzeh2t [ 2 ] = rtP .
MagnetDipoleMoment_Value_clffg1mnnq [ 2 ] ; rtB . bcudxny0c0 [ 2 ] = rtP .
MagnetDipoleMoment_Value_af4j3ouf2p [ 2 ] ; UNUSED_PARAMETER ( tid ) ; } void
MdlUpdate ( int_T tid ) { NeslSimulationData * simulationData ;
NeuDiagnosticManager * diagnosticManager ; NeuDiagnosticTree * diagnosticTree
; char * msg ; real_T tmp_p [ 168 ] ; real_T time ; int32_T tmp_i ; int_T
tmp_e [ 43 ] ; boolean_T tmp ; simulationData = ( NeslSimulationData * ) rtDW
. m2wgpmvjvy ; time = ssGetT ( rtS ) ; simulationData -> mData -> mTime . mN
= 1 ; simulationData -> mData -> mTime . mX = & time ; simulationData ->
mData -> mContStates . mN = 18 ; simulationData -> mData -> mContStates . mX
= & rtX . eouclruk0d [ 0 ] ; simulationData -> mData -> mDiscStates . mN = 0
; simulationData -> mData -> mDiscStates . mX = & rtDW . exk4rrasy3 ;
simulationData -> mData -> mModeVector . mN = 0 ; simulationData -> mData ->
mModeVector . mX = & rtDW . h2bg1hotvh ; tmp = ( ssIsMajorTimeStep ( rtS ) &&
ssGetRTWSolverInfo ( rtS ) -> foundContZcEvents ) ; simulationData -> mData
-> mFoundZcEvents = tmp ; simulationData -> mData -> mIsMajorTimeStep =
ssIsMajorTimeStep ( rtS ) ; tmp = ( ssGetMdlInfoPtr ( rtS ) -> mdlFlags .
solverAssertCheck == 1U ) ; simulationData -> mData -> mIsSolverAssertCheck =
tmp ; tmp = ssIsSolverCheckingCIC ( rtS ) ; simulationData -> mData ->
mIsSolverCheckingCIC = tmp ; tmp = ssIsSolverComputingJacobian ( rtS ) ;
simulationData -> mData -> mIsComputingJacobian = tmp ; simulationData ->
mData -> mIsEvaluatingF0 = ( ssGetEvaluatingF0ForJacobian ( rtS ) != 0 ) ;
tmp = ssIsSolverRequestingReset ( rtS ) ; simulationData -> mData ->
mIsSolverRequestingReset = tmp ; simulationData -> mData ->
mIsModeUpdateTimeStep = ssIsModeUpdateTimeStep ( rtS ) ; tmp_e [ 0 ] = 0 ;
tmp_p [ 0 ] = rtB . c513o5pg5r [ 0 ] ; tmp_p [ 1 ] = rtB . c513o5pg5r [ 1 ] ;
tmp_p [ 2 ] = rtB . c513o5pg5r [ 2 ] ; tmp_p [ 3 ] = rtB . c513o5pg5r [ 3 ] ;
tmp_e [ 1 ] = 4 ; tmp_p [ 4 ] = rtB . pfmi5qvzhf [ 0 ] ; tmp_p [ 5 ] = rtB .
pfmi5qvzhf [ 1 ] ; tmp_p [ 6 ] = rtB . pfmi5qvzhf [ 2 ] ; tmp_p [ 7 ] = rtB .
pfmi5qvzhf [ 3 ] ; tmp_e [ 2 ] = 8 ; tmp_p [ 8 ] = rtB . gi35bo1nw4 [ 0 ] ;
tmp_p [ 9 ] = rtB . gi35bo1nw4 [ 1 ] ; tmp_p [ 10 ] = rtB . gi35bo1nw4 [ 2 ]
; tmp_p [ 11 ] = rtB . gi35bo1nw4 [ 3 ] ; tmp_e [ 3 ] = 12 ; tmp_p [ 12 ] =
rtB . iii5xyl5wn [ 0 ] ; tmp_p [ 13 ] = rtB . iii5xyl5wn [ 1 ] ; tmp_p [ 14 ]
= rtB . iii5xyl5wn [ 2 ] ; tmp_p [ 15 ] = rtB . iii5xyl5wn [ 3 ] ; tmp_e [ 4
] = 16 ; tmp_p [ 16 ] = rtB . mjqwg1vruy [ 0 ] ; tmp_p [ 17 ] = rtB .
mjqwg1vruy [ 1 ] ; tmp_p [ 18 ] = rtB . mjqwg1vruy [ 2 ] ; tmp_p [ 19 ] = rtB
. mjqwg1vruy [ 3 ] ; tmp_e [ 5 ] = 20 ; tmp_p [ 20 ] = rtB . eo0yhyzcnn [ 0 ]
; tmp_p [ 21 ] = rtB . eo0yhyzcnn [ 1 ] ; tmp_p [ 22 ] = rtB . eo0yhyzcnn [ 2
] ; tmp_p [ 23 ] = rtB . eo0yhyzcnn [ 3 ] ; tmp_e [ 6 ] = 24 ; tmp_p [ 24 ] =
rtB . c5dvbg4axm [ 0 ] ; tmp_p [ 25 ] = rtB . c5dvbg4axm [ 1 ] ; tmp_p [ 26 ]
= rtB . c5dvbg4axm [ 2 ] ; tmp_p [ 27 ] = rtB . c5dvbg4axm [ 3 ] ; tmp_e [ 7
] = 28 ; tmp_p [ 28 ] = rtB . bxbu5xamaf [ 0 ] ; tmp_p [ 29 ] = rtB .
bxbu5xamaf [ 1 ] ; tmp_p [ 30 ] = rtB . bxbu5xamaf [ 2 ] ; tmp_p [ 31 ] = rtB
. bxbu5xamaf [ 3 ] ; tmp_e [ 8 ] = 32 ; tmp_p [ 32 ] = rtB . foqtyz1zpt [ 0 ]
; tmp_p [ 33 ] = rtB . foqtyz1zpt [ 1 ] ; tmp_p [ 34 ] = rtB . foqtyz1zpt [ 2
] ; tmp_p [ 35 ] = rtB . foqtyz1zpt [ 3 ] ; tmp_e [ 9 ] = 36 ; tmp_p [ 36 ] =
rtB . gy3u0julyl [ 0 ] ; tmp_p [ 37 ] = rtB . gy3u0julyl [ 1 ] ; tmp_p [ 38 ]
= rtB . gy3u0julyl [ 2 ] ; tmp_p [ 39 ] = rtB . gy3u0julyl [ 3 ] ; tmp_e [ 10
] = 40 ; tmp_p [ 40 ] = rtB . hytxiurw0s [ 0 ] ; tmp_p [ 41 ] = rtB .
hytxiurw0s [ 1 ] ; tmp_p [ 42 ] = rtB . hytxiurw0s [ 2 ] ; tmp_p [ 43 ] = rtB
. hytxiurw0s [ 3 ] ; tmp_e [ 11 ] = 44 ; tmp_p [ 44 ] = rtB . nubcraskrw [ 0
] ; tmp_p [ 45 ] = rtB . nubcraskrw [ 1 ] ; tmp_p [ 46 ] = rtB . nubcraskrw [
2 ] ; tmp_p [ 47 ] = rtB . nubcraskrw [ 3 ] ; tmp_e [ 12 ] = 48 ; tmp_p [ 48
] = rtB . lp0cxari03 [ 0 ] ; tmp_p [ 49 ] = rtB . lp0cxari03 [ 1 ] ; tmp_p [
50 ] = rtB . lp0cxari03 [ 2 ] ; tmp_p [ 51 ] = rtB . lp0cxari03 [ 3 ] ; tmp_e
[ 13 ] = 52 ; tmp_p [ 52 ] = rtB . phx5lsta5w [ 0 ] ; tmp_p [ 53 ] = rtB .
phx5lsta5w [ 1 ] ; tmp_p [ 54 ] = rtB . phx5lsta5w [ 2 ] ; tmp_p [ 55 ] = rtB
. phx5lsta5w [ 3 ] ; tmp_e [ 14 ] = 56 ; tmp_p [ 56 ] = rtB . hoap0zr3eu [ 0
] ; tmp_p [ 57 ] = rtB . hoap0zr3eu [ 1 ] ; tmp_p [ 58 ] = rtB . hoap0zr3eu [
2 ] ; tmp_p [ 59 ] = rtB . hoap0zr3eu [ 3 ] ; tmp_e [ 15 ] = 60 ; tmp_p [ 60
] = rtB . enp5f3s002 [ 0 ] ; tmp_p [ 61 ] = rtB . enp5f3s002 [ 1 ] ; tmp_p [
62 ] = rtB . enp5f3s002 [ 2 ] ; tmp_p [ 63 ] = rtB . enp5f3s002 [ 3 ] ; tmp_e
[ 16 ] = 64 ; tmp_p [ 64 ] = rtB . nsl31ekm0f [ 0 ] ; tmp_p [ 65 ] = rtB .
nsl31ekm0f [ 1 ] ; tmp_p [ 66 ] = rtB . nsl31ekm0f [ 2 ] ; tmp_p [ 67 ] = rtB
. nsl31ekm0f [ 3 ] ; tmp_e [ 17 ] = 68 ; tmp_p [ 68 ] = rtB . hrqvwki2jl [ 0
] ; tmp_p [ 69 ] = rtB . hrqvwki2jl [ 1 ] ; tmp_p [ 70 ] = rtB . hrqvwki2jl [
2 ] ; tmp_p [ 71 ] = rtB . hrqvwki2jl [ 3 ] ; tmp_e [ 18 ] = 72 ; tmp_p [ 72
] = rtB . dypcby02yn [ 0 ] ; tmp_p [ 73 ] = rtB . dypcby02yn [ 1 ] ; tmp_p [
74 ] = rtB . dypcby02yn [ 2 ] ; tmp_p [ 75 ] = rtB . dypcby02yn [ 3 ] ; tmp_e
[ 19 ] = 76 ; tmp_p [ 76 ] = rtB . a100glwqkz [ 0 ] ; tmp_p [ 77 ] = rtB .
a100glwqkz [ 1 ] ; tmp_p [ 78 ] = rtB . a100glwqkz [ 2 ] ; tmp_p [ 79 ] = rtB
. a100glwqkz [ 3 ] ; tmp_e [ 20 ] = 80 ; tmp_p [ 80 ] = rtB . c3twv1nzdj [ 0
] ; tmp_p [ 81 ] = rtB . c3twv1nzdj [ 1 ] ; tmp_p [ 82 ] = rtB . c3twv1nzdj [
2 ] ; tmp_p [ 83 ] = rtB . c3twv1nzdj [ 3 ] ; tmp_e [ 21 ] = 84 ; tmp_p [ 84
] = rtB . iurh4h1cla [ 0 ] ; tmp_p [ 85 ] = rtB . iurh4h1cla [ 1 ] ; tmp_p [
86 ] = rtB . iurh4h1cla [ 2 ] ; tmp_p [ 87 ] = rtB . iurh4h1cla [ 3 ] ; tmp_e
[ 22 ] = 88 ; tmp_p [ 88 ] = rtB . har0a2wrmi [ 0 ] ; tmp_p [ 89 ] = rtB .
har0a2wrmi [ 1 ] ; tmp_p [ 90 ] = rtB . har0a2wrmi [ 2 ] ; tmp_p [ 91 ] = rtB
. har0a2wrmi [ 3 ] ; tmp_e [ 23 ] = 92 ; tmp_p [ 92 ] = rtB . ib2ieuf2an [ 0
] ; tmp_p [ 93 ] = rtB . ib2ieuf2an [ 1 ] ; tmp_p [ 94 ] = rtB . ib2ieuf2an [
2 ] ; tmp_p [ 95 ] = rtB . ib2ieuf2an [ 3 ] ; tmp_e [ 24 ] = 96 ; tmp_p [ 96
] = rtB . ojpovrodvr [ 0 ] ; tmp_p [ 97 ] = rtB . ojpovrodvr [ 1 ] ; tmp_p [
98 ] = rtB . ojpovrodvr [ 2 ] ; tmp_p [ 99 ] = rtB . ojpovrodvr [ 3 ] ; tmp_e
[ 25 ] = 100 ; tmp_p [ 100 ] = rtB . fcgmpow5uz [ 0 ] ; tmp_p [ 101 ] = rtB .
fcgmpow5uz [ 1 ] ; tmp_p [ 102 ] = rtB . fcgmpow5uz [ 2 ] ; tmp_p [ 103 ] =
rtB . fcgmpow5uz [ 3 ] ; tmp_e [ 26 ] = 104 ; tmp_p [ 104 ] = rtB .
jylntbxfpj [ 0 ] ; tmp_p [ 105 ] = rtB . jylntbxfpj [ 1 ] ; tmp_p [ 106 ] =
rtB . jylntbxfpj [ 2 ] ; tmp_p [ 107 ] = rtB . jylntbxfpj [ 3 ] ; tmp_e [ 27
] = 108 ; tmp_p [ 108 ] = rtB . g1yo2p3pxt [ 0 ] ; tmp_p [ 109 ] = rtB .
g1yo2p3pxt [ 1 ] ; tmp_p [ 110 ] = rtB . g1yo2p3pxt [ 2 ] ; tmp_p [ 111 ] =
rtB . g1yo2p3pxt [ 3 ] ; tmp_e [ 28 ] = 112 ; tmp_p [ 112 ] = rtB .
az4eo3zh5o [ 0 ] ; tmp_p [ 113 ] = rtB . az4eo3zh5o [ 1 ] ; tmp_p [ 114 ] =
rtB . az4eo3zh5o [ 2 ] ; tmp_p [ 115 ] = rtB . az4eo3zh5o [ 3 ] ; tmp_e [ 29
] = 116 ; tmp_p [ 116 ] = rtB . ch1tbnhyon [ 0 ] ; tmp_p [ 117 ] = rtB .
ch1tbnhyon [ 1 ] ; tmp_p [ 118 ] = rtB . ch1tbnhyon [ 2 ] ; tmp_p [ 119 ] =
rtB . ch1tbnhyon [ 3 ] ; tmp_e [ 30 ] = 120 ; tmp_p [ 120 ] = rtB .
eoseoj2f3b [ 0 ] ; tmp_p [ 121 ] = rtB . eoseoj2f3b [ 1 ] ; tmp_p [ 122 ] =
rtB . eoseoj2f3b [ 2 ] ; tmp_p [ 123 ] = rtB . eoseoj2f3b [ 3 ] ; tmp_e [ 31
] = 124 ; tmp_p [ 124 ] = rtB . aigpe1li2p [ 0 ] ; tmp_p [ 125 ] = rtB .
aigpe1li2p [ 1 ] ; tmp_p [ 126 ] = rtB . aigpe1li2p [ 2 ] ; tmp_p [ 127 ] =
rtB . aigpe1li2p [ 3 ] ; tmp_e [ 32 ] = 128 ; tmp_p [ 128 ] = rtB .
pbn2nzk2k3 [ 0 ] ; tmp_p [ 129 ] = rtB . pbn2nzk2k3 [ 1 ] ; tmp_p [ 130 ] =
rtB . pbn2nzk2k3 [ 2 ] ; tmp_p [ 131 ] = rtB . pbn2nzk2k3 [ 3 ] ; tmp_e [ 33
] = 132 ; tmp_p [ 132 ] = rtB . hffhwq1imf [ 0 ] ; tmp_p [ 133 ] = rtB .
hffhwq1imf [ 1 ] ; tmp_p [ 134 ] = rtB . hffhwq1imf [ 2 ] ; tmp_p [ 135 ] =
rtB . hffhwq1imf [ 3 ] ; tmp_e [ 34 ] = 136 ; tmp_p [ 136 ] = rtB .
gbxuu4wtte [ 0 ] ; tmp_p [ 137 ] = rtB . gbxuu4wtte [ 1 ] ; tmp_p [ 138 ] =
rtB . gbxuu4wtte [ 2 ] ; tmp_p [ 139 ] = rtB . gbxuu4wtte [ 3 ] ; tmp_e [ 35
] = 140 ; tmp_p [ 140 ] = rtB . k5nlme3kns [ 0 ] ; tmp_p [ 141 ] = rtB .
k5nlme3kns [ 1 ] ; tmp_p [ 142 ] = rtB . k5nlme3kns [ 2 ] ; tmp_p [ 143 ] =
rtB . k5nlme3kns [ 3 ] ; tmp_e [ 36 ] = 144 ; tmp_p [ 144 ] = rtB .
j52cfj2d1w [ 0 ] ; tmp_p [ 145 ] = rtB . j52cfj2d1w [ 1 ] ; tmp_p [ 146 ] =
rtB . j52cfj2d1w [ 2 ] ; tmp_p [ 147 ] = rtB . j52cfj2d1w [ 3 ] ; tmp_e [ 37
] = 148 ; tmp_p [ 148 ] = rtB . ftvz2ekipg [ 0 ] ; tmp_p [ 149 ] = rtB .
ftvz2ekipg [ 1 ] ; tmp_p [ 150 ] = rtB . ftvz2ekipg [ 2 ] ; tmp_p [ 151 ] =
rtB . ftvz2ekipg [ 3 ] ; tmp_e [ 38 ] = 152 ; tmp_p [ 152 ] = rtB .
awgskjjmsb [ 0 ] ; tmp_p [ 153 ] = rtB . awgskjjmsb [ 1 ] ; tmp_p [ 154 ] =
rtB . awgskjjmsb [ 2 ] ; tmp_p [ 155 ] = rtB . awgskjjmsb [ 3 ] ; tmp_e [ 39
] = 156 ; tmp_p [ 156 ] = rtB . fl4rvtkarn [ 0 ] ; tmp_p [ 157 ] = rtB .
fl4rvtkarn [ 1 ] ; tmp_p [ 158 ] = rtB . fl4rvtkarn [ 2 ] ; tmp_p [ 159 ] =
rtB . fl4rvtkarn [ 3 ] ; tmp_e [ 40 ] = 160 ; tmp_p [ 160 ] = rtB .
jp22vgjwwm [ 0 ] ; tmp_p [ 161 ] = rtB . jp22vgjwwm [ 1 ] ; tmp_p [ 162 ] =
rtB . jp22vgjwwm [ 2 ] ; tmp_p [ 163 ] = rtB . jp22vgjwwm [ 3 ] ; tmp_e [ 41
] = 164 ; tmp_p [ 164 ] = rtB . pzpjv5s5n2 [ 0 ] ; tmp_p [ 165 ] = rtB .
pzpjv5s5n2 [ 1 ] ; tmp_p [ 166 ] = rtB . pzpjv5s5n2 [ 2 ] ; tmp_p [ 167 ] =
rtB . pzpjv5s5n2 [ 3 ] ; tmp_e [ 42 ] = 168 ; simulationData -> mData ->
mInputValues . mN = 168 ; simulationData -> mData -> mInputValues . mX = &
tmp_p [ 0 ] ; simulationData -> mData -> mInputOffsets . mN = 43 ;
simulationData -> mData -> mInputOffsets . mX = & tmp_e [ 0 ] ;
diagnosticManager = ( NeuDiagnosticManager * ) rtDW . oqiguspxfn ;
diagnosticTree = neu_diagnostic_manager_get_initial_tree ( diagnosticManager
) ; tmp_i = ne_simulator_method ( ( NeslSimulator * ) rtDW . ctycytv1sy ,
NESL_SIM_UPDATE , simulationData , diagnosticManager ) ; if ( tmp_i != 0 ) {
tmp = error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if ( tmp ) { msg =
rtw_diagnostics_msg ( diagnosticTree ) ; ssSetErrorStatus ( rtS , msg ) ; } }
UNUSED_PARAMETER ( tid ) ; } void MdlUpdateTID1 ( int_T tid ) {
UNUSED_PARAMETER ( tid ) ; } void MdlDerivatives ( void ) {
NeslSimulationData * simulationData ; NeuDiagnosticManager *
diagnosticManager ; NeuDiagnosticTree * diagnosticTree ; XDot * _rtXdot ;
char * msg ; real_T tmp_p [ 168 ] ; real_T time ; int32_T tmp_i ; int_T tmp_e
[ 43 ] ; boolean_T tmp ; _rtXdot = ( ( XDot * ) ssGetdX ( rtS ) ) ;
simulationData = ( NeslSimulationData * ) rtDW . m2wgpmvjvy ; time = ssGetT (
rtS ) ; simulationData -> mData -> mTime . mN = 1 ; simulationData -> mData
-> mTime . mX = & time ; simulationData -> mData -> mContStates . mN = 18 ;
simulationData -> mData -> mContStates . mX = & rtX . eouclruk0d [ 0 ] ;
simulationData -> mData -> mDiscStates . mN = 0 ; simulationData -> mData ->
mDiscStates . mX = & rtDW . exk4rrasy3 ; simulationData -> mData ->
mModeVector . mN = 0 ; simulationData -> mData -> mModeVector . mX = & rtDW .
h2bg1hotvh ; tmp = ( ssIsMajorTimeStep ( rtS ) && ssGetRTWSolverInfo ( rtS )
-> foundContZcEvents ) ; simulationData -> mData -> mFoundZcEvents = tmp ;
simulationData -> mData -> mIsMajorTimeStep = ssIsMajorTimeStep ( rtS ) ; tmp
= ( ssGetMdlInfoPtr ( rtS ) -> mdlFlags . solverAssertCheck == 1U ) ;
simulationData -> mData -> mIsSolverAssertCheck = tmp ; tmp =
ssIsSolverCheckingCIC ( rtS ) ; simulationData -> mData ->
mIsSolverCheckingCIC = tmp ; tmp = ssIsSolverComputingJacobian ( rtS ) ;
simulationData -> mData -> mIsComputingJacobian = tmp ; simulationData ->
mData -> mIsEvaluatingF0 = ( ssGetEvaluatingF0ForJacobian ( rtS ) != 0 ) ;
tmp = ssIsSolverRequestingReset ( rtS ) ; simulationData -> mData ->
mIsSolverRequestingReset = tmp ; simulationData -> mData ->
mIsModeUpdateTimeStep = ssIsModeUpdateTimeStep ( rtS ) ; tmp_e [ 0 ] = 0 ;
tmp_p [ 0 ] = rtB . c513o5pg5r [ 0 ] ; tmp_p [ 1 ] = rtB . c513o5pg5r [ 1 ] ;
tmp_p [ 2 ] = rtB . c513o5pg5r [ 2 ] ; tmp_p [ 3 ] = rtB . c513o5pg5r [ 3 ] ;
tmp_e [ 1 ] = 4 ; tmp_p [ 4 ] = rtB . pfmi5qvzhf [ 0 ] ; tmp_p [ 5 ] = rtB .
pfmi5qvzhf [ 1 ] ; tmp_p [ 6 ] = rtB . pfmi5qvzhf [ 2 ] ; tmp_p [ 7 ] = rtB .
pfmi5qvzhf [ 3 ] ; tmp_e [ 2 ] = 8 ; tmp_p [ 8 ] = rtB . gi35bo1nw4 [ 0 ] ;
tmp_p [ 9 ] = rtB . gi35bo1nw4 [ 1 ] ; tmp_p [ 10 ] = rtB . gi35bo1nw4 [ 2 ]
; tmp_p [ 11 ] = rtB . gi35bo1nw4 [ 3 ] ; tmp_e [ 3 ] = 12 ; tmp_p [ 12 ] =
rtB . iii5xyl5wn [ 0 ] ; tmp_p [ 13 ] = rtB . iii5xyl5wn [ 1 ] ; tmp_p [ 14 ]
= rtB . iii5xyl5wn [ 2 ] ; tmp_p [ 15 ] = rtB . iii5xyl5wn [ 3 ] ; tmp_e [ 4
] = 16 ; tmp_p [ 16 ] = rtB . mjqwg1vruy [ 0 ] ; tmp_p [ 17 ] = rtB .
mjqwg1vruy [ 1 ] ; tmp_p [ 18 ] = rtB . mjqwg1vruy [ 2 ] ; tmp_p [ 19 ] = rtB
. mjqwg1vruy [ 3 ] ; tmp_e [ 5 ] = 20 ; tmp_p [ 20 ] = rtB . eo0yhyzcnn [ 0 ]
; tmp_p [ 21 ] = rtB . eo0yhyzcnn [ 1 ] ; tmp_p [ 22 ] = rtB . eo0yhyzcnn [ 2
] ; tmp_p [ 23 ] = rtB . eo0yhyzcnn [ 3 ] ; tmp_e [ 6 ] = 24 ; tmp_p [ 24 ] =
rtB . c5dvbg4axm [ 0 ] ; tmp_p [ 25 ] = rtB . c5dvbg4axm [ 1 ] ; tmp_p [ 26 ]
= rtB . c5dvbg4axm [ 2 ] ; tmp_p [ 27 ] = rtB . c5dvbg4axm [ 3 ] ; tmp_e [ 7
] = 28 ; tmp_p [ 28 ] = rtB . bxbu5xamaf [ 0 ] ; tmp_p [ 29 ] = rtB .
bxbu5xamaf [ 1 ] ; tmp_p [ 30 ] = rtB . bxbu5xamaf [ 2 ] ; tmp_p [ 31 ] = rtB
. bxbu5xamaf [ 3 ] ; tmp_e [ 8 ] = 32 ; tmp_p [ 32 ] = rtB . foqtyz1zpt [ 0 ]
; tmp_p [ 33 ] = rtB . foqtyz1zpt [ 1 ] ; tmp_p [ 34 ] = rtB . foqtyz1zpt [ 2
] ; tmp_p [ 35 ] = rtB . foqtyz1zpt [ 3 ] ; tmp_e [ 9 ] = 36 ; tmp_p [ 36 ] =
rtB . gy3u0julyl [ 0 ] ; tmp_p [ 37 ] = rtB . gy3u0julyl [ 1 ] ; tmp_p [ 38 ]
= rtB . gy3u0julyl [ 2 ] ; tmp_p [ 39 ] = rtB . gy3u0julyl [ 3 ] ; tmp_e [ 10
] = 40 ; tmp_p [ 40 ] = rtB . hytxiurw0s [ 0 ] ; tmp_p [ 41 ] = rtB .
hytxiurw0s [ 1 ] ; tmp_p [ 42 ] = rtB . hytxiurw0s [ 2 ] ; tmp_p [ 43 ] = rtB
. hytxiurw0s [ 3 ] ; tmp_e [ 11 ] = 44 ; tmp_p [ 44 ] = rtB . nubcraskrw [ 0
] ; tmp_p [ 45 ] = rtB . nubcraskrw [ 1 ] ; tmp_p [ 46 ] = rtB . nubcraskrw [
2 ] ; tmp_p [ 47 ] = rtB . nubcraskrw [ 3 ] ; tmp_e [ 12 ] = 48 ; tmp_p [ 48
] = rtB . lp0cxari03 [ 0 ] ; tmp_p [ 49 ] = rtB . lp0cxari03 [ 1 ] ; tmp_p [
50 ] = rtB . lp0cxari03 [ 2 ] ; tmp_p [ 51 ] = rtB . lp0cxari03 [ 3 ] ; tmp_e
[ 13 ] = 52 ; tmp_p [ 52 ] = rtB . phx5lsta5w [ 0 ] ; tmp_p [ 53 ] = rtB .
phx5lsta5w [ 1 ] ; tmp_p [ 54 ] = rtB . phx5lsta5w [ 2 ] ; tmp_p [ 55 ] = rtB
. phx5lsta5w [ 3 ] ; tmp_e [ 14 ] = 56 ; tmp_p [ 56 ] = rtB . hoap0zr3eu [ 0
] ; tmp_p [ 57 ] = rtB . hoap0zr3eu [ 1 ] ; tmp_p [ 58 ] = rtB . hoap0zr3eu [
2 ] ; tmp_p [ 59 ] = rtB . hoap0zr3eu [ 3 ] ; tmp_e [ 15 ] = 60 ; tmp_p [ 60
] = rtB . enp5f3s002 [ 0 ] ; tmp_p [ 61 ] = rtB . enp5f3s002 [ 1 ] ; tmp_p [
62 ] = rtB . enp5f3s002 [ 2 ] ; tmp_p [ 63 ] = rtB . enp5f3s002 [ 3 ] ; tmp_e
[ 16 ] = 64 ; tmp_p [ 64 ] = rtB . nsl31ekm0f [ 0 ] ; tmp_p [ 65 ] = rtB .
nsl31ekm0f [ 1 ] ; tmp_p [ 66 ] = rtB . nsl31ekm0f [ 2 ] ; tmp_p [ 67 ] = rtB
. nsl31ekm0f [ 3 ] ; tmp_e [ 17 ] = 68 ; tmp_p [ 68 ] = rtB . hrqvwki2jl [ 0
] ; tmp_p [ 69 ] = rtB . hrqvwki2jl [ 1 ] ; tmp_p [ 70 ] = rtB . hrqvwki2jl [
2 ] ; tmp_p [ 71 ] = rtB . hrqvwki2jl [ 3 ] ; tmp_e [ 18 ] = 72 ; tmp_p [ 72
] = rtB . dypcby02yn [ 0 ] ; tmp_p [ 73 ] = rtB . dypcby02yn [ 1 ] ; tmp_p [
74 ] = rtB . dypcby02yn [ 2 ] ; tmp_p [ 75 ] = rtB . dypcby02yn [ 3 ] ; tmp_e
[ 19 ] = 76 ; tmp_p [ 76 ] = rtB . a100glwqkz [ 0 ] ; tmp_p [ 77 ] = rtB .
a100glwqkz [ 1 ] ; tmp_p [ 78 ] = rtB . a100glwqkz [ 2 ] ; tmp_p [ 79 ] = rtB
. a100glwqkz [ 3 ] ; tmp_e [ 20 ] = 80 ; tmp_p [ 80 ] = rtB . c3twv1nzdj [ 0
] ; tmp_p [ 81 ] = rtB . c3twv1nzdj [ 1 ] ; tmp_p [ 82 ] = rtB . c3twv1nzdj [
2 ] ; tmp_p [ 83 ] = rtB . c3twv1nzdj [ 3 ] ; tmp_e [ 21 ] = 84 ; tmp_p [ 84
] = rtB . iurh4h1cla [ 0 ] ; tmp_p [ 85 ] = rtB . iurh4h1cla [ 1 ] ; tmp_p [
86 ] = rtB . iurh4h1cla [ 2 ] ; tmp_p [ 87 ] = rtB . iurh4h1cla [ 3 ] ; tmp_e
[ 22 ] = 88 ; tmp_p [ 88 ] = rtB . har0a2wrmi [ 0 ] ; tmp_p [ 89 ] = rtB .
har0a2wrmi [ 1 ] ; tmp_p [ 90 ] = rtB . har0a2wrmi [ 2 ] ; tmp_p [ 91 ] = rtB
. har0a2wrmi [ 3 ] ; tmp_e [ 23 ] = 92 ; tmp_p [ 92 ] = rtB . ib2ieuf2an [ 0
] ; tmp_p [ 93 ] = rtB . ib2ieuf2an [ 1 ] ; tmp_p [ 94 ] = rtB . ib2ieuf2an [
2 ] ; tmp_p [ 95 ] = rtB . ib2ieuf2an [ 3 ] ; tmp_e [ 24 ] = 96 ; tmp_p [ 96
] = rtB . ojpovrodvr [ 0 ] ; tmp_p [ 97 ] = rtB . ojpovrodvr [ 1 ] ; tmp_p [
98 ] = rtB . ojpovrodvr [ 2 ] ; tmp_p [ 99 ] = rtB . ojpovrodvr [ 3 ] ; tmp_e
[ 25 ] = 100 ; tmp_p [ 100 ] = rtB . fcgmpow5uz [ 0 ] ; tmp_p [ 101 ] = rtB .
fcgmpow5uz [ 1 ] ; tmp_p [ 102 ] = rtB . fcgmpow5uz [ 2 ] ; tmp_p [ 103 ] =
rtB . fcgmpow5uz [ 3 ] ; tmp_e [ 26 ] = 104 ; tmp_p [ 104 ] = rtB .
jylntbxfpj [ 0 ] ; tmp_p [ 105 ] = rtB . jylntbxfpj [ 1 ] ; tmp_p [ 106 ] =
rtB . jylntbxfpj [ 2 ] ; tmp_p [ 107 ] = rtB . jylntbxfpj [ 3 ] ; tmp_e [ 27
] = 108 ; tmp_p [ 108 ] = rtB . g1yo2p3pxt [ 0 ] ; tmp_p [ 109 ] = rtB .
g1yo2p3pxt [ 1 ] ; tmp_p [ 110 ] = rtB . g1yo2p3pxt [ 2 ] ; tmp_p [ 111 ] =
rtB . g1yo2p3pxt [ 3 ] ; tmp_e [ 28 ] = 112 ; tmp_p [ 112 ] = rtB .
az4eo3zh5o [ 0 ] ; tmp_p [ 113 ] = rtB . az4eo3zh5o [ 1 ] ; tmp_p [ 114 ] =
rtB . az4eo3zh5o [ 2 ] ; tmp_p [ 115 ] = rtB . az4eo3zh5o [ 3 ] ; tmp_e [ 29
] = 116 ; tmp_p [ 116 ] = rtB . ch1tbnhyon [ 0 ] ; tmp_p [ 117 ] = rtB .
ch1tbnhyon [ 1 ] ; tmp_p [ 118 ] = rtB . ch1tbnhyon [ 2 ] ; tmp_p [ 119 ] =
rtB . ch1tbnhyon [ 3 ] ; tmp_e [ 30 ] = 120 ; tmp_p [ 120 ] = rtB .
eoseoj2f3b [ 0 ] ; tmp_p [ 121 ] = rtB . eoseoj2f3b [ 1 ] ; tmp_p [ 122 ] =
rtB . eoseoj2f3b [ 2 ] ; tmp_p [ 123 ] = rtB . eoseoj2f3b [ 3 ] ; tmp_e [ 31
] = 124 ; tmp_p [ 124 ] = rtB . aigpe1li2p [ 0 ] ; tmp_p [ 125 ] = rtB .
aigpe1li2p [ 1 ] ; tmp_p [ 126 ] = rtB . aigpe1li2p [ 2 ] ; tmp_p [ 127 ] =
rtB . aigpe1li2p [ 3 ] ; tmp_e [ 32 ] = 128 ; tmp_p [ 128 ] = rtB .
pbn2nzk2k3 [ 0 ] ; tmp_p [ 129 ] = rtB . pbn2nzk2k3 [ 1 ] ; tmp_p [ 130 ] =
rtB . pbn2nzk2k3 [ 2 ] ; tmp_p [ 131 ] = rtB . pbn2nzk2k3 [ 3 ] ; tmp_e [ 33
] = 132 ; tmp_p [ 132 ] = rtB . hffhwq1imf [ 0 ] ; tmp_p [ 133 ] = rtB .
hffhwq1imf [ 1 ] ; tmp_p [ 134 ] = rtB . hffhwq1imf [ 2 ] ; tmp_p [ 135 ] =
rtB . hffhwq1imf [ 3 ] ; tmp_e [ 34 ] = 136 ; tmp_p [ 136 ] = rtB .
gbxuu4wtte [ 0 ] ; tmp_p [ 137 ] = rtB . gbxuu4wtte [ 1 ] ; tmp_p [ 138 ] =
rtB . gbxuu4wtte [ 2 ] ; tmp_p [ 139 ] = rtB . gbxuu4wtte [ 3 ] ; tmp_e [ 35
] = 140 ; tmp_p [ 140 ] = rtB . k5nlme3kns [ 0 ] ; tmp_p [ 141 ] = rtB .
k5nlme3kns [ 1 ] ; tmp_p [ 142 ] = rtB . k5nlme3kns [ 2 ] ; tmp_p [ 143 ] =
rtB . k5nlme3kns [ 3 ] ; tmp_e [ 36 ] = 144 ; tmp_p [ 144 ] = rtB .
j52cfj2d1w [ 0 ] ; tmp_p [ 145 ] = rtB . j52cfj2d1w [ 1 ] ; tmp_p [ 146 ] =
rtB . j52cfj2d1w [ 2 ] ; tmp_p [ 147 ] = rtB . j52cfj2d1w [ 3 ] ; tmp_e [ 37
] = 148 ; tmp_p [ 148 ] = rtB . ftvz2ekipg [ 0 ] ; tmp_p [ 149 ] = rtB .
ftvz2ekipg [ 1 ] ; tmp_p [ 150 ] = rtB . ftvz2ekipg [ 2 ] ; tmp_p [ 151 ] =
rtB . ftvz2ekipg [ 3 ] ; tmp_e [ 38 ] = 152 ; tmp_p [ 152 ] = rtB .
awgskjjmsb [ 0 ] ; tmp_p [ 153 ] = rtB . awgskjjmsb [ 1 ] ; tmp_p [ 154 ] =
rtB . awgskjjmsb [ 2 ] ; tmp_p [ 155 ] = rtB . awgskjjmsb [ 3 ] ; tmp_e [ 39
] = 156 ; tmp_p [ 156 ] = rtB . fl4rvtkarn [ 0 ] ; tmp_p [ 157 ] = rtB .
fl4rvtkarn [ 1 ] ; tmp_p [ 158 ] = rtB . fl4rvtkarn [ 2 ] ; tmp_p [ 159 ] =
rtB . fl4rvtkarn [ 3 ] ; tmp_e [ 40 ] = 160 ; tmp_p [ 160 ] = rtB .
jp22vgjwwm [ 0 ] ; tmp_p [ 161 ] = rtB . jp22vgjwwm [ 1 ] ; tmp_p [ 162 ] =
rtB . jp22vgjwwm [ 2 ] ; tmp_p [ 163 ] = rtB . jp22vgjwwm [ 3 ] ; tmp_e [ 41
] = 164 ; tmp_p [ 164 ] = rtB . pzpjv5s5n2 [ 0 ] ; tmp_p [ 165 ] = rtB .
pzpjv5s5n2 [ 1 ] ; tmp_p [ 166 ] = rtB . pzpjv5s5n2 [ 2 ] ; tmp_p [ 167 ] =
rtB . pzpjv5s5n2 [ 3 ] ; tmp_e [ 42 ] = 168 ; simulationData -> mData ->
mInputValues . mN = 168 ; simulationData -> mData -> mInputValues . mX = &
tmp_p [ 0 ] ; simulationData -> mData -> mInputOffsets . mN = 43 ;
simulationData -> mData -> mInputOffsets . mX = & tmp_e [ 0 ] ;
simulationData -> mData -> mDx . mN = 18 ; simulationData -> mData -> mDx .
mX = & _rtXdot -> eouclruk0d [ 0 ] ; diagnosticManager = (
NeuDiagnosticManager * ) rtDW . oqiguspxfn ; diagnosticTree =
neu_diagnostic_manager_get_initial_tree ( diagnosticManager ) ; tmp_i =
ne_simulator_method ( ( NeslSimulator * ) rtDW . ctycytv1sy ,
NESL_SIM_DERIVATIVES , simulationData , diagnosticManager ) ; if ( tmp_i != 0
) { tmp = error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if ( tmp ) {
msg = rtw_diagnostics_msg ( diagnosticTree ) ; ssSetErrorStatus ( rtS , msg )
; } } } void MdlProjection ( void ) { NeslSimulationData * simulationData ;
NeuDiagnosticManager * diagnosticManager ; NeuDiagnosticTree * diagnosticTree
; char * msg ; real_T tmp_p [ 168 ] ; real_T time ; int32_T tmp_i ; int_T
tmp_e [ 43 ] ; boolean_T tmp ; simulationData = ( NeslSimulationData * ) rtDW
. m2wgpmvjvy ; time = ssGetT ( rtS ) ; simulationData -> mData -> mTime . mN
= 1 ; simulationData -> mData -> mTime . mX = & time ; simulationData ->
mData -> mContStates . mN = 18 ; simulationData -> mData -> mContStates . mX
= & rtX . eouclruk0d [ 0 ] ; simulationData -> mData -> mDiscStates . mN = 0
; simulationData -> mData -> mDiscStates . mX = & rtDW . exk4rrasy3 ;
simulationData -> mData -> mModeVector . mN = 0 ; simulationData -> mData ->
mModeVector . mX = & rtDW . h2bg1hotvh ; tmp = ( ssIsMajorTimeStep ( rtS ) &&
ssGetRTWSolverInfo ( rtS ) -> foundContZcEvents ) ; simulationData -> mData
-> mFoundZcEvents = tmp ; simulationData -> mData -> mIsMajorTimeStep =
ssIsMajorTimeStep ( rtS ) ; tmp = ( ssGetMdlInfoPtr ( rtS ) -> mdlFlags .
solverAssertCheck == 1U ) ; simulationData -> mData -> mIsSolverAssertCheck =
tmp ; tmp = ssIsSolverCheckingCIC ( rtS ) ; simulationData -> mData ->
mIsSolverCheckingCIC = tmp ; tmp = ssIsSolverComputingJacobian ( rtS ) ;
simulationData -> mData -> mIsComputingJacobian = tmp ; simulationData ->
mData -> mIsEvaluatingF0 = ( ssGetEvaluatingF0ForJacobian ( rtS ) != 0 ) ;
tmp = ssIsSolverRequestingReset ( rtS ) ; simulationData -> mData ->
mIsSolverRequestingReset = tmp ; simulationData -> mData ->
mIsModeUpdateTimeStep = ssIsModeUpdateTimeStep ( rtS ) ; tmp_e [ 0 ] = 0 ;
tmp_p [ 0 ] = rtB . c513o5pg5r [ 0 ] ; tmp_p [ 1 ] = rtB . c513o5pg5r [ 1 ] ;
tmp_p [ 2 ] = rtB . c513o5pg5r [ 2 ] ; tmp_p [ 3 ] = rtB . c513o5pg5r [ 3 ] ;
tmp_e [ 1 ] = 4 ; tmp_p [ 4 ] = rtB . pfmi5qvzhf [ 0 ] ; tmp_p [ 5 ] = rtB .
pfmi5qvzhf [ 1 ] ; tmp_p [ 6 ] = rtB . pfmi5qvzhf [ 2 ] ; tmp_p [ 7 ] = rtB .
pfmi5qvzhf [ 3 ] ; tmp_e [ 2 ] = 8 ; tmp_p [ 8 ] = rtB . gi35bo1nw4 [ 0 ] ;
tmp_p [ 9 ] = rtB . gi35bo1nw4 [ 1 ] ; tmp_p [ 10 ] = rtB . gi35bo1nw4 [ 2 ]
; tmp_p [ 11 ] = rtB . gi35bo1nw4 [ 3 ] ; tmp_e [ 3 ] = 12 ; tmp_p [ 12 ] =
rtB . iii5xyl5wn [ 0 ] ; tmp_p [ 13 ] = rtB . iii5xyl5wn [ 1 ] ; tmp_p [ 14 ]
= rtB . iii5xyl5wn [ 2 ] ; tmp_p [ 15 ] = rtB . iii5xyl5wn [ 3 ] ; tmp_e [ 4
] = 16 ; tmp_p [ 16 ] = rtB . mjqwg1vruy [ 0 ] ; tmp_p [ 17 ] = rtB .
mjqwg1vruy [ 1 ] ; tmp_p [ 18 ] = rtB . mjqwg1vruy [ 2 ] ; tmp_p [ 19 ] = rtB
. mjqwg1vruy [ 3 ] ; tmp_e [ 5 ] = 20 ; tmp_p [ 20 ] = rtB . eo0yhyzcnn [ 0 ]
; tmp_p [ 21 ] = rtB . eo0yhyzcnn [ 1 ] ; tmp_p [ 22 ] = rtB . eo0yhyzcnn [ 2
] ; tmp_p [ 23 ] = rtB . eo0yhyzcnn [ 3 ] ; tmp_e [ 6 ] = 24 ; tmp_p [ 24 ] =
rtB . c5dvbg4axm [ 0 ] ; tmp_p [ 25 ] = rtB . c5dvbg4axm [ 1 ] ; tmp_p [ 26 ]
= rtB . c5dvbg4axm [ 2 ] ; tmp_p [ 27 ] = rtB . c5dvbg4axm [ 3 ] ; tmp_e [ 7
] = 28 ; tmp_p [ 28 ] = rtB . bxbu5xamaf [ 0 ] ; tmp_p [ 29 ] = rtB .
bxbu5xamaf [ 1 ] ; tmp_p [ 30 ] = rtB . bxbu5xamaf [ 2 ] ; tmp_p [ 31 ] = rtB
. bxbu5xamaf [ 3 ] ; tmp_e [ 8 ] = 32 ; tmp_p [ 32 ] = rtB . foqtyz1zpt [ 0 ]
; tmp_p [ 33 ] = rtB . foqtyz1zpt [ 1 ] ; tmp_p [ 34 ] = rtB . foqtyz1zpt [ 2
] ; tmp_p [ 35 ] = rtB . foqtyz1zpt [ 3 ] ; tmp_e [ 9 ] = 36 ; tmp_p [ 36 ] =
rtB . gy3u0julyl [ 0 ] ; tmp_p [ 37 ] = rtB . gy3u0julyl [ 1 ] ; tmp_p [ 38 ]
= rtB . gy3u0julyl [ 2 ] ; tmp_p [ 39 ] = rtB . gy3u0julyl [ 3 ] ; tmp_e [ 10
] = 40 ; tmp_p [ 40 ] = rtB . hytxiurw0s [ 0 ] ; tmp_p [ 41 ] = rtB .
hytxiurw0s [ 1 ] ; tmp_p [ 42 ] = rtB . hytxiurw0s [ 2 ] ; tmp_p [ 43 ] = rtB
. hytxiurw0s [ 3 ] ; tmp_e [ 11 ] = 44 ; tmp_p [ 44 ] = rtB . nubcraskrw [ 0
] ; tmp_p [ 45 ] = rtB . nubcraskrw [ 1 ] ; tmp_p [ 46 ] = rtB . nubcraskrw [
2 ] ; tmp_p [ 47 ] = rtB . nubcraskrw [ 3 ] ; tmp_e [ 12 ] = 48 ; tmp_p [ 48
] = rtB . lp0cxari03 [ 0 ] ; tmp_p [ 49 ] = rtB . lp0cxari03 [ 1 ] ; tmp_p [
50 ] = rtB . lp0cxari03 [ 2 ] ; tmp_p [ 51 ] = rtB . lp0cxari03 [ 3 ] ; tmp_e
[ 13 ] = 52 ; tmp_p [ 52 ] = rtB . phx5lsta5w [ 0 ] ; tmp_p [ 53 ] = rtB .
phx5lsta5w [ 1 ] ; tmp_p [ 54 ] = rtB . phx5lsta5w [ 2 ] ; tmp_p [ 55 ] = rtB
. phx5lsta5w [ 3 ] ; tmp_e [ 14 ] = 56 ; tmp_p [ 56 ] = rtB . hoap0zr3eu [ 0
] ; tmp_p [ 57 ] = rtB . hoap0zr3eu [ 1 ] ; tmp_p [ 58 ] = rtB . hoap0zr3eu [
2 ] ; tmp_p [ 59 ] = rtB . hoap0zr3eu [ 3 ] ; tmp_e [ 15 ] = 60 ; tmp_p [ 60
] = rtB . enp5f3s002 [ 0 ] ; tmp_p [ 61 ] = rtB . enp5f3s002 [ 1 ] ; tmp_p [
62 ] = rtB . enp5f3s002 [ 2 ] ; tmp_p [ 63 ] = rtB . enp5f3s002 [ 3 ] ; tmp_e
[ 16 ] = 64 ; tmp_p [ 64 ] = rtB . nsl31ekm0f [ 0 ] ; tmp_p [ 65 ] = rtB .
nsl31ekm0f [ 1 ] ; tmp_p [ 66 ] = rtB . nsl31ekm0f [ 2 ] ; tmp_p [ 67 ] = rtB
. nsl31ekm0f [ 3 ] ; tmp_e [ 17 ] = 68 ; tmp_p [ 68 ] = rtB . hrqvwki2jl [ 0
] ; tmp_p [ 69 ] = rtB . hrqvwki2jl [ 1 ] ; tmp_p [ 70 ] = rtB . hrqvwki2jl [
2 ] ; tmp_p [ 71 ] = rtB . hrqvwki2jl [ 3 ] ; tmp_e [ 18 ] = 72 ; tmp_p [ 72
] = rtB . dypcby02yn [ 0 ] ; tmp_p [ 73 ] = rtB . dypcby02yn [ 1 ] ; tmp_p [
74 ] = rtB . dypcby02yn [ 2 ] ; tmp_p [ 75 ] = rtB . dypcby02yn [ 3 ] ; tmp_e
[ 19 ] = 76 ; tmp_p [ 76 ] = rtB . a100glwqkz [ 0 ] ; tmp_p [ 77 ] = rtB .
a100glwqkz [ 1 ] ; tmp_p [ 78 ] = rtB . a100glwqkz [ 2 ] ; tmp_p [ 79 ] = rtB
. a100glwqkz [ 3 ] ; tmp_e [ 20 ] = 80 ; tmp_p [ 80 ] = rtB . c3twv1nzdj [ 0
] ; tmp_p [ 81 ] = rtB . c3twv1nzdj [ 1 ] ; tmp_p [ 82 ] = rtB . c3twv1nzdj [
2 ] ; tmp_p [ 83 ] = rtB . c3twv1nzdj [ 3 ] ; tmp_e [ 21 ] = 84 ; tmp_p [ 84
] = rtB . iurh4h1cla [ 0 ] ; tmp_p [ 85 ] = rtB . iurh4h1cla [ 1 ] ; tmp_p [
86 ] = rtB . iurh4h1cla [ 2 ] ; tmp_p [ 87 ] = rtB . iurh4h1cla [ 3 ] ; tmp_e
[ 22 ] = 88 ; tmp_p [ 88 ] = rtB . har0a2wrmi [ 0 ] ; tmp_p [ 89 ] = rtB .
har0a2wrmi [ 1 ] ; tmp_p [ 90 ] = rtB . har0a2wrmi [ 2 ] ; tmp_p [ 91 ] = rtB
. har0a2wrmi [ 3 ] ; tmp_e [ 23 ] = 92 ; tmp_p [ 92 ] = rtB . ib2ieuf2an [ 0
] ; tmp_p [ 93 ] = rtB . ib2ieuf2an [ 1 ] ; tmp_p [ 94 ] = rtB . ib2ieuf2an [
2 ] ; tmp_p [ 95 ] = rtB . ib2ieuf2an [ 3 ] ; tmp_e [ 24 ] = 96 ; tmp_p [ 96
] = rtB . ojpovrodvr [ 0 ] ; tmp_p [ 97 ] = rtB . ojpovrodvr [ 1 ] ; tmp_p [
98 ] = rtB . ojpovrodvr [ 2 ] ; tmp_p [ 99 ] = rtB . ojpovrodvr [ 3 ] ; tmp_e
[ 25 ] = 100 ; tmp_p [ 100 ] = rtB . fcgmpow5uz [ 0 ] ; tmp_p [ 101 ] = rtB .
fcgmpow5uz [ 1 ] ; tmp_p [ 102 ] = rtB . fcgmpow5uz [ 2 ] ; tmp_p [ 103 ] =
rtB . fcgmpow5uz [ 3 ] ; tmp_e [ 26 ] = 104 ; tmp_p [ 104 ] = rtB .
jylntbxfpj [ 0 ] ; tmp_p [ 105 ] = rtB . jylntbxfpj [ 1 ] ; tmp_p [ 106 ] =
rtB . jylntbxfpj [ 2 ] ; tmp_p [ 107 ] = rtB . jylntbxfpj [ 3 ] ; tmp_e [ 27
] = 108 ; tmp_p [ 108 ] = rtB . g1yo2p3pxt [ 0 ] ; tmp_p [ 109 ] = rtB .
g1yo2p3pxt [ 1 ] ; tmp_p [ 110 ] = rtB . g1yo2p3pxt [ 2 ] ; tmp_p [ 111 ] =
rtB . g1yo2p3pxt [ 3 ] ; tmp_e [ 28 ] = 112 ; tmp_p [ 112 ] = rtB .
az4eo3zh5o [ 0 ] ; tmp_p [ 113 ] = rtB . az4eo3zh5o [ 1 ] ; tmp_p [ 114 ] =
rtB . az4eo3zh5o [ 2 ] ; tmp_p [ 115 ] = rtB . az4eo3zh5o [ 3 ] ; tmp_e [ 29
] = 116 ; tmp_p [ 116 ] = rtB . ch1tbnhyon [ 0 ] ; tmp_p [ 117 ] = rtB .
ch1tbnhyon [ 1 ] ; tmp_p [ 118 ] = rtB . ch1tbnhyon [ 2 ] ; tmp_p [ 119 ] =
rtB . ch1tbnhyon [ 3 ] ; tmp_e [ 30 ] = 120 ; tmp_p [ 120 ] = rtB .
eoseoj2f3b [ 0 ] ; tmp_p [ 121 ] = rtB . eoseoj2f3b [ 1 ] ; tmp_p [ 122 ] =
rtB . eoseoj2f3b [ 2 ] ; tmp_p [ 123 ] = rtB . eoseoj2f3b [ 3 ] ; tmp_e [ 31
] = 124 ; tmp_p [ 124 ] = rtB . aigpe1li2p [ 0 ] ; tmp_p [ 125 ] = rtB .
aigpe1li2p [ 1 ] ; tmp_p [ 126 ] = rtB . aigpe1li2p [ 2 ] ; tmp_p [ 127 ] =
rtB . aigpe1li2p [ 3 ] ; tmp_e [ 32 ] = 128 ; tmp_p [ 128 ] = rtB .
pbn2nzk2k3 [ 0 ] ; tmp_p [ 129 ] = rtB . pbn2nzk2k3 [ 1 ] ; tmp_p [ 130 ] =
rtB . pbn2nzk2k3 [ 2 ] ; tmp_p [ 131 ] = rtB . pbn2nzk2k3 [ 3 ] ; tmp_e [ 33
] = 132 ; tmp_p [ 132 ] = rtB . hffhwq1imf [ 0 ] ; tmp_p [ 133 ] = rtB .
hffhwq1imf [ 1 ] ; tmp_p [ 134 ] = rtB . hffhwq1imf [ 2 ] ; tmp_p [ 135 ] =
rtB . hffhwq1imf [ 3 ] ; tmp_e [ 34 ] = 136 ; tmp_p [ 136 ] = rtB .
gbxuu4wtte [ 0 ] ; tmp_p [ 137 ] = rtB . gbxuu4wtte [ 1 ] ; tmp_p [ 138 ] =
rtB . gbxuu4wtte [ 2 ] ; tmp_p [ 139 ] = rtB . gbxuu4wtte [ 3 ] ; tmp_e [ 35
] = 140 ; tmp_p [ 140 ] = rtB . k5nlme3kns [ 0 ] ; tmp_p [ 141 ] = rtB .
k5nlme3kns [ 1 ] ; tmp_p [ 142 ] = rtB . k5nlme3kns [ 2 ] ; tmp_p [ 143 ] =
rtB . k5nlme3kns [ 3 ] ; tmp_e [ 36 ] = 144 ; tmp_p [ 144 ] = rtB .
j52cfj2d1w [ 0 ] ; tmp_p [ 145 ] = rtB . j52cfj2d1w [ 1 ] ; tmp_p [ 146 ] =
rtB . j52cfj2d1w [ 2 ] ; tmp_p [ 147 ] = rtB . j52cfj2d1w [ 3 ] ; tmp_e [ 37
] = 148 ; tmp_p [ 148 ] = rtB . ftvz2ekipg [ 0 ] ; tmp_p [ 149 ] = rtB .
ftvz2ekipg [ 1 ] ; tmp_p [ 150 ] = rtB . ftvz2ekipg [ 2 ] ; tmp_p [ 151 ] =
rtB . ftvz2ekipg [ 3 ] ; tmp_e [ 38 ] = 152 ; tmp_p [ 152 ] = rtB .
awgskjjmsb [ 0 ] ; tmp_p [ 153 ] = rtB . awgskjjmsb [ 1 ] ; tmp_p [ 154 ] =
rtB . awgskjjmsb [ 2 ] ; tmp_p [ 155 ] = rtB . awgskjjmsb [ 3 ] ; tmp_e [ 39
] = 156 ; tmp_p [ 156 ] = rtB . fl4rvtkarn [ 0 ] ; tmp_p [ 157 ] = rtB .
fl4rvtkarn [ 1 ] ; tmp_p [ 158 ] = rtB . fl4rvtkarn [ 2 ] ; tmp_p [ 159 ] =
rtB . fl4rvtkarn [ 3 ] ; tmp_e [ 40 ] = 160 ; tmp_p [ 160 ] = rtB .
jp22vgjwwm [ 0 ] ; tmp_p [ 161 ] = rtB . jp22vgjwwm [ 1 ] ; tmp_p [ 162 ] =
rtB . jp22vgjwwm [ 2 ] ; tmp_p [ 163 ] = rtB . jp22vgjwwm [ 3 ] ; tmp_e [ 41
] = 164 ; tmp_p [ 164 ] = rtB . pzpjv5s5n2 [ 0 ] ; tmp_p [ 165 ] = rtB .
pzpjv5s5n2 [ 1 ] ; tmp_p [ 166 ] = rtB . pzpjv5s5n2 [ 2 ] ; tmp_p [ 167 ] =
rtB . pzpjv5s5n2 [ 3 ] ; tmp_e [ 42 ] = 168 ; simulationData -> mData ->
mInputValues . mN = 168 ; simulationData -> mData -> mInputValues . mX = &
tmp_p [ 0 ] ; simulationData -> mData -> mInputOffsets . mN = 43 ;
simulationData -> mData -> mInputOffsets . mX = & tmp_e [ 0 ] ;
diagnosticManager = ( NeuDiagnosticManager * ) rtDW . oqiguspxfn ;
diagnosticTree = neu_diagnostic_manager_get_initial_tree ( diagnosticManager
) ; tmp_i = ne_simulator_method ( ( NeslSimulator * ) rtDW . ctycytv1sy ,
NESL_SIM_PROJECTION , simulationData , diagnosticManager ) ; if ( tmp_i != 0
) { tmp = error_buffer_is_empty ( ssGetErrorStatus ( rtS ) ) ; if ( tmp ) {
msg = rtw_diagnostics_msg ( diagnosticTree ) ; ssSetErrorStatus ( rtS , msg )
; } } } void MdlZeroCrossings ( void ) { ZCV * _rtZCSV ; _rtZCSV = ( ( ZCV *
) ssGetSolverZcSignalVector ( rtS ) ) ; _rtZCSV -> gxjau55nwb = rtB .
oywwaak55z - rtP . NormalizeVector_maxzero ; _rtZCSV -> n4sn2x5fn0 = rtB .
i1wwep0quu - rtP . NormalizeVector_maxzero_ab1bjfgz5c ; _rtZCSV -> nxyrsmdlot
= rtB . j0nzmiwajs - rtP . NormalizeVector_maxzero_jpxndxsk21 ; _rtZCSV ->
jvkjfqj1ug = rtB . d4sl1kvckp - rtP . NormalizeVector1_maxzero ; _rtZCSV ->
lge40yyb1l = rtB . jknqpcyts5 - rtP . NormalizeVector_maxzero_hach1vuaw5 ;
_rtZCSV -> bhqttm0cw1 = rtB . eded4se0e2 - rtP .
NormalizeVector_maxzero_nzh4nf0nvr ; _rtZCSV -> lphrnpyiuf = rtB . ltesgv5ex5
- rtP . NormalizeVector_maxzero_pjti1ch2qi ; _rtZCSV -> js0gaa5lhi = rtB .
bkxjvl2wwl - rtP . NormalizeVector1_maxzero_hfyr12og40 ; _rtZCSV ->
eldmkhybny = rtB . n2cw123fsj - rtP . NormalizeVector_maxzero_jvlgwkux2b ;
_rtZCSV -> j13ih0l24o = rtB . lvldddhyjb - rtP .
NormalizeVector_maxzero_k3zfgywmre ; _rtZCSV -> c0cdw0j0oq = rtB . ij2rjc05nz
- rtP . NormalizeVector_maxzero_jceokergvn ; _rtZCSV -> hii3txsdrd = rtB .
idwieepnv4 - rtP . NormalizeVector1_maxzero_edxcmuydlr ; _rtZCSV ->
ohmdhysjjv = rtB . fh04srz5g0 - rtP . NormalizeVector_maxzero_levzkkbbsj ;
_rtZCSV -> bqtyhf5qky = rtB . hlnjhscgik - rtP .
NormalizeVector_maxzero_kt51ehrppa ; _rtZCSV -> iylzqyo3it = rtB . my44vuapi0
- rtP . NormalizeVector_maxzero_hbfzo5qsh2 ; _rtZCSV -> hxfby2z4m1 = rtB .
jqzltcnumt - rtP . NormalizeVector1_maxzero_os1qknjcxt ; _rtZCSV ->
gbmmjctguf = rtB . d4it13dar5 - rtP . NormalizeVector_maxzero_anlse3xi4h ;
_rtZCSV -> bkoti100ma = rtB . b4b10et3s3 - rtP .
NormalizeVector_maxzero_dfin15pekf ; _rtZCSV -> ek2vlosxs0 = rtB . ntt11gc034
- rtP . NormalizeVector_maxzero_nfiirapz00 ; _rtZCSV -> chg25xty1o = rtB .
mp054zeo2w - rtP . NormalizeVector1_maxzero_pbowfqqdx5 ; _rtZCSV ->
hrhezrd0m5 = rtB . efmfmwthqv - rtP . NormalizeVector_maxzero_ap3h53mgxx ;
_rtZCSV -> gk4jezh4kc = rtB . cjrpkdpszw - rtP .
NormalizeVector_maxzero_kytveykcqy ; _rtZCSV -> mccj2kywvs = rtB . ieyrepmb05
- rtP . NormalizeVector_maxzero_ihmtmp1vz1 ; _rtZCSV -> njgojdbb1c = rtB .
btylsma43j - rtP . NormalizeVector1_maxzero_byir5pwee4 ; _rtZCSV ->
ipeaoktcgx = rtB . jgkbpstejl - rtP . NormalizeVector_maxzero_fqwaltgljo ;
_rtZCSV -> ijm0bfdiwn = rtB . ifl2ongx1d - rtP .
NormalizeVector_maxzero_hkbi1qbiut ; _rtZCSV -> osjkfrzf0x = rtB . c00wpmk13d
- rtP . NormalizeVector_maxzero_dnlmszn2ri ; _rtZCSV -> myfa0e3bis = rtB .
e0gxuicxzo - rtP . NormalizeVector1_maxzero_j5bp1hk5bd ; } void MdlTerminate
( void ) { neu_destroy_diagnostic_manager ( ( NeuDiagnosticManager * ) rtDW .
oqiguspxfn ) ; nesl_destroy_simulation_data ( ( NeslSimulationData * ) rtDW .
m2wgpmvjvy ) ; nesl_erase_simulator (
"MagneticBasket_Simscape_Optimizer/Solver Configuration_1" ) ;
nesl_destroy_registry ( ) ; neu_destroy_diagnostic_manager ( (
NeuDiagnosticManager * ) rtDW . dcyq5ws3qs ) ; nesl_destroy_simulation_data (
( NeslSimulationData * ) rtDW . fqoysoyi51 ) ; nesl_erase_simulator (
"MagneticBasket_Simscape_Optimizer/Solver Configuration_1" ) ;
nesl_destroy_registry ( ) ; { if ( rtDW . l4r0r31ubh . AQHandles ) {
sdiTerminateStreaming ( & rtDW . l4r0r31ubh . AQHandles ) ; } } } static void
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( mxArray * destArray
, mwIndex i , int j , const void * srcData , size_t numBytes ) ; static void
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( mxArray * destArray
, mwIndex i , int j , const void * srcData , size_t numBytes ) { mxArray *
newArray = mxCreateUninitNumericMatrix ( ( size_t ) 1 , numBytes ,
mxUINT8_CLASS , mxREAL ) ; memcpy ( ( uint8_T * ) mxGetData ( newArray ) , (
const uint8_T * ) srcData , numBytes ) ; mxSetFieldByNumber ( destArray , i ,
j , newArray ) ; } static void
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( void * destData
, const mxArray * srcArray , mwIndex i , int j , size_t numBytes ) ; static
void mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( void *
destData , const mxArray * srcArray , mwIndex i , int j , size_t numBytes ) {
memcpy ( ( uint8_T * ) destData , ( const uint8_T * ) mxGetData (
mxGetFieldByNumber ( srcArray , i , j ) ) , numBytes ) ; } static void
mr_MagneticBasket_Simscape_Optimizer_cacheBitFieldToMxArray ( mxArray *
destArray , mwIndex i , int j , uint_T bitVal ) ; static void
mr_MagneticBasket_Simscape_Optimizer_cacheBitFieldToMxArray ( mxArray *
destArray , mwIndex i , int j , uint_T bitVal ) { mxSetFieldByNumber (
destArray , i , j , mxCreateDoubleScalar ( ( real_T ) bitVal ) ) ; } static
uint_T mr_MagneticBasket_Simscape_Optimizer_extractBitFieldFromMxArray (
const mxArray * srcArray , mwIndex i , int j , uint_T numBits ) ; static
uint_T mr_MagneticBasket_Simscape_Optimizer_extractBitFieldFromMxArray (
const mxArray * srcArray , mwIndex i , int j , uint_T numBits ) { const
uint_T varVal = ( uint_T ) mxGetScalar ( mxGetFieldByNumber ( srcArray , i ,
j ) ) ; return varVal & ( ( 1u << numBits ) - 1u ) ; } static void
mr_MagneticBasket_Simscape_Optimizer_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) ; static void
mr_MagneticBasket_Simscape_Optimizer_cacheDataToMxArrayWithOffset ( mxArray *
destArray , mwIndex i , int j , mwIndex offset , const void * srcData ,
size_t numBytes ) { uint8_T * varData = ( uint8_T * ) mxGetData (
mxGetFieldByNumber ( destArray , i , j ) ) ; memcpy ( ( uint8_T * ) & varData
[ offset * numBytes ] , ( const uint8_T * ) srcData , numBytes ) ; } static
void mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArrayWithOffset (
void * destData , const mxArray * srcArray , mwIndex i , int j , mwIndex
offset , size_t numBytes ) ; static void
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArrayWithOffset ( void
* destData , const mxArray * srcArray , mwIndex i , int j , mwIndex offset ,
size_t numBytes ) { const uint8_T * varData = ( const uint8_T * ) mxGetData (
mxGetFieldByNumber ( srcArray , i , j ) ) ; memcpy ( ( uint8_T * ) destData ,
( const uint8_T * ) & varData [ offset * numBytes ] , numBytes ) ; } static
void mr_MagneticBasket_Simscape_Optimizer_cacheBitFieldToCellArrayWithOffset
( mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal
) ; static void
mr_MagneticBasket_Simscape_Optimizer_cacheBitFieldToCellArrayWithOffset (
mxArray * destArray , mwIndex i , int j , mwIndex offset , uint_T fieldVal )
{ mxSetCell ( mxGetFieldByNumber ( destArray , i , j ) , offset ,
mxCreateDoubleScalar ( ( real_T ) fieldVal ) ) ; } static uint_T
mr_MagneticBasket_Simscape_Optimizer_extractBitFieldFromCellArrayWithOffset (
const mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T
numBits ) ; static uint_T
mr_MagneticBasket_Simscape_Optimizer_extractBitFieldFromCellArrayWithOffset (
const mxArray * srcArray , mwIndex i , int j , mwIndex offset , uint_T
numBits ) { const uint_T fieldVal = ( uint_T ) mxGetScalar ( mxGetCell (
mxGetFieldByNumber ( srcArray , i , j ) , offset ) ) ; return fieldVal & ( (
1u << numBits ) - 1u ) ; } mxArray *
mr_MagneticBasket_Simscape_Optimizer_GetDWork ( ) { static const char_T *
ssDWFieldNames [ 3 ] = { "rtB" , "rtDW" , "NULL_PrevZCX" , } ; mxArray * ssDW
= mxCreateStructMatrix ( 1 , 1 , 3 , ssDWFieldNames ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( ssDW , 0 , 0 , (
const void * ) & ( rtB ) , sizeof ( rtB ) ) ; { static const char_T *
rtdwDataFieldNames [ 118 ] = { "rtDW.ik1trwr3gq" , "rtDW.m0rr22lxax" ,
"rtDW.jhm0cqabg3" , "rtDW.oq35dv23vw" , "rtDW.lwhmmfv2ps" , "rtDW.iyivz4czyo"
, "rtDW.ls5rwnwukt" , "rtDW.mebb30zlfb" , "rtDW.kmluc3bnxm" ,
"rtDW.ardpmlcavq" , "rtDW.o2te45t5iz" , "rtDW.ceunniigow" , "rtDW.by4chuqxsd"
, "rtDW.pjygrikwyo" , "rtDW.dbasyf4pfi" , "rtDW.nj2tur0kwz" ,
"rtDW.f1cbvtmv22" , "rtDW.mh3gjthtbl" , "rtDW.jwgslfbtic" , "rtDW.ezstggpf2j"
, "rtDW.pd45haa5mh" , "rtDW.lngqcivqqg" , "rtDW.ozmcvi5b31" ,
"rtDW.hdadv1kifl" , "rtDW.dhds3swxi0" , "rtDW.oqc0b5does" , "rtDW.ieohtob00g"
, "rtDW.eiyeel4sy0" , "rtDW.lgu0i43xpl" , "rtDW.n01cjpi1wr" ,
"rtDW.bviufihvoa" , "rtDW.i0bfhtiky2" , "rtDW.lq23ppo5t1" , "rtDW.ffetch5fn3"
, "rtDW.aj3ho3s0dk" , "rtDW.jjxvwteapb" , "rtDW.mqctzxdg0o" ,
"rtDW.nnmjponkr5" , "rtDW.ky3lis0a5h" , "rtDW.iftqygn0yl" , "rtDW.ndljvxaegf"
, "rtDW.k5hl2nqzw5" , "rtDW.exk4rrasy3" , "rtDW.oz0f0tgc3u" ,
"rtDW.h2bg1hotvh" , "rtDW.glrxfbh4lh" , "rtDW.oay4r2t35k" , "rtDW.itbalhm53p"
, "rtDW.o3hz2jqijw" , "rtDW.a1ctm05tly" , "rtDW.gib2w1lke4" ,
"rtDW.ncvj32mtr4" , "rtDW.ojwdvmmbzv" , "rtDW.blmbxswg4c" , "rtDW.eiqm0ndfyg"
, "rtDW.pibyff0eql" , "rtDW.l2ec5jammp" , "rtDW.iuadxosr40" ,
"rtDW.iqaxcni1cl" , "rtDW.o4w412lyhz" , "rtDW.l4r1nqhb1g" , "rtDW.gq3uhkidnv"
, "rtDW.emuswyjmhs" , "rtDW.eqmbwqv4uu" , "rtDW.bxfixvxfdw" ,
"rtDW.dei1252yrl" , "rtDW.ogqlv3h0ke" , "rtDW.cpgojctm31" , "rtDW.hgiqtwygy1"
, "rtDW.ddu0m0coos" , "rtDW.plosltfnse" , "rtDW.i2wzmz255s" ,
"rtDW.a1er4qtcuo" , "rtDW.chs3ptplyf" , "rtDW.jtymq0r3wm" , "rtDW.avqe5nr4ga"
, "rtDW.h53r4zk0cj" , "rtDW.gktfawr2st" , "rtDW.bj4ejpvtiv" ,
"rtDW.eyhr14ycpb" , "rtDW.i0q150xsvv" , "rtDW.inbjvmxdis" , "rtDW.lam1o4zpo4"
, "rtDW.ot4rduhsnb" , "rtDW.gmlynrr0kl" , "rtDW.dfembtpulp" ,
"rtDW.mg55yfehmm" , "rtDW.knhfhd3fs2" , "rtDW.l4gngyaono" , "rtDW.b0eijwagye"
, "rtDW.oxke0tvvta" , "rtDW.kr4a3cr0mf" , "rtDW.gdblfyh2ib" ,
"rtDW.mmq5ttry5y" , "rtDW.agwigalrld" , "rtDW.fudyiowmus" , "rtDW.ftrjzsptcb"
, "rtDW.ciu1ngvzgz" , "rtDW.csdsboi20g" , "rtDW.pxsktqyoas" ,
"rtDW.bya2xiv5j3" , "rtDW.hnl0z3aet1" , "rtDW.ecvugwr44v" , "rtDW.dfdi1srria"
, "rtDW.ibwargsoog" , "rtDW.cb3tpkjx5q" , "rtDW.fg4y55crvv" ,
"rtDW.h4wxpouwz2" , "rtDW.loncelqhha" , "rtDW.afd0kmksnl" , "rtDW.pv0jmcqzdg"
, "rtDW.lbzvytb4d4" , "rtDW.pagreemrym" , "rtDW.kkiofr3le1" ,
"rtDW.pqnrnlrnbt" , "rtDW.d1ahpk1qqr" , "rtDW.hw5n3aggbh" , "rtDW.kgduzxmrs5"
, } ; mxArray * rtdwData = mxCreateStructMatrix ( 1 , 1 , 118 ,
rtdwDataFieldNames ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 0 ,
( const void * ) & ( rtDW . ik1trwr3gq ) , sizeof ( rtDW . ik1trwr3gq ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 1 ,
( const void * ) & ( rtDW . m0rr22lxax ) , sizeof ( rtDW . m0rr22lxax ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 2 ,
( const void * ) & ( rtDW . jhm0cqabg3 ) , sizeof ( rtDW . jhm0cqabg3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 3 ,
( const void * ) & ( rtDW . oq35dv23vw ) , sizeof ( rtDW . oq35dv23vw ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 4 ,
( const void * ) & ( rtDW . lwhmmfv2ps ) , sizeof ( rtDW . lwhmmfv2ps ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 5 ,
( const void * ) & ( rtDW . iyivz4czyo ) , sizeof ( rtDW . iyivz4czyo ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 6 ,
( const void * ) & ( rtDW . ls5rwnwukt ) , sizeof ( rtDW . ls5rwnwukt ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 7 ,
( const void * ) & ( rtDW . mebb30zlfb ) , sizeof ( rtDW . mebb30zlfb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 8 ,
( const void * ) & ( rtDW . kmluc3bnxm ) , sizeof ( rtDW . kmluc3bnxm ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 9 ,
( const void * ) & ( rtDW . ardpmlcavq ) , sizeof ( rtDW . ardpmlcavq ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 10 ,
( const void * ) & ( rtDW . o2te45t5iz ) , sizeof ( rtDW . o2te45t5iz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 11 ,
( const void * ) & ( rtDW . ceunniigow ) , sizeof ( rtDW . ceunniigow ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 12 ,
( const void * ) & ( rtDW . by4chuqxsd ) , sizeof ( rtDW . by4chuqxsd ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 13 ,
( const void * ) & ( rtDW . pjygrikwyo ) , sizeof ( rtDW . pjygrikwyo ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 14 ,
( const void * ) & ( rtDW . dbasyf4pfi ) , sizeof ( rtDW . dbasyf4pfi ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 15 ,
( const void * ) & ( rtDW . nj2tur0kwz ) , sizeof ( rtDW . nj2tur0kwz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 16 ,
( const void * ) & ( rtDW . f1cbvtmv22 ) , sizeof ( rtDW . f1cbvtmv22 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 17 ,
( const void * ) & ( rtDW . mh3gjthtbl ) , sizeof ( rtDW . mh3gjthtbl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 18 ,
( const void * ) & ( rtDW . jwgslfbtic ) , sizeof ( rtDW . jwgslfbtic ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 19 ,
( const void * ) & ( rtDW . ezstggpf2j ) , sizeof ( rtDW . ezstggpf2j ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 20 ,
( const void * ) & ( rtDW . pd45haa5mh ) , sizeof ( rtDW . pd45haa5mh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 21 ,
( const void * ) & ( rtDW . lngqcivqqg ) , sizeof ( rtDW . lngqcivqqg ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 22 ,
( const void * ) & ( rtDW . ozmcvi5b31 ) , sizeof ( rtDW . ozmcvi5b31 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 23 ,
( const void * ) & ( rtDW . hdadv1kifl ) , sizeof ( rtDW . hdadv1kifl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 24 ,
( const void * ) & ( rtDW . dhds3swxi0 ) , sizeof ( rtDW . dhds3swxi0 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 25 ,
( const void * ) & ( rtDW . oqc0b5does ) , sizeof ( rtDW . oqc0b5does ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 26 ,
( const void * ) & ( rtDW . ieohtob00g ) , sizeof ( rtDW . ieohtob00g ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 27 ,
( const void * ) & ( rtDW . eiyeel4sy0 ) , sizeof ( rtDW . eiyeel4sy0 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 28 ,
( const void * ) & ( rtDW . lgu0i43xpl ) , sizeof ( rtDW . lgu0i43xpl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 29 ,
( const void * ) & ( rtDW . n01cjpi1wr ) , sizeof ( rtDW . n01cjpi1wr ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 30 ,
( const void * ) & ( rtDW . bviufihvoa ) , sizeof ( rtDW . bviufihvoa ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 31 ,
( const void * ) & ( rtDW . i0bfhtiky2 ) , sizeof ( rtDW . i0bfhtiky2 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 32 ,
( const void * ) & ( rtDW . lq23ppo5t1 ) , sizeof ( rtDW . lq23ppo5t1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 33 ,
( const void * ) & ( rtDW . ffetch5fn3 ) , sizeof ( rtDW . ffetch5fn3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 34 ,
( const void * ) & ( rtDW . aj3ho3s0dk ) , sizeof ( rtDW . aj3ho3s0dk ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 35 ,
( const void * ) & ( rtDW . jjxvwteapb ) , sizeof ( rtDW . jjxvwteapb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 36 ,
( const void * ) & ( rtDW . mqctzxdg0o ) , sizeof ( rtDW . mqctzxdg0o ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 37 ,
( const void * ) & ( rtDW . nnmjponkr5 ) , sizeof ( rtDW . nnmjponkr5 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 38 ,
( const void * ) & ( rtDW . ky3lis0a5h ) , sizeof ( rtDW . ky3lis0a5h ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 39 ,
( const void * ) & ( rtDW . iftqygn0yl ) , sizeof ( rtDW . iftqygn0yl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 40 ,
( const void * ) & ( rtDW . ndljvxaegf ) , sizeof ( rtDW . ndljvxaegf ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 41 ,
( const void * ) & ( rtDW . k5hl2nqzw5 ) , sizeof ( rtDW . k5hl2nqzw5 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 42 ,
( const void * ) & ( rtDW . exk4rrasy3 ) , sizeof ( rtDW . exk4rrasy3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 43 ,
( const void * ) & ( rtDW . oz0f0tgc3u ) , sizeof ( rtDW . oz0f0tgc3u ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 44 ,
( const void * ) & ( rtDW . h2bg1hotvh ) , sizeof ( rtDW . h2bg1hotvh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 45 ,
( const void * ) & ( rtDW . glrxfbh4lh ) , sizeof ( rtDW . glrxfbh4lh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 46 ,
( const void * ) & ( rtDW . oay4r2t35k ) , sizeof ( rtDW . oay4r2t35k ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 47 ,
( const void * ) & ( rtDW . itbalhm53p ) , sizeof ( rtDW . itbalhm53p ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 48 ,
( const void * ) & ( rtDW . o3hz2jqijw ) , sizeof ( rtDW . o3hz2jqijw ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 49 ,
( const void * ) & ( rtDW . a1ctm05tly ) , sizeof ( rtDW . a1ctm05tly ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 50 ,
( const void * ) & ( rtDW . gib2w1lke4 ) , sizeof ( rtDW . gib2w1lke4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 51 ,
( const void * ) & ( rtDW . ncvj32mtr4 ) , sizeof ( rtDW . ncvj32mtr4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 52 ,
( const void * ) & ( rtDW . ojwdvmmbzv ) , sizeof ( rtDW . ojwdvmmbzv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 53 ,
( const void * ) & ( rtDW . blmbxswg4c ) , sizeof ( rtDW . blmbxswg4c ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 54 ,
( const void * ) & ( rtDW . eiqm0ndfyg ) , sizeof ( rtDW . eiqm0ndfyg ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 55 ,
( const void * ) & ( rtDW . pibyff0eql ) , sizeof ( rtDW . pibyff0eql ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 56 ,
( const void * ) & ( rtDW . l2ec5jammp ) , sizeof ( rtDW . l2ec5jammp ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 57 ,
( const void * ) & ( rtDW . iuadxosr40 ) , sizeof ( rtDW . iuadxosr40 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 58 ,
( const void * ) & ( rtDW . iqaxcni1cl ) , sizeof ( rtDW . iqaxcni1cl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 59 ,
( const void * ) & ( rtDW . o4w412lyhz ) , sizeof ( rtDW . o4w412lyhz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 60 ,
( const void * ) & ( rtDW . l4r1nqhb1g ) , sizeof ( rtDW . l4r1nqhb1g ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 61 ,
( const void * ) & ( rtDW . gq3uhkidnv ) , sizeof ( rtDW . gq3uhkidnv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 62 ,
( const void * ) & ( rtDW . emuswyjmhs ) , sizeof ( rtDW . emuswyjmhs ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 63 ,
( const void * ) & ( rtDW . eqmbwqv4uu ) , sizeof ( rtDW . eqmbwqv4uu ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 64 ,
( const void * ) & ( rtDW . bxfixvxfdw ) , sizeof ( rtDW . bxfixvxfdw ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 65 ,
( const void * ) & ( rtDW . dei1252yrl ) , sizeof ( rtDW . dei1252yrl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 66 ,
( const void * ) & ( rtDW . ogqlv3h0ke ) , sizeof ( rtDW . ogqlv3h0ke ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 67 ,
( const void * ) & ( rtDW . cpgojctm31 ) , sizeof ( rtDW . cpgojctm31 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 68 ,
( const void * ) & ( rtDW . hgiqtwygy1 ) , sizeof ( rtDW . hgiqtwygy1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 69 ,
( const void * ) & ( rtDW . ddu0m0coos ) , sizeof ( rtDW . ddu0m0coos ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 70 ,
( const void * ) & ( rtDW . plosltfnse ) , sizeof ( rtDW . plosltfnse ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 71 ,
( const void * ) & ( rtDW . i2wzmz255s ) , sizeof ( rtDW . i2wzmz255s ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 72 ,
( const void * ) & ( rtDW . a1er4qtcuo ) , sizeof ( rtDW . a1er4qtcuo ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 73 ,
( const void * ) & ( rtDW . chs3ptplyf ) , sizeof ( rtDW . chs3ptplyf ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 74 ,
( const void * ) & ( rtDW . jtymq0r3wm ) , sizeof ( rtDW . jtymq0r3wm ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 75 ,
( const void * ) & ( rtDW . avqe5nr4ga ) , sizeof ( rtDW . avqe5nr4ga ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 76 ,
( const void * ) & ( rtDW . h53r4zk0cj ) , sizeof ( rtDW . h53r4zk0cj ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 77 ,
( const void * ) & ( rtDW . gktfawr2st ) , sizeof ( rtDW . gktfawr2st ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 78 ,
( const void * ) & ( rtDW . bj4ejpvtiv ) , sizeof ( rtDW . bj4ejpvtiv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 79 ,
( const void * ) & ( rtDW . eyhr14ycpb ) , sizeof ( rtDW . eyhr14ycpb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 80 ,
( const void * ) & ( rtDW . i0q150xsvv ) , sizeof ( rtDW . i0q150xsvv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 81 ,
( const void * ) & ( rtDW . inbjvmxdis ) , sizeof ( rtDW . inbjvmxdis ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 82 ,
( const void * ) & ( rtDW . lam1o4zpo4 ) , sizeof ( rtDW . lam1o4zpo4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 83 ,
( const void * ) & ( rtDW . ot4rduhsnb ) , sizeof ( rtDW . ot4rduhsnb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 84 ,
( const void * ) & ( rtDW . gmlynrr0kl ) , sizeof ( rtDW . gmlynrr0kl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 85 ,
( const void * ) & ( rtDW . dfembtpulp ) , sizeof ( rtDW . dfembtpulp ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 86 ,
( const void * ) & ( rtDW . mg55yfehmm ) , sizeof ( rtDW . mg55yfehmm ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 87 ,
( const void * ) & ( rtDW . knhfhd3fs2 ) , sizeof ( rtDW . knhfhd3fs2 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 88 ,
( const void * ) & ( rtDW . l4gngyaono ) , sizeof ( rtDW . l4gngyaono ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 89 ,
( const void * ) & ( rtDW . b0eijwagye ) , sizeof ( rtDW . b0eijwagye ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 90 ,
( const void * ) & ( rtDW . oxke0tvvta ) , sizeof ( rtDW . oxke0tvvta ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 91 ,
( const void * ) & ( rtDW . kr4a3cr0mf ) , sizeof ( rtDW . kr4a3cr0mf ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 92 ,
( const void * ) & ( rtDW . gdblfyh2ib ) , sizeof ( rtDW . gdblfyh2ib ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 93 ,
( const void * ) & ( rtDW . mmq5ttry5y ) , sizeof ( rtDW . mmq5ttry5y ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 94 ,
( const void * ) & ( rtDW . agwigalrld ) , sizeof ( rtDW . agwigalrld ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 95 ,
( const void * ) & ( rtDW . fudyiowmus ) , sizeof ( rtDW . fudyiowmus ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 96 ,
( const void * ) & ( rtDW . ftrjzsptcb ) , sizeof ( rtDW . ftrjzsptcb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 97 ,
( const void * ) & ( rtDW . ciu1ngvzgz ) , sizeof ( rtDW . ciu1ngvzgz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 98 ,
( const void * ) & ( rtDW . csdsboi20g ) , sizeof ( rtDW . csdsboi20g ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 99 ,
( const void * ) & ( rtDW . pxsktqyoas ) , sizeof ( rtDW . pxsktqyoas ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 100
, ( const void * ) & ( rtDW . bya2xiv5j3 ) , sizeof ( rtDW . bya2xiv5j3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 101
, ( const void * ) & ( rtDW . hnl0z3aet1 ) , sizeof ( rtDW . hnl0z3aet1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 102
, ( const void * ) & ( rtDW . ecvugwr44v ) , sizeof ( rtDW . ecvugwr44v ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 103
, ( const void * ) & ( rtDW . dfdi1srria ) , sizeof ( rtDW . dfdi1srria ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 104
, ( const void * ) & ( rtDW . ibwargsoog ) , sizeof ( rtDW . ibwargsoog ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 105
, ( const void * ) & ( rtDW . cb3tpkjx5q ) , sizeof ( rtDW . cb3tpkjx5q ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 106
, ( const void * ) & ( rtDW . fg4y55crvv ) , sizeof ( rtDW . fg4y55crvv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 107
, ( const void * ) & ( rtDW . h4wxpouwz2 ) , sizeof ( rtDW . h4wxpouwz2 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 108
, ( const void * ) & ( rtDW . loncelqhha ) , sizeof ( rtDW . loncelqhha ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 109
, ( const void * ) & ( rtDW . afd0kmksnl ) , sizeof ( rtDW . afd0kmksnl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 110
, ( const void * ) & ( rtDW . pv0jmcqzdg ) , sizeof ( rtDW . pv0jmcqzdg ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 111
, ( const void * ) & ( rtDW . lbzvytb4d4 ) , sizeof ( rtDW . lbzvytb4d4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 112
, ( const void * ) & ( rtDW . pagreemrym ) , sizeof ( rtDW . pagreemrym ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 113
, ( const void * ) & ( rtDW . kkiofr3le1 ) , sizeof ( rtDW . kkiofr3le1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 114
, ( const void * ) & ( rtDW . pqnrnlrnbt ) , sizeof ( rtDW . pqnrnlrnbt ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 115
, ( const void * ) & ( rtDW . d1ahpk1qqr ) , sizeof ( rtDW . d1ahpk1qqr ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 116
, ( const void * ) & ( rtDW . hw5n3aggbh ) , sizeof ( rtDW . hw5n3aggbh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_cacheDataAsMxArray ( rtdwData , 0 , 117
, ( const void * ) & ( rtDW . kgduzxmrs5 ) , sizeof ( rtDW . kgduzxmrs5 ) ) ;
mxSetFieldByNumber ( ssDW , 0 , 1 , rtdwData ) ; } return ssDW ; } void
mr_MagneticBasket_Simscape_Optimizer_SetDWork ( const mxArray * ssDW ) { (
void ) ssDW ; mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( (
void * ) & ( rtB ) , ssDW , 0 , 0 , sizeof ( rtB ) ) ; { const mxArray *
rtdwData = mxGetFieldByNumber ( ssDW , 0 , 1 ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ik1trwr3gq ) , rtdwData , 0 , 0 , sizeof ( rtDW . ik1trwr3gq ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . m0rr22lxax ) , rtdwData , 0 , 1 , sizeof ( rtDW . m0rr22lxax ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . jhm0cqabg3 ) , rtdwData , 0 , 2 , sizeof ( rtDW . jhm0cqabg3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . oq35dv23vw ) , rtdwData , 0 , 3 , sizeof ( rtDW . oq35dv23vw ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . lwhmmfv2ps ) , rtdwData , 0 , 4 , sizeof ( rtDW . lwhmmfv2ps ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . iyivz4czyo ) , rtdwData , 0 , 5 , sizeof ( rtDW . iyivz4czyo ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ls5rwnwukt ) , rtdwData , 0 , 6 , sizeof ( rtDW . ls5rwnwukt ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . mebb30zlfb ) , rtdwData , 0 , 7 , sizeof ( rtDW . mebb30zlfb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . kmluc3bnxm ) , rtdwData , 0 , 8 , sizeof ( rtDW . kmluc3bnxm ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ardpmlcavq ) , rtdwData , 0 , 9 , sizeof ( rtDW . ardpmlcavq ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . o2te45t5iz ) , rtdwData , 0 , 10 , sizeof ( rtDW . o2te45t5iz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ceunniigow ) , rtdwData , 0 , 11 , sizeof ( rtDW . ceunniigow ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . by4chuqxsd ) , rtdwData , 0 , 12 , sizeof ( rtDW . by4chuqxsd ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pjygrikwyo ) , rtdwData , 0 , 13 , sizeof ( rtDW . pjygrikwyo ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . dbasyf4pfi ) , rtdwData , 0 , 14 , sizeof ( rtDW . dbasyf4pfi ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . nj2tur0kwz ) , rtdwData , 0 , 15 , sizeof ( rtDW . nj2tur0kwz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . f1cbvtmv22 ) , rtdwData , 0 , 16 , sizeof ( rtDW . f1cbvtmv22 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . mh3gjthtbl ) , rtdwData , 0 , 17 , sizeof ( rtDW . mh3gjthtbl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . jwgslfbtic ) , rtdwData , 0 , 18 , sizeof ( rtDW . jwgslfbtic ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ezstggpf2j ) , rtdwData , 0 , 19 , sizeof ( rtDW . ezstggpf2j ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pd45haa5mh ) , rtdwData , 0 , 20 , sizeof ( rtDW . pd45haa5mh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . lngqcivqqg ) , rtdwData , 0 , 21 , sizeof ( rtDW . lngqcivqqg ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ozmcvi5b31 ) , rtdwData , 0 , 22 , sizeof ( rtDW . ozmcvi5b31 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . hdadv1kifl ) , rtdwData , 0 , 23 , sizeof ( rtDW . hdadv1kifl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . dhds3swxi0 ) , rtdwData , 0 , 24 , sizeof ( rtDW . dhds3swxi0 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . oqc0b5does ) , rtdwData , 0 , 25 , sizeof ( rtDW . oqc0b5does ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ieohtob00g ) , rtdwData , 0 , 26 , sizeof ( rtDW . ieohtob00g ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . eiyeel4sy0 ) , rtdwData , 0 , 27 , sizeof ( rtDW . eiyeel4sy0 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . lgu0i43xpl ) , rtdwData , 0 , 28 , sizeof ( rtDW . lgu0i43xpl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . n01cjpi1wr ) , rtdwData , 0 , 29 , sizeof ( rtDW . n01cjpi1wr ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . bviufihvoa ) , rtdwData , 0 , 30 , sizeof ( rtDW . bviufihvoa ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . i0bfhtiky2 ) , rtdwData , 0 , 31 , sizeof ( rtDW . i0bfhtiky2 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . lq23ppo5t1 ) , rtdwData , 0 , 32 , sizeof ( rtDW . lq23ppo5t1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ffetch5fn3 ) , rtdwData , 0 , 33 , sizeof ( rtDW . ffetch5fn3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . aj3ho3s0dk ) , rtdwData , 0 , 34 , sizeof ( rtDW . aj3ho3s0dk ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . jjxvwteapb ) , rtdwData , 0 , 35 , sizeof ( rtDW . jjxvwteapb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . mqctzxdg0o ) , rtdwData , 0 , 36 , sizeof ( rtDW . mqctzxdg0o ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . nnmjponkr5 ) , rtdwData , 0 , 37 , sizeof ( rtDW . nnmjponkr5 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ky3lis0a5h ) , rtdwData , 0 , 38 , sizeof ( rtDW . ky3lis0a5h ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . iftqygn0yl ) , rtdwData , 0 , 39 , sizeof ( rtDW . iftqygn0yl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ndljvxaegf ) , rtdwData , 0 , 40 , sizeof ( rtDW . ndljvxaegf ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . k5hl2nqzw5 ) , rtdwData , 0 , 41 , sizeof ( rtDW . k5hl2nqzw5 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . exk4rrasy3 ) , rtdwData , 0 , 42 , sizeof ( rtDW . exk4rrasy3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . oz0f0tgc3u ) , rtdwData , 0 , 43 , sizeof ( rtDW . oz0f0tgc3u ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . h2bg1hotvh ) , rtdwData , 0 , 44 , sizeof ( rtDW . h2bg1hotvh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . glrxfbh4lh ) , rtdwData , 0 , 45 , sizeof ( rtDW . glrxfbh4lh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . oay4r2t35k ) , rtdwData , 0 , 46 , sizeof ( rtDW . oay4r2t35k ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . itbalhm53p ) , rtdwData , 0 , 47 , sizeof ( rtDW . itbalhm53p ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . o3hz2jqijw ) , rtdwData , 0 , 48 , sizeof ( rtDW . o3hz2jqijw ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . a1ctm05tly ) , rtdwData , 0 , 49 , sizeof ( rtDW . a1ctm05tly ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . gib2w1lke4 ) , rtdwData , 0 , 50 , sizeof ( rtDW . gib2w1lke4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ncvj32mtr4 ) , rtdwData , 0 , 51 , sizeof ( rtDW . ncvj32mtr4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ojwdvmmbzv ) , rtdwData , 0 , 52 , sizeof ( rtDW . ojwdvmmbzv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . blmbxswg4c ) , rtdwData , 0 , 53 , sizeof ( rtDW . blmbxswg4c ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . eiqm0ndfyg ) , rtdwData , 0 , 54 , sizeof ( rtDW . eiqm0ndfyg ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pibyff0eql ) , rtdwData , 0 , 55 , sizeof ( rtDW . pibyff0eql ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . l2ec5jammp ) , rtdwData , 0 , 56 , sizeof ( rtDW . l2ec5jammp ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . iuadxosr40 ) , rtdwData , 0 , 57 , sizeof ( rtDW . iuadxosr40 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . iqaxcni1cl ) , rtdwData , 0 , 58 , sizeof ( rtDW . iqaxcni1cl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . o4w412lyhz ) , rtdwData , 0 , 59 , sizeof ( rtDW . o4w412lyhz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . l4r1nqhb1g ) , rtdwData , 0 , 60 , sizeof ( rtDW . l4r1nqhb1g ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . gq3uhkidnv ) , rtdwData , 0 , 61 , sizeof ( rtDW . gq3uhkidnv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . emuswyjmhs ) , rtdwData , 0 , 62 , sizeof ( rtDW . emuswyjmhs ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . eqmbwqv4uu ) , rtdwData , 0 , 63 , sizeof ( rtDW . eqmbwqv4uu ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . bxfixvxfdw ) , rtdwData , 0 , 64 , sizeof ( rtDW . bxfixvxfdw ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . dei1252yrl ) , rtdwData , 0 , 65 , sizeof ( rtDW . dei1252yrl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ogqlv3h0ke ) , rtdwData , 0 , 66 , sizeof ( rtDW . ogqlv3h0ke ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . cpgojctm31 ) , rtdwData , 0 , 67 , sizeof ( rtDW . cpgojctm31 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . hgiqtwygy1 ) , rtdwData , 0 , 68 , sizeof ( rtDW . hgiqtwygy1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ddu0m0coos ) , rtdwData , 0 , 69 , sizeof ( rtDW . ddu0m0coos ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . plosltfnse ) , rtdwData , 0 , 70 , sizeof ( rtDW . plosltfnse ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . i2wzmz255s ) , rtdwData , 0 , 71 , sizeof ( rtDW . i2wzmz255s ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . a1er4qtcuo ) , rtdwData , 0 , 72 , sizeof ( rtDW . a1er4qtcuo ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . chs3ptplyf ) , rtdwData , 0 , 73 , sizeof ( rtDW . chs3ptplyf ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . jtymq0r3wm ) , rtdwData , 0 , 74 , sizeof ( rtDW . jtymq0r3wm ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . avqe5nr4ga ) , rtdwData , 0 , 75 , sizeof ( rtDW . avqe5nr4ga ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . h53r4zk0cj ) , rtdwData , 0 , 76 , sizeof ( rtDW . h53r4zk0cj ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . gktfawr2st ) , rtdwData , 0 , 77 , sizeof ( rtDW . gktfawr2st ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . bj4ejpvtiv ) , rtdwData , 0 , 78 , sizeof ( rtDW . bj4ejpvtiv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . eyhr14ycpb ) , rtdwData , 0 , 79 , sizeof ( rtDW . eyhr14ycpb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . i0q150xsvv ) , rtdwData , 0 , 80 , sizeof ( rtDW . i0q150xsvv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . inbjvmxdis ) , rtdwData , 0 , 81 , sizeof ( rtDW . inbjvmxdis ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . lam1o4zpo4 ) , rtdwData , 0 , 82 , sizeof ( rtDW . lam1o4zpo4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ot4rduhsnb ) , rtdwData , 0 , 83 , sizeof ( rtDW . ot4rduhsnb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . gmlynrr0kl ) , rtdwData , 0 , 84 , sizeof ( rtDW . gmlynrr0kl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . dfembtpulp ) , rtdwData , 0 , 85 , sizeof ( rtDW . dfembtpulp ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . mg55yfehmm ) , rtdwData , 0 , 86 , sizeof ( rtDW . mg55yfehmm ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . knhfhd3fs2 ) , rtdwData , 0 , 87 , sizeof ( rtDW . knhfhd3fs2 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . l4gngyaono ) , rtdwData , 0 , 88 , sizeof ( rtDW . l4gngyaono ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . b0eijwagye ) , rtdwData , 0 , 89 , sizeof ( rtDW . b0eijwagye ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . oxke0tvvta ) , rtdwData , 0 , 90 , sizeof ( rtDW . oxke0tvvta ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . kr4a3cr0mf ) , rtdwData , 0 , 91 , sizeof ( rtDW . kr4a3cr0mf ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . gdblfyh2ib ) , rtdwData , 0 , 92 , sizeof ( rtDW . gdblfyh2ib ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . mmq5ttry5y ) , rtdwData , 0 , 93 , sizeof ( rtDW . mmq5ttry5y ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . agwigalrld ) , rtdwData , 0 , 94 , sizeof ( rtDW . agwigalrld ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . fudyiowmus ) , rtdwData , 0 , 95 , sizeof ( rtDW . fudyiowmus ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ftrjzsptcb ) , rtdwData , 0 , 96 , sizeof ( rtDW . ftrjzsptcb ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ciu1ngvzgz ) , rtdwData , 0 , 97 , sizeof ( rtDW . ciu1ngvzgz ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . csdsboi20g ) , rtdwData , 0 , 98 , sizeof ( rtDW . csdsboi20g ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pxsktqyoas ) , rtdwData , 0 , 99 , sizeof ( rtDW . pxsktqyoas ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . bya2xiv5j3 ) , rtdwData , 0 , 100 , sizeof ( rtDW . bya2xiv5j3 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . hnl0z3aet1 ) , rtdwData , 0 , 101 , sizeof ( rtDW . hnl0z3aet1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ecvugwr44v ) , rtdwData , 0 , 102 , sizeof ( rtDW . ecvugwr44v ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . dfdi1srria ) , rtdwData , 0 , 103 , sizeof ( rtDW . dfdi1srria ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . ibwargsoog ) , rtdwData , 0 , 104 , sizeof ( rtDW . ibwargsoog ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . cb3tpkjx5q ) , rtdwData , 0 , 105 , sizeof ( rtDW . cb3tpkjx5q ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . fg4y55crvv ) , rtdwData , 0 , 106 , sizeof ( rtDW . fg4y55crvv ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . h4wxpouwz2 ) , rtdwData , 0 , 107 , sizeof ( rtDW . h4wxpouwz2 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . loncelqhha ) , rtdwData , 0 , 108 , sizeof ( rtDW . loncelqhha ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . afd0kmksnl ) , rtdwData , 0 , 109 , sizeof ( rtDW . afd0kmksnl ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pv0jmcqzdg ) , rtdwData , 0 , 110 , sizeof ( rtDW . pv0jmcqzdg ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . lbzvytb4d4 ) , rtdwData , 0 , 111 , sizeof ( rtDW . lbzvytb4d4 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pagreemrym ) , rtdwData , 0 , 112 , sizeof ( rtDW . pagreemrym ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . kkiofr3le1 ) , rtdwData , 0 , 113 , sizeof ( rtDW . kkiofr3le1 ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . pqnrnlrnbt ) , rtdwData , 0 , 114 , sizeof ( rtDW . pqnrnlrnbt ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . d1ahpk1qqr ) , rtdwData , 0 , 115 , sizeof ( rtDW . d1ahpk1qqr ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . hw5n3aggbh ) , rtdwData , 0 , 116 , sizeof ( rtDW . hw5n3aggbh ) ) ;
mr_MagneticBasket_Simscape_Optimizer_restoreDataFromMxArray ( ( void * ) & (
rtDW . kgduzxmrs5 ) , rtdwData , 0 , 117 , sizeof ( rtDW . kgduzxmrs5 ) ) ; }
} mxArray * mr_MagneticBasket_Simscape_Optimizer_GetSimStateDisallowedBlocks
( ) { mxArray * data = mxCreateCellMatrix ( 3 , 3 ) ; mwIndex subs [ 2 ] ,
offset ; { static const char_T * blockType [ 3 ] = { "SimscapeExecutionBlock"
, "SimscapeExecutionBlock" , "SimscapeSinkBlock" , } ; static const char_T *
blockPath [ 3 ] = {
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/STATE_1" ,
"MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/OUTPUT_1_0"
, "MagneticBasket_Simscape_Optimizer/Solver Configuration/EVAL_KEY/SINK_1" ,
} ; static const int reason [ 3 ] = { 0 , 0 , 0 , } ; for ( subs [ 0 ] = 0 ;
subs [ 0 ] < 3 ; ++ ( subs [ 0 ] ) ) { subs [ 1 ] = 0 ; offset =
mxCalcSingleSubscript ( data , 2 , subs ) ; mxSetCell ( data , offset ,
mxCreateString ( blockType [ subs [ 0 ] ] ) ) ; subs [ 1 ] = 1 ; offset =
mxCalcSingleSubscript ( data , 2 , subs ) ; mxSetCell ( data , offset ,
mxCreateString ( blockPath [ subs [ 0 ] ] ) ) ; subs [ 1 ] = 2 ; offset =
mxCalcSingleSubscript ( data , 2 , subs ) ; mxSetCell ( data , offset ,
mxCreateDoubleScalar ( ( real_T ) reason [ subs [ 0 ] ] ) ) ; } } return data
; } void MdlInitializeSizes ( void ) { ssSetNumContStates ( rtS , 18 ) ;
ssSetNumPeriodicContStates ( rtS , 0 ) ; ssSetNumY ( rtS , 0 ) ; ssSetNumU (
rtS , 0 ) ; ssSetDirectFeedThrough ( rtS , 0 ) ; ssSetNumSampleTimes ( rtS ,
1 ) ; ssSetNumBlocks ( rtS , 820 ) ; ssSetNumBlockIO ( rtS , 333 ) ;
ssSetNumBlockParams ( rtS , 186 ) ; } void MdlInitializeSampleTimes ( void )
{ ssSetSampleTime ( rtS , 0 , 0.0 ) ; ssSetOffsetTime ( rtS , 0 , 0.0 ) ; }
void raccel_set_checksum ( ) { ssSetChecksumVal ( rtS , 0 , 343361676U ) ;
ssSetChecksumVal ( rtS , 1 , 2573916448U ) ; ssSetChecksumVal ( rtS , 2 ,
664120902U ) ; ssSetChecksumVal ( rtS , 3 , 3813395425U ) ; }
#if defined(_MSC_VER)
#pragma optimize( "", off )
#endif
SimStruct * raccel_register_model ( ssExecutionInfo * executionInfo ) {
static struct _ssMdlInfo mdlInfo ; static struct _ssBlkInfo2 blkInfo2 ;
static struct _ssBlkInfoSLSize blkInfoSLSize ; rt_modelMapInfoPtr = & (
rt_dataMapInfo . mmi ) ; executionInfo -> gblObjects_ . numToFiles = 0 ;
executionInfo -> gblObjects_ . numFrFiles = 0 ; executionInfo -> gblObjects_
. numFrWksBlocks = 0 ; executionInfo -> gblObjects_ . numModelInputs = 0 ;
executionInfo -> gblObjects_ . numRootInportBlks = 0 ; executionInfo ->
gblObjects_ . inportDataTypeIdx = NULL ; executionInfo -> gblObjects_ .
inportDims = NULL ; executionInfo -> gblObjects_ . inportComplex = NULL ;
executionInfo -> gblObjects_ . inportInterpoFlag = NULL ; executionInfo ->
gblObjects_ . inportContinuous = NULL ; ( void ) memset ( ( char_T * ) rtS ,
0 , sizeof ( SimStruct ) ) ; ( void ) memset ( ( char_T * ) & mdlInfo , 0 ,
sizeof ( struct _ssMdlInfo ) ) ; ( void ) memset ( ( char_T * ) & blkInfo2 ,
0 , sizeof ( struct _ssBlkInfo2 ) ) ; ( void ) memset ( ( char_T * ) &
blkInfoSLSize , 0 , sizeof ( struct _ssBlkInfoSLSize ) ) ; ssSetBlkInfo2Ptr (
rtS , & blkInfo2 ) ; ssSetBlkInfoSLSizePtr ( rtS , & blkInfoSLSize ) ;
ssSetMdlInfoPtr ( rtS , & mdlInfo ) ; ssSetExecutionInfo ( rtS ,
executionInfo ) ; slsaAllocOPModelData ( rtS ) ; { static time_T mdlPeriod [
NSAMPLE_TIMES ] ; static time_T mdlOffset [ NSAMPLE_TIMES ] ; static time_T
mdlTaskTimes [ NSAMPLE_TIMES ] ; static int_T mdlTsMap [ NSAMPLE_TIMES ] ;
static int_T mdlSampleHits [ NSAMPLE_TIMES ] ; static boolean_T
mdlTNextWasAdjustedPtr [ NSAMPLE_TIMES ] ; static int_T mdlPerTaskSampleHits
[ NSAMPLE_TIMES * NSAMPLE_TIMES ] ; static time_T mdlTimeOfNextSampleHit [
NSAMPLE_TIMES ] ; { int_T i ; for ( i = 0 ; i < NSAMPLE_TIMES ; i ++ ) {
mdlPeriod [ i ] = 0.0 ; mdlOffset [ i ] = 0.0 ; mdlTaskTimes [ i ] = 0.0 ;
mdlTsMap [ i ] = i ; mdlSampleHits [ i ] = 1 ; } } ssSetSampleTimePtr ( rtS ,
& mdlPeriod [ 0 ] ) ; ssSetOffsetTimePtr ( rtS , & mdlOffset [ 0 ] ) ;
ssSetSampleTimeTaskIDPtr ( rtS , & mdlTsMap [ 0 ] ) ; ssSetTPtr ( rtS , &
mdlTaskTimes [ 0 ] ) ; ssSetSampleHitPtr ( rtS , & mdlSampleHits [ 0 ] ) ;
ssSetTNextWasAdjustedPtr ( rtS , & mdlTNextWasAdjustedPtr [ 0 ] ) ;
ssSetPerTaskSampleHitsPtr ( rtS , & mdlPerTaskSampleHits [ 0 ] ) ;
ssSetTimeOfNextSampleHitPtr ( rtS , & mdlTimeOfNextSampleHit [ 0 ] ) ; }
ssSetSolverMode ( rtS , SOLVER_MODE_SINGLETASKING ) ; { ssSetBlockIO ( rtS ,
( ( void * ) & rtB ) ) ; ( void ) memset ( ( ( void * ) & rtB ) , 0 , sizeof
( B ) ) ; } { real_T * x = ( real_T * ) & rtX ; ssSetContStates ( rtS , x ) ;
( void ) memset ( ( void * ) x , 0 , sizeof ( X ) ) ; } { void * dwork = (
void * ) & rtDW ; ssSetRootDWork ( rtS , dwork ) ; ( void ) memset ( dwork ,
0 , sizeof ( DW ) ) ; } { static DataTypeTransInfo dtInfo ; ( void ) memset (
( char_T * ) & dtInfo , 0 , sizeof ( dtInfo ) ) ; ssSetModelMappingInfo ( rtS
, & dtInfo ) ; dtInfo . numDataTypes = 25 ; dtInfo . dataTypeSizes = &
rtDataTypeSizes [ 0 ] ; dtInfo . dataTypeNames = & rtDataTypeNames [ 0 ] ;
dtInfo . BTransTable = & rtBTransTable ; dtInfo . PTransTable = &
rtPTransTable ; dtInfo . dataTypeInfoTable = rtDataTypeInfoTable ; }
MagneticBasket_Simscape_Optimizer_InitializeDataMapInfo ( ) ;
ssSetIsRapidAcceleratorActive ( rtS , true ) ; ssSetRootSS ( rtS , rtS ) ;
ssSetVersion ( rtS , SIMSTRUCT_VERSION_LEVEL2 ) ; ssSetModelName ( rtS ,
"MagneticBasket_Simscape_Optimizer" ) ; ssSetPath ( rtS ,
"MagneticBasket_Simscape_Optimizer" ) ; ssSetTStart ( rtS , 0.0 ) ;
ssSetTFinal ( rtS , 1.0 ) ; { static RTWLogInfo rt_DataLoggingInfo ;
rt_DataLoggingInfo . loggingInterval = ( NULL ) ; ssSetRTWLogInfo ( rtS , &
rt_DataLoggingInfo ) ; } { { static int_T rt_LoggedStateWidths [ ] = { 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 ,
2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2
, 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 }
; static int_T rt_LoggedStateNumDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,
1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1
, 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 } ; static int_T
rt_LoggedStateDimensions [ ] = { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 ,
1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2
, 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 ,
2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 , 2 } ; static boolean_T
rt_LoggedStateIsVarDims [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static BuiltInDTypeId
rt_LoggedStateDataTypeIds [ ] = { SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE , SS_DOUBLE ,
SS_DOUBLE , SS_DOUBLE , SS_DOUBLE } ; static int_T
rt_LoggedStateComplexSignals [ ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static RTWPreprocessingFcnPtr
rt_LoggingStatePreprocessingFcnPtrs [ ] = { ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) , ( NULL ) ,
( NULL ) } ; static const char_T * rt_LoggedStateLabels [ ] = { "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" , "CSTATE" ,
"CSTATE" , "CSTATE" , "CSTATE" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" } ; static const char_T *
rt_LoggedStateBlockNames [ ] = {
"MagneticBasket_Simscape_Optimizer/Prismatic Joint1" ,
"MagneticBasket_Simscape_Optimizer/Prismatic Joint1" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint8" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint8" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint7" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint7" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint6" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint6" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint1" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint1" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint2" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint2" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint3" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint3" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint4" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint4" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint5" ,
"MagneticBasket_Simscape_Optimizer/Subsystem1/Revolute Joint5" ,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_1_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_1_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_1_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_2_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_2_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_2_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_3_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_3_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_3_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_4_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_4_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_4_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_5_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_5_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_5_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_6_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_6_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_6_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_7_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_7_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_7_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_8_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_8_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_8_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_9_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_9_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_9_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_10_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_10_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_10_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_11_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_11_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_11_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_12_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_12_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_12_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_13_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_13_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_13_1_3"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_14_1_1"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_14_1_2"
,
"MagneticBasket_Simscape_Optimizer/Solver\nConfiguration/EVAL_KEY/INPUT_14_1_3"
} ; static const char_T * rt_LoggedStateNames [ ] = {
"MagneticBasket_Simscape_Optimizer.Prismatic_Joint1.Pz.p" ,
"MagneticBasket_Simscape_Optimizer.Prismatic_Joint1.Pz.v" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint8.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint8.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint7.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint7.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint6.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint6.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint1.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint1.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint2.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint2.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint3.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint3.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint4.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint4.Rz.w" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint5.Rz.q" ,
"MagneticBasket_Simscape_Optimizer.Subsystem1.Revolute_Joint5.Rz.w" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" ,
"Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" , "Discrete" }
; static boolean_T rt_LoggedStateCrossMdlRef [ ] = { 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0
, 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ; static
RTWLogDataTypeConvert rt_RTWLogDataTypeConvert [ ] = { { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE ,
SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0
, 0 , 1.0 , 0 , 0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 ,
0.0 } , { 0 , SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } , { 0 ,
SS_DOUBLE , SS_DOUBLE , 0 , 0 , 0 , 1.0 , 0 , 0.0 } } ; static int_T
rt_LoggedStateIdxList [ ] = { 0 , 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10
, 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 , 24 , 25 ,
26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 , 37 , 38 , 39 , 40 , 41
} ; static RTWLogSignalInfo rt_LoggedStateSignalInfo = { 60 ,
rt_LoggedStateWidths , rt_LoggedStateNumDimensions , rt_LoggedStateDimensions
, rt_LoggedStateIsVarDims , ( NULL ) , ( NULL ) , rt_LoggedStateDataTypeIds ,
rt_LoggedStateComplexSignals , ( NULL ) , rt_LoggingStatePreprocessingFcnPtrs
, { rt_LoggedStateLabels } , ( NULL ) , ( NULL ) , ( NULL ) , {
rt_LoggedStateBlockNames } , { rt_LoggedStateNames } ,
rt_LoggedStateCrossMdlRef , rt_RTWLogDataTypeConvert , rt_LoggedStateIdxList
} ; static void * rt_LoggedStateSignalPtrs [ 60 ] ; rtliSetLogXSignalPtrs (
ssGetRTWLogInfo ( rtS ) , ( LogSignalPtrsType ) rt_LoggedStateSignalPtrs ) ;
rtliSetLogXSignalInfo ( ssGetRTWLogInfo ( rtS ) , & rt_LoggedStateSignalInfo
) ; rt_LoggedStateSignalPtrs [ 0 ] = ( void * ) & rtX . eouclruk0d [ 0 ] ;
rt_LoggedStateSignalPtrs [ 1 ] = ( void * ) & rtX . eouclruk0d [ 1 ] ;
rt_LoggedStateSignalPtrs [ 2 ] = ( void * ) & rtX . eouclruk0d [ 2 ] ;
rt_LoggedStateSignalPtrs [ 3 ] = ( void * ) & rtX . eouclruk0d [ 3 ] ;
rt_LoggedStateSignalPtrs [ 4 ] = ( void * ) & rtX . eouclruk0d [ 4 ] ;
rt_LoggedStateSignalPtrs [ 5 ] = ( void * ) & rtX . eouclruk0d [ 5 ] ;
rt_LoggedStateSignalPtrs [ 6 ] = ( void * ) & rtX . eouclruk0d [ 6 ] ;
rt_LoggedStateSignalPtrs [ 7 ] = ( void * ) & rtX . eouclruk0d [ 7 ] ;
rt_LoggedStateSignalPtrs [ 8 ] = ( void * ) & rtX . eouclruk0d [ 8 ] ;
rt_LoggedStateSignalPtrs [ 9 ] = ( void * ) & rtX . eouclruk0d [ 9 ] ;
rt_LoggedStateSignalPtrs [ 10 ] = ( void * ) & rtX . eouclruk0d [ 10 ] ;
rt_LoggedStateSignalPtrs [ 11 ] = ( void * ) & rtX . eouclruk0d [ 11 ] ;
rt_LoggedStateSignalPtrs [ 12 ] = ( void * ) & rtX . eouclruk0d [ 12 ] ;
rt_LoggedStateSignalPtrs [ 13 ] = ( void * ) & rtX . eouclruk0d [ 13 ] ;
rt_LoggedStateSignalPtrs [ 14 ] = ( void * ) & rtX . eouclruk0d [ 14 ] ;
rt_LoggedStateSignalPtrs [ 15 ] = ( void * ) & rtX . eouclruk0d [ 15 ] ;
rt_LoggedStateSignalPtrs [ 16 ] = ( void * ) & rtX . eouclruk0d [ 16 ] ;
rt_LoggedStateSignalPtrs [ 17 ] = ( void * ) & rtX . eouclruk0d [ 17 ] ;
rt_LoggedStateSignalPtrs [ 18 ] = ( void * ) rtDW . ik1trwr3gq ;
rt_LoggedStateSignalPtrs [ 19 ] = ( void * ) rtDW . m0rr22lxax ;
rt_LoggedStateSignalPtrs [ 20 ] = ( void * ) rtDW . jhm0cqabg3 ;
rt_LoggedStateSignalPtrs [ 21 ] = ( void * ) rtDW . oq35dv23vw ;
rt_LoggedStateSignalPtrs [ 22 ] = ( void * ) rtDW . lwhmmfv2ps ;
rt_LoggedStateSignalPtrs [ 23 ] = ( void * ) rtDW . iyivz4czyo ;
rt_LoggedStateSignalPtrs [ 24 ] = ( void * ) rtDW . ls5rwnwukt ;
rt_LoggedStateSignalPtrs [ 25 ] = ( void * ) rtDW . mebb30zlfb ;
rt_LoggedStateSignalPtrs [ 26 ] = ( void * ) rtDW . kmluc3bnxm ;
rt_LoggedStateSignalPtrs [ 27 ] = ( void * ) rtDW . ardpmlcavq ;
rt_LoggedStateSignalPtrs [ 28 ] = ( void * ) rtDW . o2te45t5iz ;
rt_LoggedStateSignalPtrs [ 29 ] = ( void * ) rtDW . ceunniigow ;
rt_LoggedStateSignalPtrs [ 30 ] = ( void * ) rtDW . by4chuqxsd ;
rt_LoggedStateSignalPtrs [ 31 ] = ( void * ) rtDW . pjygrikwyo ;
rt_LoggedStateSignalPtrs [ 32 ] = ( void * ) rtDW . dbasyf4pfi ;
rt_LoggedStateSignalPtrs [ 33 ] = ( void * ) rtDW . nj2tur0kwz ;
rt_LoggedStateSignalPtrs [ 34 ] = ( void * ) rtDW . f1cbvtmv22 ;
rt_LoggedStateSignalPtrs [ 35 ] = ( void * ) rtDW . mh3gjthtbl ;
rt_LoggedStateSignalPtrs [ 36 ] = ( void * ) rtDW . jwgslfbtic ;
rt_LoggedStateSignalPtrs [ 37 ] = ( void * ) rtDW . ezstggpf2j ;
rt_LoggedStateSignalPtrs [ 38 ] = ( void * ) rtDW . pd45haa5mh ;
rt_LoggedStateSignalPtrs [ 39 ] = ( void * ) rtDW . lngqcivqqg ;
rt_LoggedStateSignalPtrs [ 40 ] = ( void * ) rtDW . ozmcvi5b31 ;
rt_LoggedStateSignalPtrs [ 41 ] = ( void * ) rtDW . hdadv1kifl ;
rt_LoggedStateSignalPtrs [ 42 ] = ( void * ) rtDW . dhds3swxi0 ;
rt_LoggedStateSignalPtrs [ 43 ] = ( void * ) rtDW . oqc0b5does ;
rt_LoggedStateSignalPtrs [ 44 ] = ( void * ) rtDW . ieohtob00g ;
rt_LoggedStateSignalPtrs [ 45 ] = ( void * ) rtDW . eiyeel4sy0 ;
rt_LoggedStateSignalPtrs [ 46 ] = ( void * ) rtDW . lgu0i43xpl ;
rt_LoggedStateSignalPtrs [ 47 ] = ( void * ) rtDW . n01cjpi1wr ;
rt_LoggedStateSignalPtrs [ 48 ] = ( void * ) rtDW . bviufihvoa ;
rt_LoggedStateSignalPtrs [ 49 ] = ( void * ) rtDW . i0bfhtiky2 ;
rt_LoggedStateSignalPtrs [ 50 ] = ( void * ) rtDW . lq23ppo5t1 ;
rt_LoggedStateSignalPtrs [ 51 ] = ( void * ) rtDW . ffetch5fn3 ;
rt_LoggedStateSignalPtrs [ 52 ] = ( void * ) rtDW . aj3ho3s0dk ;
rt_LoggedStateSignalPtrs [ 53 ] = ( void * ) rtDW . jjxvwteapb ;
rt_LoggedStateSignalPtrs [ 54 ] = ( void * ) rtDW . mqctzxdg0o ;
rt_LoggedStateSignalPtrs [ 55 ] = ( void * ) rtDW . nnmjponkr5 ;
rt_LoggedStateSignalPtrs [ 56 ] = ( void * ) rtDW . ky3lis0a5h ;
rt_LoggedStateSignalPtrs [ 57 ] = ( void * ) rtDW . iftqygn0yl ;
rt_LoggedStateSignalPtrs [ 58 ] = ( void * ) rtDW . ndljvxaegf ;
rt_LoggedStateSignalPtrs [ 59 ] = ( void * ) rtDW . k5hl2nqzw5 ; }
rtliSetLogT ( ssGetRTWLogInfo ( rtS ) , "tout" ) ; rtliSetLogX (
ssGetRTWLogInfo ( rtS ) , "" ) ; rtliSetLogXFinal ( ssGetRTWLogInfo ( rtS ) ,
"xFinal" ) ; rtliSetLogVarNameModifier ( ssGetRTWLogInfo ( rtS ) , "none" ) ;
rtliSetLogFormat ( ssGetRTWLogInfo ( rtS ) , 4 ) ; rtliSetLogMaxRows (
ssGetRTWLogInfo ( rtS ) , 0 ) ; rtliSetLogDecimation ( ssGetRTWLogInfo ( rtS
) , 1 ) ; rtliSetLogY ( ssGetRTWLogInfo ( rtS ) , "" ) ;
rtliSetLogYSignalInfo ( ssGetRTWLogInfo ( rtS ) , ( NULL ) ) ;
rtliSetLogYSignalPtrs ( ssGetRTWLogInfo ( rtS ) , ( NULL ) ) ; } { static
struct _ssStatesInfo2 statesInfo2 ; ssSetStatesInfo2 ( rtS , & statesInfo2 )
; } { static ssPeriodicStatesInfo periodicStatesInfo ;
ssSetPeriodicStatesInfo ( rtS , & periodicStatesInfo ) ; } { static
ssJacobianPerturbationBounds jacobianPerturbationBounds ;
ssSetJacobianPerturbationBounds ( rtS , & jacobianPerturbationBounds ) ; } {
static ssSolverInfo slvrInfo ; static struct _ssSFcnModelMethods3 mdlMethods3
; static struct _ssSFcnModelMethods2 mdlMethods2 ; static boolean_T
contStatesDisabled [ 18 ] ; static real_T absTol [ 18 ] = { 1.0E-6 , 1.0E-6 ,
1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 ,
1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 , 1.0E-6 } ;
static uint8_T absTolControl [ 18 ] = { 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U
, 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U , 0U } ; static real_T
contStateJacPerturbBoundMinVec [ 18 ] ; static real_T
contStateJacPerturbBoundMaxVec [ 18 ] ; static uint8_T zcAttributes [ 28 ] =
{ ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) ,
( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , (
ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) , ( ZC_EVENT_ALL ) } ; {
int i ; for ( i = 0 ; i < 18 ; ++ i ) { contStateJacPerturbBoundMinVec [ i ]
= 0 ; contStateJacPerturbBoundMaxVec [ i ] = rtGetInf ( ) ; } }
ssSetSolverRelTol ( rtS , 0.001 ) ; ssSetStepSize ( rtS , 0.0 ) ;
ssSetMinStepSize ( rtS , 0.0 ) ; ssSetMaxNumMinSteps ( rtS , - 1 ) ;
ssSetMinStepViolatedError ( rtS , 0 ) ; ssSetMaxStepSize ( rtS , 0.001 ) ;
ssSetSolverMaxOrder ( rtS , - 1 ) ; ssSetSolverRefineFactor ( rtS , 1 ) ;
ssSetOutputTimes ( rtS , ( NULL ) ) ; ssSetNumOutputTimes ( rtS , 0 ) ;
ssSetOutputTimesOnly ( rtS , 0 ) ; ssSetOutputTimesIndex ( rtS , 0 ) ;
ssSetZCCacheNeedsReset ( rtS , 0 ) ; ssSetDerivCacheNeedsReset ( rtS , 0 ) ;
ssSetNumNonContDerivSigInfos ( rtS , 0 ) ; ssSetNonContDerivSigInfos ( rtS ,
( NULL ) ) ; ssSetSolverInfo ( rtS , & slvrInfo ) ; ssSetSolverName ( rtS ,
"VariableStepAuto" ) ; ssSetVariableStepSolver ( rtS , 1 ) ;
ssSetSolverConsistencyChecking ( rtS , 0 ) ; ssSetSolverAdaptiveZcDetection (
rtS , 0 ) ; ssSetSolverRobustResetMethod ( rtS , 0 ) ;
_ssSetSolverUpdateJacobianAtReset ( rtS , true ) ; ssSetAbsTolVector ( rtS ,
absTol ) ; ssSetAbsTolControlVector ( rtS , absTolControl ) ;
ssSetSolverAbsTol_Obsolete ( rtS , absTol ) ;
ssSetSolverAbsTolControl_Obsolete ( rtS , absTolControl ) ;
ssSetJacobianPerturbationBoundsMinVec ( rtS , contStateJacPerturbBoundMinVec
) ; ssSetJacobianPerturbationBoundsMaxVec ( rtS ,
contStateJacPerturbBoundMaxVec ) ; ssSetSolverStateProjection ( rtS , 1 ) ; (
void ) memset ( ( void * ) & mdlMethods2 , 0 , sizeof ( mdlMethods2 ) ) ;
ssSetModelMethods2 ( rtS , & mdlMethods2 ) ; ( void ) memset ( ( void * ) &
mdlMethods3 , 0 , sizeof ( mdlMethods3 ) ) ; ssSetModelMethods3 ( rtS , &
mdlMethods3 ) ; ssSetModelProjection ( rtS , MdlProjection ) ;
ssSetSolverMassMatrixType ( rtS , ( ssMatrixType ) 0 ) ;
ssSetSolverMassMatrixNzMax ( rtS , 0 ) ; ssSetModelOutputs ( rtS , MdlOutputs
) ; ssSetModelUpdate ( rtS , MdlUpdate ) ; ssSetModelDerivatives ( rtS ,
MdlDerivatives ) ; ssSetSolverZcSignalAttrib ( rtS , zcAttributes ) ;
ssSetSolverNumZcSignals ( rtS , 28 ) ; ssSetModelZeroCrossings ( rtS ,
MdlZeroCrossings ) ; ssSetSolverConsecutiveZCsStepRelTol ( rtS ,
2.8421709430404007E-13 ) ; ssSetSolverMaxConsecutiveZCs ( rtS , 1000 ) ;
ssSetSolverConsecutiveZCsError ( rtS , 2 ) ; ssSetSolverMaskedZcDiagnostic (
rtS , 1 ) ; ssSetSolverIgnoredZcDiagnostic ( rtS , 1 ) ;
ssSetSolverMaxConsecutiveMinStep ( rtS , 1 ) ;
ssSetSolverShapePreserveControl ( rtS , 2 ) ; ssSetTNextTid ( rtS , INT_MIN )
; ssSetTNext ( rtS , rtMinusInf ) ; ssSetSolverNeedsReset ( rtS ) ;
ssSetNumNonsampledZCs ( rtS , 28 ) ; ssSetContStateDisabled ( rtS ,
contStatesDisabled ) ; ssSetSolverMaxConsecutiveMinStep ( rtS , 1 ) ; }
ssSetChecksumVal ( rtS , 0 , 343361676U ) ; ssSetChecksumVal ( rtS , 1 ,
2573916448U ) ; ssSetChecksumVal ( rtS , 2 , 664120902U ) ; ssSetChecksumVal
( rtS , 3 , 3813395425U ) ; { static const sysRanDType rtAlwaysEnabled =
SUBSYS_RAN_BC_ENABLE ; static RTWExtModeInfo rt_ExtModeInfo ; static const
sysRanDType * systemRan [ 29 ] ; gblRTWExtModeInfo = & rt_ExtModeInfo ;
ssSetRTWExtModeInfo ( rtS , & rt_ExtModeInfo ) ;
rteiSetSubSystemActiveVectorAddresses ( & rt_ExtModeInfo , systemRan ) ;
systemRan [ 0 ] = & rtAlwaysEnabled ; systemRan [ 1 ] = & rtAlwaysEnabled ;
systemRan [ 2 ] = & rtAlwaysEnabled ; systemRan [ 3 ] = & rtAlwaysEnabled ;
systemRan [ 4 ] = & rtAlwaysEnabled ; systemRan [ 5 ] = & rtAlwaysEnabled ;
systemRan [ 6 ] = & rtAlwaysEnabled ; systemRan [ 7 ] = & rtAlwaysEnabled ;
systemRan [ 8 ] = & rtAlwaysEnabled ; systemRan [ 9 ] = & rtAlwaysEnabled ;
systemRan [ 10 ] = & rtAlwaysEnabled ; systemRan [ 11 ] = & rtAlwaysEnabled ;
systemRan [ 12 ] = & rtAlwaysEnabled ; systemRan [ 13 ] = & rtAlwaysEnabled ;
systemRan [ 14 ] = & rtAlwaysEnabled ; systemRan [ 15 ] = & rtAlwaysEnabled ;
systemRan [ 16 ] = & rtAlwaysEnabled ; systemRan [ 17 ] = & rtAlwaysEnabled ;
systemRan [ 18 ] = & rtAlwaysEnabled ; systemRan [ 19 ] = & rtAlwaysEnabled ;
systemRan [ 20 ] = & rtAlwaysEnabled ; systemRan [ 21 ] = & rtAlwaysEnabled ;
systemRan [ 22 ] = & rtAlwaysEnabled ; systemRan [ 23 ] = & rtAlwaysEnabled ;
systemRan [ 24 ] = & rtAlwaysEnabled ; systemRan [ 25 ] = & rtAlwaysEnabled ;
systemRan [ 26 ] = & rtAlwaysEnabled ; systemRan [ 27 ] = & rtAlwaysEnabled ;
systemRan [ 28 ] = & rtAlwaysEnabled ; rteiSetModelMappingInfoPtr (
ssGetRTWExtModeInfo ( rtS ) , & ssGetModelMappingInfo ( rtS ) ) ;
rteiSetChecksumsPtr ( ssGetRTWExtModeInfo ( rtS ) , ssGetChecksums ( rtS ) )
; rteiSetTPtr ( ssGetRTWExtModeInfo ( rtS ) , ssGetTPtr ( rtS ) ) ; }
slsaDisallowedBlocksForSimTargetOP ( rtS ,
mr_MagneticBasket_Simscape_Optimizer_GetSimStateDisallowedBlocks ) ;
slsaGetWorkFcnForSimTargetOP ( rtS ,
mr_MagneticBasket_Simscape_Optimizer_GetDWork ) ;
slsaSetWorkFcnForSimTargetOP ( rtS ,
mr_MagneticBasket_Simscape_Optimizer_SetDWork ) ;
rt_RapidReadMatFileAndUpdateParams ( rtS ) ; if ( ssGetErrorStatus ( rtS ) )
{ return rtS ; } return rtS ; }
#if defined(_MSC_VER)
#pragma optimize( "", on )
#endif
void MdlOutputsParameterSampleTime ( int_T tid ) { MdlOutputsTID1 ( tid ) ; }
