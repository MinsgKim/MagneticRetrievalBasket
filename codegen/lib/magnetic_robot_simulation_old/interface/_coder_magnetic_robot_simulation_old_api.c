/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_magnetic_robot_simulation_old_api.c
 *
 * MATLAB Coder version            : 5.5
 * C/C++ source code generated on  : 06-Jan-2025 18:33:03
 */

/* Include Files */
#include "_coder_magnetic_robot_simulation_old_api.h"
#include "_coder_magnetic_robot_simulation_old_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131627U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "magnetic_robot_simulation_old",                      /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void c_magnetic_robot_simulation_old(void)
{
  /* Invoke the target function */
  magnetic_robot_simulation_old();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void magnetic_robot_simulation_old_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  magnetic_robot_simulation_old_xil_terminate();
  magnetic_robot_simulation_old_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void magnetic_robot_simulation_old_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void magnetic_robot_simulation_old_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_magnetic_robot_simulation_old_api.c
 *
 * [EOF]
 */
