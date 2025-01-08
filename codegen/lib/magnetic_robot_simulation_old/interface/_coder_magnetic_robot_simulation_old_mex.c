/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_magnetic_robot_simulation_old_mex.c
 *
 * MATLAB Coder version            : 5.5
 * C/C++ source code generated on  : 06-Jan-2025 18:33:03
 */

/* Include Files */
#include "_coder_magnetic_robot_simulation_old_mex.h"
#include "_coder_magnetic_robot_simulation_old_api.h"

/* Function Definitions */
/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  (void)plhs;
  (void)prhs;
  mexAtExit(&magnetic_robot_simulation_old_atexit);
  /* Module initialization. */
  magnetic_robot_simulation_old_initialize();
  /* Dispatch the entry-point. */
  unsafe_magnetic_robot_simulation_old_mexFunction(nlhs, nrhs);
  /* Module termination. */
  magnetic_robot_simulation_old_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-949", true);
  return emlrtRootTLSGlobal;
}

/*
 * Arguments    : int32_T nlhs
 *                int32_T nrhs
 * Return Type  : void
 */
void unsafe_magnetic_robot_simulation_old_mexFunction(int32_T nlhs,
                                                      int32_T nrhs)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 0, 4,
                        29, "magnetic_robot_simulation_old");
  }
  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 29,
                        "magnetic_robot_simulation_old");
  }
  /* Call the function. */
  c_magnetic_robot_simulation_old();
}

/*
 * File trailer for _coder_magnetic_robot_simulation_old_mex.c
 *
 * [EOF]
 */
