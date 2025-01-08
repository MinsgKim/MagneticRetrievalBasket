/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: magnetic_robot_simulation_old_emxutil.h
 *
 * MATLAB Coder version            : 5.5
 * C/C++ source code generated on  : 06-Jan-2025 18:33:03
 */

#ifndef MAGNETIC_ROBOT_SIMULATION_OLD_EMXUTIL_H
#define MAGNETIC_ROBOT_SIMULATION_OLD_EMXUTIL_H

/* Include Files */
#include "magnetic_robot_simulation_old_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);

extern void emxFree_char_T(emxArray_char_T **pEmxArray);

extern void emxFree_real_T(emxArray_real_T **pEmxArray);

extern void emxInit_char_T(emxArray_char_T **pEmxArray);

extern void emxInit_real_T(emxArray_real_T **pEmxArray);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for magnetic_robot_simulation_old_emxutil.h
 *
 * [EOF]
 */
