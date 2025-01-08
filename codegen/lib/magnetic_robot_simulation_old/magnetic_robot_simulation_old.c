/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: magnetic_robot_simulation_old.c
 *
 * MATLAB Coder version            : 5.5
 * C/C++ source code generated on  : 06-Jan-2025 18:33:03
 */

/* Include Files */
#include "magnetic_robot_simulation_old.h"
#include "magnetic_robot_simulation_old_emxutil.h"
#include "magnetic_robot_simulation_old_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Declarations */
static void c_magnetic_robot_simulation_old(const double y[6],
                                            double varargout_1[6]);

static int div_nde_s32_floor(int numerator);

static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
/*
 * Arguments    : const double y[6]
 *                double varargout_1[6]
 * Return Type  : void
 */
static void c_magnetic_robot_simulation_old(const double y[6],
                                            double varargout_1[6])
{
  static const double dv[3] = {0.0, 1.5707963267948966, 0.0};
  int i;
  /*  State vector Y = [theta1, theta2, theta3, omega1, omega2, omega3] */
  /*  Initialize acceleration vector */
  /*  Calculate torques for each link */
  for (i = 0; i < 3; i++) {
    double tau_spring;
    double varargout_1_tmp;
    /*  Magnetic torque */
    /*  Spring torques */
    tau_spring = 0.0;
    if (i + 1 > 1) {
      tau_spring = 0.0 - 0.01 * (y[i] - y[i - 1]);
    }
    if (i + 1 < 3) {
      tau_spring += 0.01 * (y[i + 1] - y[i]);
    }
    /*  Simple damping */
    /*  Sum all torques and calculate acceleration */
    /*  Simple moment of inertia for a rod */
    varargout_1_tmp = y[i + 3];
    varargout_1[i] = varargout_1_tmp;
    varargout_1[i + 3] =
        ((0.01 * sin((y[i] + dv[i]) - 0.78539816339744817) + tau_spring) +
         -0.01 * varargout_1_tmp) /
        8.3333333333333337E-6;
  }
}

/*
 * Arguments    : int numerator
 * Return Type  : int
 */
static int div_nde_s32_floor(int numerator)
{
  int i;
  if ((numerator < 0) && (numerator % 6 != 0)) {
    i = -1;
  } else {
    i = 0;
  }
  return numerator / 6 + i;
}

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = fabs(u0);
    d1 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }
  return y;
}

/*
 * Parameters
 *
 * Arguments    : void
 * Return Type  : void
 */
void magnetic_robot_simulation_old(void)
{
  static const double x[21] = {0.2,
                               0.075,
                               0.225,
                               0.97777777777777775,
                               -3.7333333333333334,
                               3.5555555555555554,
                               2.9525986892242035,
                               -11.595793324188385,
                               9.8228928516994358,
                               -0.29080932784636487,
                               2.8462752525252526,
                               -10.757575757575758,
                               8.9064227177434727,
                               0.27840909090909088,
                               -0.2735313036020583,
                               0.091145833333333329,
                               0.0,
                               0.44923629829290207,
                               0.65104166666666663,
                               -0.322376179245283,
                               0.13095238095238096};
  static const double b[7] = {0.0012326388888888888,
                              0.0,
                              -0.0042527702905061394,
                              0.036979166666666667,
                              -0.05086379716981132,
                              0.0419047619047619,
                              -0.025};
  static const double f0[6] = {
      0.0, 0.0, 0.0, -848.528137423857, 848.52813742385706, -848.528137423857};
  emxArray_char_T *str;
  emxArray_real_T *tout;
  double f[42];
  double y[6];
  double absh;
  double t;
  double *tout_data;
  int Bcolidx;
  int b_i;
  int exponent;
  int i;
  int i1;
  int ia;
  int iac;
  int j;
  int nout;
  char *str_data;
  boolean_T Done;
  boolean_T MinStepExit;
  /*  Main script for magnetic robot simulation */
  /*  1cm links */
  /*  Spring constants between links [N*m/rad] */
  /*  Magnetization magnitude for each link [A*m^2] */
  /*  Initial magnetization angles relative to link orientation [rad] */
  /*  Initial configuration */
  /*  Initial link angles [rad] */
  /*  Time parameters */
  /*  External magnetic field parameters */
  /*  Tesla */
  /*  rad */
  /*  Solve equations of motion using ode45 */
  emxInit_real_T(&tout);
  i = tout->size[0] * tout->size[1];
  tout->size[0] = 1;
  tout->size[1] = 200;
  emxEnsureCapacity_real_T(tout, i);
  tout_data = tout->data;
  for (i = 0; i < 200; i++) {
    tout_data[i] = 0.0;
  }
  nout = 1;
  tout_data[0] = 0.0;
  absh = 5.9487228922853474E-5;
  t = 0.0;
  for (b_i = 0; b_i < 6; b_i++) {
    y[b_i] = 0.0;
  }
  memset(&f[0], 0, 42U * sizeof(double));
  for (i = 0; i < 6; i++) {
    f[i] = f0[i];
  }
  MinStepExit = false;
  Done = false;
  double absx;
  int exitg1;
  do {
    double ystage[6];
    double h;
    double hmin;
    double tnew;
    boolean_T NoFailedAttempts;
    exitg1 = 0;
    absx = fabs(t);
    if (rtIsInf(absx) || rtIsNaN(absx)) {
      absx = rtNaN;
    } else if (absx < 4.4501477170144028E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
    hmin = 16.0 * absx;
    absh = fmin(1.0, fmax(hmin, absh));
    h = absh;
    absx = fabs(10.0 - t);
    if (1.1 * absh >= absx) {
      h = 10.0 - t;
      absh = absx;
      Done = true;
    }
    NoFailedAttempts = true;
    int exitg2;
    do {
      double fE[6];
      exitg2 = 0;
      Bcolidx = 0;
      for (j = 0; j < 5; j++) {
        Bcolidx += j;
        for (b_i = 0; b_i < 6; b_i++) {
          ystage[b_i] = y[b_i];
        }
        if (!(h == 0.0)) {
          i = 6 * j + 1;
          for (iac = 1; iac <= i; iac += 6) {
            absx = h * x[Bcolidx + div_nde_s32_floor(iac - 1)];
            i1 = iac + 5;
            for (ia = iac; ia <= i1; ia++) {
              b_i = ia - iac;
              ystage[b_i] += f[ia - 1] * absx;
            }
          }
        }
        c_magnetic_robot_simulation_old(ystage,
                                        *(double(*)[6]) & f[6 * (j + 1)]);
      }
      tnew = t + h;
      for (b_i = 0; b_i < 6; b_i++) {
        ystage[b_i] = y[b_i];
      }
      if (!(h == 0.0)) {
        for (iac = 0; iac <= 30; iac += 6) {
          absx = h * x[(Bcolidx + div_nde_s32_floor(iac)) + 5];
          i = iac + 6;
          for (ia = iac + 1; ia <= i; ia++) {
            b_i = (ia - iac) - 1;
            ystage[b_i] += f[ia - 1] * absx;
          }
        }
      }
      c_magnetic_robot_simulation_old(ystage, *(double(*)[6]) & f[36]);
      for (i = 0; i < 6; i++) {
        absx = 0.0;
        for (i1 = 0; i1 < 7; i1++) {
          absx += f[i + 6 * i1] * b[i1];
        }
        fE[i] = absx;
      }
      if (Done) {
        tnew = 10.0;
      }
      absx = 0.0;
      for (b_i = 0; b_i < 6; b_i++) {
        double d1;
        double d2;
        h = fabs(fE[b_i]);
        d1 = fabs(y[b_i]);
        d2 = fabs(ystage[b_i]);
        if ((d1 > d2) || rtIsNaN(d2)) {
          if (d1 > 1.0) {
            h /= d1;
          }
        } else if (d2 > 1.0) {
          h /= d2;
        }
        if ((h > absx) || rtIsNaN(h)) {
          absx = h;
        }
      }
      h = absh * absx;
      if (!(h <= 1.0E-6)) {
        if (absh <= hmin) {
          MinStepExit = true;
          exitg2 = 1;
        } else {
          if (NoFailedAttempts) {
            NoFailedAttempts = false;
            absh = fmax(hmin,
                        absh * fmax(0.1, 0.8 * rt_powd_snf(1.0E-6 / h, 0.2)));
          } else {
            absh = fmax(hmin, 0.5 * absh);
          }
          h = absh;
          Done = false;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);
    if (MinStepExit) {
      exitg1 = 1;
    } else {
      b_i = nout;
      absx = tnew - t;
      nout += 4;
      if (nout > tout->size[1]) {
        Bcolidx = tout->size[1];
        i = tout->size[0] * tout->size[1];
        tout->size[0] = 1;
        tout->size[1] += 200;
        emxEnsureCapacity_real_T(tout, i);
        tout_data = tout->data;
        for (j = 0; j < 200; j++) {
          tout_data[Bcolidx + j] = 0.0;
        }
      }
      tout_data[b_i] = t + absx * 0.25;
      tout_data[b_i + 1] = t + absx * 0.5;
      tout_data[b_i + 2] = t + absx * 0.75;
      tout_data[b_i + 3] = tnew;
      if (Done) {
        exitg1 = 1;
      } else {
        if (NoFailedAttempts) {
          absx = 1.25 * rt_powd_snf(h / 1.0E-6, 0.2);
          if (absx > 0.2) {
            absh /= absx;
          } else {
            absh *= 5.0;
          }
        }
        t = tnew;
        for (b_i = 0; b_i < 6; b_i++) {
          y[b_i] = ystage[b_i];
          absx = f[b_i + 36];
          ystage[b_i] = absx;
          f[b_i] = absx;
        }
      }
    }
  } while (exitg1 == 0);
  if (nout < 1) {
    Bcolidx = 0;
  } else {
    Bcolidx = nout;
  }
  /*  Animate results */
  if (Bcolidx >= 1) {
    b_i = Bcolidx;
  } else {
    b_i = 1;
  }
  if (Bcolidx == 0) {
    b_i = 0;
  }
  i = (int)(((double)b_i + 99.0) / 100.0);
  emxInit_char_T(&str);
  for (b_i = 0; b_i < i; b_i++) {
    unsigned int c_i;
    c_i = (unsigned int)b_i * 100U + 1U;
    /*  Plot robot links */
    /*  Plot magnetic field direction */
    /*  Set axes limits */
    Bcolidx =
        (int)snprintf(NULL, 0, "Time: %.2f s", tout_data[(int)c_i - 1]) + 1;
    i1 = str->size[0] * str->size[1];
    str->size[0] = 1;
    str->size[1] = Bcolidx;
    emxEnsureCapacity_char_T(str, i1);
    str_data = str->data;
    snprintf(&str_data[0], (size_t)Bcolidx, "Time: %.2f s",
             tout_data[(int)c_i - 1]);
    /*  pause(0.01); */
  }
  emxFree_char_T(&str);
  emxFree_real_T(&tout);
}

/*
 * File trailer for magnetic_robot_simulation_old.c
 *
 * [EOF]
 */
