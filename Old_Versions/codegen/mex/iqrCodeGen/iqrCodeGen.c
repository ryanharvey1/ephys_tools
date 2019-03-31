/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * iqrCodeGen.c
 *
 * Code generation for function 'iqrCodeGen'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "iqrCodeGen.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 5, "iqrCodeGen",
  "/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analyses/Old Versions/iqrCodeGen.m"
};

static emlrtRSInfo b_emlrtRSI = { 9, "iqr",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/iqr.m" };

static emlrtRSInfo c_emlrtRSI = { 54, "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m" };

static emlrtRSInfo d_emlrtRSI = { 134, "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m" };

static emlrtRSInfo e_emlrtRSI = { 136, "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m" };

static emlrtBCInfo emlrtBCI = { 1, 100, 166, 15, "", "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m", 0 };

static emlrtBCInfo b_emlrtBCI = { 1, 100, 151, 36, "", "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m", 0 };

static emlrtBCInfo c_emlrtBCI = { 1, 100, 151, 58, "", "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m", 0 };

static emlrtBCInfo d_emlrtBCI = { 1, 100, 148, 26, "", "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m", 0 };

static emlrtBCInfo e_emlrtBCI = { 1, 100, 140, 18, "", "prctile",
  "/Applications/MATLAB_R2015a.app/toolbox/stats/eml/prctile.m", 0 };

/* Function Definitions */
real_T iqrCodeGen(const emlrtStack *sp, const real_T x[100])
{
  real_T r;
  int32_T iwork[100];
  int32_T idx[100];
  int32_T k;
  boolean_T p;
  int32_T i;
  int32_T nj;
  int32_T j;
  int32_T pEnd;
  int32_T b_p;
  int32_T q;
  int32_T qEnd;
  int32_T kEnd;
  boolean_T exitg1;
  real_T y[2];
  int32_T b_idx;
  int32_T b_j;
  int32_T c_j;
  int32_T c_idx;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;

  /* IQRCODEGEN Estimate interquartile range  */
  /*    iqrCodeGen returns the interquartile range of the data x, a single- */
  /*    or double-precision vector. */
  st.site = &emlrtRSI;
  b_st.site = &b_emlrtRSI;
  c_st.site = &c_emlrtRSI;
  d_st.site = &d_emlrtRSI;
  for (k = 0; k < 100; k += 2) {
    if ((x[k] <= x[k + 1]) || muDoubleScalarIsNaN(x[k + 1])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  i = 2;
  while (i < 100) {
    nj = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 101; pEnd = qEnd + i) {
      b_p = j;
      q = pEnd - 1;
      qEnd = j + nj;
      if (qEnd > 101) {
        qEnd = 101;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if ((x[idx[b_p - 1] - 1] <= x[idx[q] - 1]) || muDoubleScalarIsNaN
            (x[idx[q] - 1])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (q + 1 < qEnd) {
              k++;
              iwork[k] = idx[q];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q];
          q++;
          if (q + 1 == qEnd) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }

      j = qEnd;
    }

    i = nj;
  }

  d_st.site = &e_emlrtRSI;
  nj = 100;
  exitg1 = false;
  while ((!exitg1) && (nj > 0)) {
    j = idx[nj - 1];
    emlrtDynamicBoundsCheckR2012b(j, 1, 100, &emlrtBCI, &d_st);
    if (!muDoubleScalarIsNaN(x[idx[nj - 1] - 1])) {
      exitg1 = true;
    } else {
      nj--;
    }
  }

  if (nj < 1) {
    for (i = 0; i < 2; i++) {
      y[i] = rtNaN;
    }
  } else if (nj == 1) {
    for (i = 0; i < 2; i++) {
      if ((idx[0] >= 1) && (idx[0] < 100)) {
        b_idx = idx[0];
      } else {
        b_idx = emlrtDynamicBoundsCheckR2012b(idx[0], 1, 100, &e_emlrtBCI, &c_st);
      }

      y[i] = x[b_idx - 1];
    }
  } else {
    for (k = 0; k < 2; k++) {
      r = (25.0 + 50.0 * (real_T)k) / 100.0 * (real_T)nj;
      i = (int32_T)muDoubleScalarRound(r);
      if (nj <= i) {
        j = idx[nj - 1];
        if ((j >= 1) && (j < 100)) {
          b_j = j;
        } else {
          b_j = emlrtDynamicBoundsCheckR2012b(j, 1, 100, &d_emlrtBCI, &c_st);
        }

        y[k] = x[b_j - 1];
      } else {
        r -= (real_T)i;
        j = idx[i - 1];
        if ((j >= 1) && (j < 100)) {
          c_j = j;
        } else {
          c_j = emlrtDynamicBoundsCheckR2012b(j, 1, 100, &b_emlrtBCI, &c_st);
        }

        if ((idx[i] >= 1) && (idx[i] < 100)) {
          c_idx = idx[i];
        } else {
          c_idx = emlrtDynamicBoundsCheckR2012b(idx[i], 1, 100, &c_emlrtBCI,
            &c_st);
        }

        y[k] = (0.5 - r) * x[c_j - 1] + (0.5 + r) * x[c_idx - 1];
      }
    }
  }

  return y[1] - y[0];
}

/* End of code generation (iqrCodeGen.c) */
