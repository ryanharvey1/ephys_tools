/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * iqrCodeGen_initialize.c
 *
 * Code generation for function 'iqrCodeGen_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "iqrCodeGen.h"
#include "iqrCodeGen_initialize.h"
#include "iqrCodeGen_data.h"

/* Function Definitions */
void iqrCodeGen_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtLicenseCheckR2012b(&st, "Statistics_Toolbox", 2);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (iqrCodeGen_initialize.c) */
