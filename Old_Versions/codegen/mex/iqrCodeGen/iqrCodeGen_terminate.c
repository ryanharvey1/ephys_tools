/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * iqrCodeGen_terminate.c
 *
 * Code generation for function 'iqrCodeGen_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "iqrCodeGen.h"
#include "iqrCodeGen_terminate.h"
#include "iqrCodeGen_data.h"

/* Function Definitions */
void iqrCodeGen_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void iqrCodeGen_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (iqrCodeGen_terminate.c) */
