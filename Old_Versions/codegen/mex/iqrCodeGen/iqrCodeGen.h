/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * iqrCodeGen.h
 *
 * Code generation for function 'iqrCodeGen'
 *
 */

#ifndef __IQRCODEGEN_H__
#define __IQRCODEGEN_H__

/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "iqrCodeGen_types.h"

/* Function Declarations */
extern real_T iqrCodeGen(const emlrtStack *sp, const real_T x[100]);

#ifdef __WATCOMC__

#pragma aux iqrCodeGen value [8087];

#endif
#endif

/* End of code generation (iqrCodeGen.h) */
