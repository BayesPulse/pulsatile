///-----------------------------------------------------------------------------
///
/// FILE: test_inout.h
///
/// DESCRIPTION: 
///   Function definitions for test_inout.c
/// 
///-----------------------------------------------------------------------------

#include <R.h>
#include <Rinternals.h>

// Don't recall what this is for
//#ifndef TEST_INOUT_H 
//#define TEST_INOUT_H 

// Function Definitions
SEXP testc(SEXP indata);
SEXP testspec(SEXP inspec);
SEXP getListElement(SEXP list, const char *str);
SEXP showArgs1(SEXP largs);

