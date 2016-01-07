/**
 * Type definitions for the Real data type, which is either single or double
 * precision depending on the build.
 */

#ifndef METROPOLIS_DATA_TYPES_H
#define METROPOLIS_DATA_TYPES_H

#ifdef SINGLE_PRECISION
typedef float Real;
#else
typedef double Real;
#endif

#endif
