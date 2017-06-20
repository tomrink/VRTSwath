/* Generated automatically from netcdf3.h.in by configure. */
/*
 *	Copyright 1993, University Corporation for Atmospheric Research
 *
 *  Permission to use, copy, modify, and distribute this software and its
 * documentation for any purpose without fee is hereby granted, provided
 * that the above copyright notice appear in all copies, that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of UCAR/Unidata not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.  UCAR makes no
 * representations about the suitability of this software for any purpose.
 * It is provided "as is" without express or implied warranty.  It is
 * provided with no support and without obligation on the part of UCAR
 * Unidata, to assist in its use, correction, modification, or enhancement.
 *
 */
/* "$Id: netcdf3.h,v 1.1 2008/09/04 16:14:28 mmerritt Exp $" */

#ifndef _NETCDF3_
#define _NETCDF3_

#ifdef __MWERKS__
#ifndef HDF
#define HDF
#endif
#endif /* __MWERKS__ */

/*
 *   If xdr_enum works properly on your system, you can define 
 * USE_ENUM so that nc_type is an enum. 
 * Otherwise, delete this definition so that the nc_type is
 * an int and the valid values are #defined.
 */
#ifndef __MWERKS__
#define USE_ENUM
#endif



#ifndef HDF
/*
 * This can be as large as the maximum number of stdio streams
 * you can have open on your system.
 */
#define MAX_NC_OPEN 32

/*
 * These maximums are enforced by the interface, to facilitate writing
 * applications and utilities.  However, nothing is statically allocated to
 * these sizes internally.
 */
#define MAX_NC_DIMS 5000	 /* max dimensions per file */
#define MAX_NC_ATTRS 3000	 /* max global or per variable attributes */
#define MAX_NC_VARS 5000	 /* max variables per file */
#define MAX_NC_NAME 256		 /* max length of a name */
#define MAX_VAR_DIMS 32          /* max per variable dimensions */

/*
 * Added feature. 
 * If you wish a variable to use a different value than the above
 * defaults, create an attribute with the same type as the variable
 * and the following reserved name. The value you give the attribute
 * will be used as the fill value for that variable.
 */
#define _FillValue	"_FillValue"

#else /* HDF */

#include "hlimits.h"  /* Hard coded constants for HDF library */

#endif /* HDF */













/*
 * NB: The following feature-test line is too long in order to accomodate a 
 * bug in the VMS 5.3 C compiler.
 */
#ifndef HAVE_PROTOTYPES
#   if defined(__STDC__) || defined(__GNUC__) || defined(__cplusplus) || defined(c_plusplus)
#       define	HAVE_PROTOTYPES
#   endif
#endif

#undef PROTO
#ifdef HAVE_PROTOTYPES 
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif


#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif /* _NETCDF_ */
