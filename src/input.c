/*
!C****************************************************************************

!File: input.c
  
!Description: Functions reading data from the input data file.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.1 2002/05/02
 Robert Wolfe
 Added handling for SDS's with ranks greater than 2.

 Revision 1.5 2002/12/02
 Gail Schmidt
 Added support for INT8 data types.

 Revision 2.0 2003/10/15
 Gail Schmidt
 GRID types are not supported for MRTSwath

 Revision 2.0 2003/12/13
 Gail Schmidt
 Modified MIN_LS_DIM_SIZE to 250 from 1000, since some of the direct
   broadcast products have fewer than 1000 lines.

 Revision 2.0a 2004/08/21
 Gail Schmidt
 Read the _FillValue and pass it through for resampling.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Laboratory for Terrestrial Physics (Code 922) at the 
  National Aeronautics and Space Administration, Goddard Space Flight 
  Center, under NASA Task 92-012-00.

 ! References and Credits:
  ! MODIS Science Team Member:
      Christopher O. Justice
      MODIS Land Science Team           University of Maryland
      justice@hermes.geog.umd.edu       Dept. of Geography
      phone: 301-405-1600               1113 LeFrak Hall
                                        College Park, MD, 20742

  ! Developers:
      Robert E. Wolfe (Code 922)
      MODIS Land Team Support Group     Raytheon ITSS
      robert.e.wolfe.1@gsfc.nasa.gov    4400 Forbes Blvd.
      phone: 301-614-5508               Lanham, MD 20770  
  
 ! Design Notes:
   1. The following public functions handle the input data:

	OpenInput - Setup 'input' data structure and open file for access.
	CloseInput - Close the input file.
	FreeOutput - Free the 'input' data structure memory.
        FindInputDim - Determines the line/sample dimension values.

   2. 'OpenInput' must be called before any of the other routines.  
   3. 'FreeInput' should be used to free the 'input' data structure.
   4. The only input file type supported is HDF.

!END****************************************************************************
*/

#include <stdlib.h>
#include "input.h"
#include "myerror.h"
#include "mystring.h"
#include "hdf.h"
#include "mfhdf.h"
#include "HdfEosDef.h"
#include "geoloc.h"

#include "hdf5.h"
#include "resamp.h"

/* Constants */
#define FILL_VALUE_NAME "_FillValue"
#define VALID_MIN_NAME "valid_min"
#define VALID_MAX_NAME "valid_max"


Input_t *OpenInput(char *file_name, char *geoProductName, char *sds_name, int iband, int rank,
                   int *dim, char *errstr)
/* 
!C******************************************************************************

!Description: 'OpenInput' sets up the 'input' data structure, opens the
 input file for read access.
 
!Input Parameters:
 file_name      input file name
 sds_name       name of sds to be read
 iband          band number for application of band offset
 rank           rank of the input SDS
 dim            dimension flags; the line and sample dimensions are 
                indicated by a -1 and -2 in this array, respectively;
		the index in the other dimensions are indicated by a value 
		of zero or greater

!Output Parameters:
 dim            dimension flags
 errstr         error string if an error occurred
 (returns)      'input' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. When 'OpenInput' returns, the file is open for HDF access and the 
      SDS is open for access.
   2. For an input space type of 'GRID_SPACE', a band number of -1 should be 
      given.
   3. The only input HDF data types supported are CHAR8, INT8, UINT8, INT16,
      UINT16, INT32, and UINT32.
   4. An error status is returned when:
       a. the SDS rank is less than 2 or greater than 'MYHDF_MAX_RANK'
       b. the band number is less than -1 or greater than or equal to
          'NBAND_OFFSET'
       c. either none or more than one dimension is given for the line 
          or sample dimension
       d. an invalid dimension field is given
       e. duplicating strings is not successful
       f. errors occur when opening the input HDF file
       g. errors occur when reading SDS dimensions or attributes
       h. errors occur when opening the SDS for read access
       i. the given SDS rank or dimensions do not match the input file
       j. for an input space type of SWATH_SPACE, the dimensions of a swath 
          are not 1, 2 or 4 times the nominal size of a MODIS swath
       k. for an input space type of SWATH_SPACE, the number of lines is not 
          an integral multiple of the size of the scan at the given resolution
       l. memory allocation is not successful
       m. an invalid input data type is not supported.
   5. Error messages are handled with the 'LOG_RETURN_ERROR' macro.
   6. 'FreeInput' should be called to deallocate memory used by the 
      'input' data structures.
   7. 'CloseInput' should be called after all of the data is written and 
      before the 'input' data structure memory is released.

!END****************************************************************************
*/
{
  Input_t *this;
  char *error_string = (char *)NULL;
  int ir, ir1;
  char tmperrstr[M_MSG_LEN+1];
  double fill[MYHDF_MAX_NATTR_VAL];
  Myhdf_attr_t attr;
  
  /* Check parameters */
  
  if (rank < 2  ||  rank > MYHDF_MAX_RANK) {
    strcpy (errstr, "OpenInput: invalid rank");
    return (Input_t *)NULL;
  }
  
  if (iband < -1  ||  iband >= NBAND_OFFSET) {
    strcpy (errstr, "OpenInput: invalid band");
    return (Input_t *)NULL;
  }

  /* Create the Input data structure */

  this = (Input_t *)malloc(sizeof(Input_t));
  if (this == (Input_t *)NULL) {
    strcpy (errstr, "OpenInput: allocating Input data structure");
    return (Input_t *)NULL;
  }

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    strcpy (errstr, "OpenInput: duplicating file name");
    return (Input_t *)NULL;
  }

  this->sds.name = DupString(sds_name);
  if (this->sds.name == (char *)NULL) {
    free(this->file_name);
    free(this);
    strcpy (errstr, "OpenInput: duplicating sds name");
    return (Input_t *)NULL;
  }

  /* Open file for SD access */

  this->sds_file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (this->sds_file_id < 0) {
    free(this->sds.name);
    free(this->file_name);
    free(this);  
    strcpy (errstr, "OpenInput: opening input file");
    return (Input_t *)NULL;
  }
  this->open = true;
  
  /* Get SDS information and start SDS access */

  if (!GetSDSInfoV(this->sds_file_id, &this->sds)) {
    free(this->sds.name);
    free(this->file_name);
    free(this);
    strcpy (errstr, "OpenInput: getting sds info");
    return (Input_t *)NULL;
  }

  /* Get dimensions */

  for (ir = 0; ir < this->sds.rank; ir++) {
    if (!GetSDSDimInfoV(this->sds.id, &this->sds.dim[ir], ir)) {
      for (ir1 = 0; ir1 < ir; ir1++) free(this->sds.dim[ir1].name);
      free(this->sds.name);
      free(this->file_name);
      free(this);
      strcpy (errstr, "OpenInput: getting dimension");
      return (Input_t *)NULL;
    }
  }

  /* Check the rank and dimensions */

  if (this->sds.rank != rank) 
    error_string = "expected rank does not match";
    
  if (error_string == (char *)NULL)
  {
    if (!FindInputDimV(rank, dim, this->sds.dim, this->extra_dim, &this->dim,
        tmperrstr)) {
      for (ir1 = 0; ir1 < ir; ir1++) free(this->sds.dim[ir1].name);
      free(this->sds.name);
      free(this->file_name);
      free(this);
      sprintf (errstr, "%s\nOpenInput: unable to determine input line and "
        "sample dimensions", tmperrstr);
      return (Input_t *)NULL;
    }
  }

  /* Check the line and sample dimensions */

  if (error_string == (char *)NULL) {
    this->size.l = this->sds.dim[this->dim.l].nval;
    this->size.s = this->sds.dim[this->dim.s].nval;
    this->ires = 1; // always 
  }

  /* Check for valid range */
  
  double attrVal[1];
  this->hasValidRange = false;
  if (GetSingleValueAttrAsDouble(this->sds.id, VALID_MIN_NAME, attrVal)) {
      this->valid_min = attrVal[0];
  
      if (GetSingleValueAttrAsDouble(this->sds.id, VALID_MAX_NAME, attrVal)) {
          this->valid_max = attrVal[0];
          this->hasValidRange = true;
      }
  }
  
  if (this->hasValidRange) {
      printf("Dataset valid range: %f, %f\n", this->valid_min, this->valid_max);
  }
  
  /* Get fill value or use 0.0 as the fill */
  this->fill_value = 0.0;
  if (GetSingleValueAttrAsDouble(this->sds.id, FILL_VALUE_NAME, attrVal)) {
      this->fill_value = attrVal[0];
      printf("Dataset fill value:  %f\n", this->fill_value);
  }
  else {
      printf("Using default fill value: %f\n", this->fill_value);
  }


  /* Calculate number of scans, and for swath space, check for 
     an integral number of scans */

  this->scan_size.l = getNDETinScanFromGeoProductName(geoProductName);
  this->scan_size.l *= this->ires;
  this->scan_size.s = this->size.s;
  this->nscan = (this->size.l - 1) / this->scan_size.l + 1;
  if ((this->nscan * this->scan_size.l) != this->size.l) {
    for (ir1 = 0; ir1 < ir; ir1++) free(this->sds.dim[ir1].name);
    free(this->sds.name);
    free(this->file_name);
    free(this);
    strcpy (errstr, "OpenInput: not an integral number of scans");
    return (Input_t *)NULL;
  }

  /* Allocate input buffer */
  
  hid_t type = getNativeType(this->sds.typeh5);
  if (type == NULL) {
      return NULL;
  }
  
  this->data_type_size = sizeof(this->sds.datasize);
  void *buf = (void *) calloc(this->size.s, this->data_type_size);
  if (buf == NULL) {
      error_string = "allocating input i/o buffer";
      return NULL;
  }
  
  if (type == H5T_NATIVE_SCHAR) {
      this->sds.type = DFNT_INT8;
      this->buf.val_int8 = (int8 *) buf;
  }
  else if (type == H5T_NATIVE_UCHAR) {
      this->sds.type = DFNT_UINT8;
      this->buf.val_uint8 = (uint8 *) buf;
  }
  else if (type == H5T_NATIVE_USHORT) {
      this->sds.type = DFNT_UINT16;
      this->buf.val_uint16 = (uint16 *) buf;
  }
  else if (type == H5T_NATIVE_SHORT) {
      this->sds.type = DFNT_INT16;
      this->buf.val_int16 = (int16 *) buf;
  }
  else if (type == H5T_NATIVE_UINT) {
      this->sds.type = DFNT_UINT32;
      this->buf.val_uint32 = (uint32 *) buf;
  }
  else if (type == H5T_NATIVE_INT) {
      this->sds.type = DFNT_INT32;
      this->buf.val_int32 = (int32 *) buf;
  }
  else if (type == H5T_NATIVE_FLOAT) {
      this->sds.type = DFNT_FLOAT32;
      this->buf.val_float32 = (float32 *) buf;
  }
  else if (type == H5T_NATIVE_DOUBLE) {
      this->sds.type = DFNT_FLOAT64;
      this->buf.val_double = (float64 *) buf;
  }
    
  return this;
}


#define MIN_LS_DIM_SIZE (250)

bool FindInputDimV(int rank, int *param_dim, Myhdf_dim_t *sds_dim, 
                  int *extra_dim, Img_coord_int_t *dim, char *errstr)
/* 
!C******************************************************************************

!Description: 'FindInputDim' distinguishes between the line/sample and other 
              dimensions from the user parameters and the size of the
	      dimensions.
 
!Input Parameters:
 rank           rank of the input SDS
 param_dim      dimension flags from the user parameters; the line and 
                sample dimensions are 
                indicated by a -1 and -2 in this array, respectively;
		the index in the other dimensions are indicated by a value 
		of zero or greater
 file_dim       actual SDS dimensions

!Output Parameters:
 param_dim      updated dimensions
 extra_dim      extra SDS dimensions; other than line and sample
 dim            SDS line and sample dimension locations
 errstr         error string
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the dimensions from the user parameters are not valid
       b. there are either too many large (not line/sample) or small
          (line/sample) dimensions.
   2. Error messages are returned as part of the function call. This
      allows the calling routine (OpenInput) to print those error
      messages or ignore them.

!END****************************************************************************
*/
{
  int ir, ils, iextra;
  int temp_dim[MYHDF_MAX_RANK];

  /* Set up the default values to be returned */

  for (ir = 0; ir < rank; ir++) extra_dim[ir] = 0;
  dim->l = dim->s = -1;

  /* Check to make sure the line/sample dimensions from the user parameters 
     are the first two dimesions */

  if (param_dim[0] >= 0  ||  
      param_dim[1] >= 0  ||
      (param_dim[0] * param_dim[1]) != 2) {
    strcpy (errstr, "FindInputDim: invalid line/sample dimensions");
    return false;
  }

  /* Check to make sure that the remaining dimensions from the user parameters
     are valid */

  for (ir = 2; ir < rank; ir++) {
    if (param_dim[ir] < 0) {
      strcpy (errstr, "FindInputDim: invalid remaining dimesions");
      return false;
    }
  }

  /* Figure out which are the line/sample dimensions and which are the
     extra dimensions.  The line and sample dimensions are expected to be
     greater than MIN_LS_DIM_SIZE.  If too many dimensions are "line/sample"
     dimensions, then it is an error.  If not enough dimensions are
     available for "line/sample" dimensions, then it is also an error. */

  // TDR for VIIRS:
  param_dim[0] = -1;
  param_dim[1] = -2;

  dim->l = 0;
  dim->s = 1;  

  /* Update the user parameters */

  for (ir = 0; ir < rank; ir++)
    param_dim[ir] = temp_dim[ir];
  
  return true;
}


bool CloseInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'CloseInput' ends SDS access and closes the input file.
 
!Input Parameters:
 this           'input' data structure; the following fields are input:
                   open, sds.id, sds_file_id

!Output Parameters:
 this           'input' data structure; the following fields are modified:
                   open
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. an error occurs when closing access to the SDS.
   2. Error messages are handled with the 'LOG_RETURN_ERROR' macro.
   3. 'OpenInput' must be called before this routine is called.
   4. 'FreeInput' should be called to deallocate memory used by the 
      'input' data structure.

!END****************************************************************************
*/
{
 Myhdf_sds_t *sds;

 if (!this->open)
    LOG_RETURN_ERROR("file not open", "CloseInput", false);
 
    sds = &this->sds;
    H5Tclose(sds->typeh5);
    H5Dclose(sds->id);
    H5Fclose(this->sds_file_id); 

    this->open = false;

  return true;
}


bool FreeInput(Input_t *this)
/* 
!C******************************************************************************

!Description: 'FreeInput' frees the 'input' data structure memory.
 
!Input Parameters:
 this           'input' data structure; the following fields are input:
                   sds.rank, sds.dim[*].name, sds.name, file_name

!Output Parameters:
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'OpenInput' and 'CloseInput' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  int ir;

  if (this != (Input_t *)NULL) {
    for (ir = 0; ir < this->sds.rank; ir++) {
      if (this->sds.dim[ir].name != (char *)NULL) 
        free(this->sds.dim[ir].name);
    }
    if (this->sds.name != (char *)NULL) free(this->sds.name);
    if (this->file_name != (char *)NULL) free(this->file_name);
    free(this);
  }

  return true;
}

