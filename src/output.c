/*
!C****************************************************************************

!File: output.c
  
!Description: Functions creating and writting data to the HDF output file.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.5 2002/12/02
 Gail Schmidt
 Added support for INT8 data types.

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
   1. The following public functions handle the HDF output files:

        CreateOutput - Create new output file.
	OutputFile - Setup 'output' data structure.
	CloseOutput - Close the output file.
	FreeOutput - Free the 'output' data structure memory.
	WriteOutput - Write a line of data to the output product file.

   2. 'OutputFile' must be called before any of the other routines (except for 
      'CreateOutput').  
   3. 'FreeOutput' should be used to free the 'output' data structure.
   4. The only output file type supported is HDF.

!END****************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include "output.h"
#include "myerror.h"
#include "myproj.h"
#include "const.h"
#include <netcdf.h>

#define UNITS "units"
#define STANDARD_NAME "standard_name"
#define PROJECTION_X_COORDINATE "projection_x_coordinate"
#define PROJECTION_Y_COORDINATE "projection_y_coordinate"
#define GRID_MAPPING "grid_mapping"
#define GRID_PROJECTION "grid_projection"
#define SEMI_MAJOR_AXIS "semi_major_axis"
#define SEMI_MINOR_AXIS "semi_minor_axis"
#define GRID_MAPPING_NAME "grid_mapping_name"

#define CENTER_LON "longitude_of_projection_origin"
#define CENTER_LAT "latitude_of_projection_origin"
#define FALSE_EASTING "false_easting"
#define FALSE_NORTHING "false_northing"
#define LAMAZ "lambert_azimuthal_equal_area"

#define GENERIC "generic"
#define METER "meter"

#define CONVENTIONS "Conventions"
#define CF_1_7 "CF-1.7"

bool CreateOutput(char *file_name)
/* 
!C******************************************************************************

!Description: 'CreateOuptut' creates a new HDF output file.
 
!Input Parameters:
 file_name      output file name

!Output Parameters:
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the creation of the output file failes.
   2. Error messages are handled with the 'LOG_RETURN_ERROR' macro.
   3. The output file is in HDF format and closed after it is created.

!END****************************************************************************
*/
{
  int ncid, retval;

  retval = nc_create(file_name, 0x1000, &ncid);
  if (retval != 0) {
     printf("Error: %s\n", nc_strerror(retval));
     LOG_RETURN_ERROR("creating output file", "CreateOutput", false);
  }

  /* Close the file */

  retval = nc_close(ncid);
  if (retval != 0) {
     printf("Error: %s\n", nc_strerror(retval));
     LOG_RETURN_ERROR("creating output file", "CreateOutput", false);
  }

  return true;
}

Output_t *OutputFile(char *file_name, char *sds_name, 
                     int output_data_type, Space_def_t *space_def)
/* 
!C******************************************************************************

!Description: 'OutputFile' sets up the 'output' data structure, opens the
 output file for write access, and creates the output Science Data Set (SDS).
 
!Input Parameters:
 file_name      output file name
 sds_name       name of sds to be created
 output_data_type  output HDF data type; data types currently supported are
                     CHAR8, UINT8, INT8, INT16, UINT16, INT32, and UINT32.
 space_def      output space definition; the following fields are input:
                   img_size

!Output Parameters:
 (returns)      'output' data structure or NULL when an error occurs

!Team Unique Header:

 ! Design Notes:
   1. When 'OutputFile' returns sucessfully, the file is open for HDF access 
      and the SDS is open for access.
   2. An error status is returned when:
       a. the output image dimensions are zero or negative
       b. an invalid output data type is passed
       c. memory allocation is not successful
       d. duplicating strings is not successful
       e. errors occur when opening the output HDF file
       f. errors occur when creating and initializing the SDS.
   3. Error messages are handled with the 'LOG_RETURN_ERROR' macro.
   4. 'FreeOutput' should be called to deallocate memory used by the 
      'output' data structures.
   5. 'CloseFile' should be called after all of the data is written and 
      before the 'output' data structure memory is released.

!END****************************************************************************
*/
{
  Output_t *this;
  int ir, k;
  char *error_string = (char *)NULL;
  char tmpstr[1024];
  
  int retval, ncid, x_dimid, y_dimid, varid, x_varid, y_varid, proj_varid;
  int dimids[2], xdimids[1], ydimids[1], dimid_1d[1];
  int dimid;
  int start[1];
  int nval[1];

  /* Check parameters */
  
  if (space_def->img_size.l < 1)
    LOG_RETURN_ERROR("invalid number of output lines", 
                 "OutputFile", (Output_t *)NULL);

  if (space_def->img_size.s < 1)
    LOG_RETURN_ERROR("invalid number of samples per output line", 
                 "OutputFile", (Output_t *)NULL);

  if (output_data_type != DFNT_CHAR8  &&
      output_data_type != DFNT_UINT8  &&
      output_data_type != DFNT_INT8  &&
      output_data_type != DFNT_INT16  &&
      output_data_type != DFNT_UINT16 &&
      output_data_type != DFNT_INT32  &&
      output_data_type != DFNT_UINT32 &&
      output_data_type != DFNT_FLOAT32 &&
      output_data_type != DFNT_FLOAT64)
    LOG_RETURN_ERROR("output data type not supported", "OpenOutput", 
                 (Output_t *)NULL);

  /* Create the Output data structure */

  this = (Output_t *)malloc(sizeof(Output_t));
  if (this == (Output_t *)NULL) 
    LOG_RETURN_ERROR("allocating Output data structure", "OpenOutput", 
                 (Output_t *)NULL);

  /* Populate the data structure */

  this->file_name = DupString(file_name);
  if (this->file_name == (char *)NULL) {
    free(this);
    LOG_RETURN_ERROR("duplicating file name", "OutputFile", (Output_t *)NULL);
  }
  
  /* keep only leaf of original dataset path name*/
  const char *delim = "/";
  char *name = strdup(sds_name);
  char *tok;
  tok = strtok(name, delim);
  char *dataset_name;
  while (tok != NULL) {
      dataset_name = strdup(tok);
      tok = strtok(NULL, delim);
  }

  this->sds.name = dataset_name;
  if (this->sds.name == (char *)NULL) {
    free(this->file_name);
    free(this);
    LOG_RETURN_ERROR("duplicating sds name", "OutputFile", (Output_t *)NULL);
  }

  this->size.l = space_def->img_size.l;
  this->size.s = space_def->img_size.s;

  /* Open file for SD access */

  retval = nc_open(file_name, NC_WRITE, &ncid);
  if (retval != 0) {
    free(this->sds.name);
    free(this->file_name);
    free(this);  
    LOG_RETURN_ERROR("opening output file access", "OutputFile", (Output_t *)NULL);
  }
  else {
     this->sds_file_id = ncid;
     this->open = true;
  } 

  /* Set up SDS */

  this->sds.type = output_data_type;
  this->sds.rank = 2;
  this->sds.dim[0].nval = this->size.l;
  this->sds.dim[1].nval = this->size.s;
  this->sds.dim[0].name = DupString("y");
  this->sds.dim[1].name = DupString("x");
  
  /* Define the dimensions. NetCDF will hand back an ID for each. */
  retval = nc_def_dim(ncid, this->sds.dim[1].name, this->size.s, &x_dimid);
  retval = nc_def_dim(ncid, this->sds.dim[0].name, this->size.l, &y_dimid);
  retval = nc_def_dim(ncid, GENERIC, 1, &dimid);

  /* The dimids array is used to pass the IDs of the dimensions of
   * the variable. */
  dimids[0] = y_dimid;
  dimids[1] = x_dimid;
  dimid_1d[0] = dimid;
  
  nc_type nctype;

  if (output_data_type == DFNT_INT8) {
      nctype = NC_BYTE;
  }
  else if (output_data_type == DFNT_CHAR8) {
      nctype = NC_CHAR;
  }
  else if (output_data_type == DFNT_UINT16) {
      nctype = NC_USHORT;
  }
  else if (output_data_type == DFNT_INT16) {
      nctype = NC_SHORT;
  }
  else if (output_data_type == DFNT_UINT32) {
      nctype = NC_UINT;
  }
  else if (output_data_type == DFNT_INT32) {
      nctype = NC_INT;
  }
  else if (output_data_type == DFNT_FLOAT32) {
      nctype = NC_FLOAT;
  }
  else if (output_data_type == DFNT_FLOAT64) {
      nctype = NC_DOUBLE;
  }
  else {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Unrecognized output data type", "OutputFile", (Output_t *)NULL);      
  }
   
  retval = nc_def_var(ncid, this->sds.name, nctype, 2, dimids, &varid);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem creating variable", "OutputFile", (Output_t *)NULL);
      
  }
  else {
     this->sds.id = varid;
  }
  retval = nc_put_att_text(ncid, varid, GRID_MAPPING, strlen(GRID_PROJECTION), GRID_PROJECTION);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem putting attribute", "OutputFile", (Output_t *)NULL);      
  }
  
  /* Define spatial coordinate variables */
  xdimids[0] = dimids[1];
  ydimids[0] = dimids[0];
  
  retval = nc_def_var(ncid, this->sds.dim[1].name, NC_FLOAT, 1, xdimids, &x_varid);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem creating coordinate variable", "OutputFile", (Output_t *)NULL);
      
  }
  
  retval = nc_put_att_text(ncid, x_varid, STANDARD_NAME, strlen(PROJECTION_X_COORDINATE), PROJECTION_X_COORDINATE);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem putting attribute", "OutputFile", (Output_t *)NULL);      
  }
  
  retval = nc_put_att_text(ncid, x_varid, UNITS, strlen(METER), METER);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem putting units attribute", "OutputFile", (Output_t *)NULL);      
  }
  
  retval = nc_def_var(ncid, this->sds.dim[0].name, NC_FLOAT, 1, ydimids, &y_varid);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem creating coordinate variable", "OutputFile", (Output_t *)NULL);
      
  }
  
  retval = nc_put_att_text(ncid, y_varid, STANDARD_NAME, strlen(PROJECTION_Y_COORDINATE), PROJECTION_Y_COORDINATE);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem putting attribute", "OutputFile", (Output_t *)NULL);      
  }
  
  retval = nc_put_att_text(ncid, y_varid, UNITS, strlen(METER), METER);  
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem putting attribute", "OutputFile", (Output_t *)NULL);      
  }
  
  /* Define the Projection meta-data variable */
  retval = nc_def_var(ncid, GRID_PROJECTION, NC_INT, 0, dimid_1d, &proj_varid);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem creating projection variable", "OutputFile", (Output_t *)NULL);
  }

  float32 false_easting = 0.0;
  float32 false_northing = 0.0;
  float32 center_lon = space_def->orig_proj_param[4];
  float32 center_lat = space_def->orig_proj_param[5];
  
  retval = nc_put_att_float(ncid, proj_varid, FALSE_EASTING, NC_FLOAT, 1, &false_easting);
  retval = nc_put_att_float(ncid, proj_varid, FALSE_NORTHING, NC_FLOAT, 1, &false_northing);
  retval = nc_put_att_text(ncid, proj_varid, GRID_MAPPING_NAME, strlen(LAMAZ), LAMAZ);
  retval = nc_put_att_float(ncid, proj_varid, CENTER_LON, NC_FLOAT, 1, &center_lon);
  retval = nc_put_att_float(ncid, proj_varid, CENTER_LAT, NC_FLOAT, 1, &center_lat);
  
  
  
  /* Conventions global attribute */
  retval = nc_put_att_text(ncid, NC_GLOBAL, CONVENTIONS, strlen(CF_1_7), CF_1_7);
  
  /* End define mode. This tells netCDF we are done defining metadata */
  retval = nc_enddef(ncid);
  
  /* populate spatial coordinate variables x and y*/
  
  float32 *x_values = (float32 *) calloc(this->size.s, sizeof(float32));
  float32 *y_values = (float32 *) calloc(this->size.l, sizeof(float32));
  
  
  for (k=0; k<this->size.s; k++) {
      x_values[k] = space_def->ul_corner.x + (float32) (k*space_def->pixel_size);
  }
  
  for (k=0; k<this->size.l; k++) {
      y_values[k] = space_def->ul_corner.y - (float32) (k*space_def->pixel_size);
  }

  start[0] = 0;
  nval[0] = this->size.s;
  retval = nc_put_var_float(ncid, x_varid, x_values);
  if (retval != 0) {
      printf("retval: %d\n", retval);
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem writing values to coordinate variable", "OutputFile", (Output_t *)NULL);      
  }
  
  start[0] = 0;
  nval[0] = this->size.l;
  retval = nc_put_var_float(ncid, y_varid, y_values);
  if (retval != 0) {
      free(this->sds.name);
      free(this->file_name);
      free(this);
      LOG_RETURN_ERROR("Problem writing values to coordinate variable", "OutputFile", (Output_t *)NULL);      
  }

  if (error_string != (char *)NULL) {
    free(this->sds.name);
    free(this->file_name);
    free(this);  
    LOG_RETURN_ERROR(error_string, "OutputFile", 
                 (Output_t *)NULL); 
  }

  return this;
}


bool CloseOutput(Output_t *this)
/* 
!C******************************************************************************

!Description: 'CloseOutput' ends SDS access and closes the output file.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   open, sds.id, sds_file_id

!Output Parameters:
 this           'output' data structure; the following fields are modified:
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
   3. 'OutputFile' must be called before this routine is called.
   4. 'FreeOutput' should be called to deallocate memory used by the 
      'output' data structures.

!END****************************************************************************
*/
{
    int retval;

    if (!this->open) {
      LOG_RETURN_ERROR("file not open", "CloseOutput", false);
    }


    retval = nc_close(this->sds_file_id);
    this->open = false;

  return true;
}


bool FreeOutput(Output_t *this)
/* 
!C******************************************************************************

!Description: 'FreeOutput' frees the 'output' data structure memory.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   sds.rank, sds.dim[*].name, sds.name, file_name

!Output Parameters:
 this           'output' data structure; the following fields are modified:
                   sds.dim[*].name, sds.name, file_name
 (returns)      status:
                  'true' = okay (always returned)

!Team Unique Header:

 ! Design Notes:
   1. 'OutputFile' must be called before this routine is called.
   2. An error status is never returned.

!END****************************************************************************
*/
{
  int ir;

  if (this != (Output_t *)NULL) {
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

bool WriteOutput(Output_t *this, int iline, void *buf)
/* 
!C******************************************************************************

!Description: 'WriteOutput' writes a line of data to the output HDF file.
 
!Input Parameters:
 this           'output' data structure; the following fields are input:
                   open, size, sds.id
 iline          output line number
 buf            buffer of data to be written

!Output Parameters:
 this           'output' data structure; the following fields are modified:
 (returns)      status:
                  'true' = okay
		  'false' = error return

!Team Unique Header:

 ! Design Notes:
   1. An error status is returned when:
       a. the file is not open for access
       b. the line number is invalid (< 0; >= 'this->size.l')
       b. an error occurs when writting to the SDS.
   2. Error messages are handled with the 'LOG_RETURN_ERROR' macro.
   3. 'OutputFile' must be called before this routine is called.

!END****************************************************************************
*/
{
  size_t start[2], nval[2];
  int retval;

  /* Check the parameters */

  if (!this->open)
    LOG_RETURN_ERROR("file not open", "WriteOutput", false);

  if (iline < 0  ||  iline >= this->size.l)
    LOG_RETURN_ERROR("invalid line number", "WriteOutput", false);

  /* Write the data */

  start[0] = iline;
  start[1] = 0;
  nval[0] = 1;
  nval[1] = this->size.s;
  

  int32 type = this->sds.type;
  

  if (type == DFNT_INT8) {
      retval = nc_put_vara_schar(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_CHAR) {
      retval = nc_put_vara_uchar(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_UINT16) {
      retval = nc_put_vara_ushort(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_INT16) {
      retval = nc_put_vara_short(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_UINT32) {
      retval = nc_put_vara_uint(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_INT32) {
      retval = nc_put_vara_int(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_FLOAT32) {
      retval = nc_put_vara_float(this->sds_file_id, this->sds.id, start, nval, buf);
  }
  else if (type == DFNT_FLOAT64) {
      retval = nc_put_vara_double(this->sds_file_id, this->sds.id, start, nval, buf);      
  }
  else {
      LOG_RETURN_ERROR("unrecognized type", "WriteOutput", false);
  }
  
  if (retval != 0) {
      LOG_RETURN_ERROR("problem writing scan line", "WriteOutput", false);     
  }
  else {
      return true;
  }
}
