/*
!C****************************************************************************

!File: myhdf.c
  
!Description: Functions for handling HDF files.

!Revision History:
 Revision 1.0 2001/05/08
 Robert Wolfe
 Original Version.

 Revision 1.5 2002/11/02
 Gail Schmidt
 Added support for INT8 data types.

 Revision 2.0 2003/10/23
 Gail Schmidt
 Added routine to read the bounding coords from the metadata.
 Added support for multiple pixel sizes and number of lines/samples.

 Revision 2.0a 2004/09/23
 Gail Schmidt
 Modified the routine that reads the bounding coords to look for the
 bounding coords in the main metadata section, if it is not found in the
 Struct, Archive, or Core metadata structures.

!Team Unique Header:
  This software was developed by the MODIS Land Science Team Support 
  Group for the Labatory for Terrestrial Physics (Code 922) at the 
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
   1. The following functions handle Science Data Sets (SDSs) and attributes
      in HDF files:

       GetSDSInfo - Read SDS information.
       GetSDSDimInfo - Read SDS dimension information.
       PutSDSInfo - Create an SDS and write information.
       PutSDSDimInfo - Write SDS dimension information.
       GetAttrDouble - Get an HDF attribute's value.
       ReadBoundCoords - Read the bounding coordinates from the metadata.
       ReadMetadata - Read the specified attribute from the metadata.
       DetermineResolution - Determine the resolution of the specified SDSs.
       DeterminePixelSize - Determine the pixel size of the specified SDSs.

!END****************************************************************************
*/

#include <stdlib.h>
#include <math.h>
#include "myhdf.h"
#include "myerror.h"
#include "hdf.h"
#include "mfhdf.h"
#include "mystring.h"
#include "input.h"
#include "geoloc.h"
#include "myproj.h"
#include "const.h"
#include "hdf5.h"

/* Constants */

#define DIM_MAX_NCHAR (80)  /* Maximum size of a dimension name */


hid_t       file, datasetHandle;   /* handles */
hid_t       datatype, dataspace;
hid_t       memspace;
herr_t       status_n;
H5T_class_t t_class;               /* data type class */
H5T_order_t order;                 /* data order */
size_t      size;                  /*
                                    * size of the data element
                                    * stored in file
                                    */ 
hsize_t     datasetDimensions[2];

bool GetSDSInfoV(int32 sds_file_id, Myhdf_sds_t *sds) {
    sds->id = H5Dopen(sds_file_id, sds->name, H5P_DEFAULT);
    hid_t datasetHandle = (hid_t) sds->id;
    
    /*
     * Get datatype and dataspace handles and then query
     * dataset class, order, size, rank and dimensions.
     */
    datatype  = H5Dget_type(datasetHandle); /* datatype handle */
    t_class     = H5Tget_class(datatype);
    order     = H5Tget_order(datatype);

    size  = H5Tget_size(datatype);

    dataspace = H5Dget_space(datasetHandle); /* dataspace handle */
    int32 rank      = H5Sget_simple_extent_ndims(dataspace);
    long *datasetDimensions = (long *) calloc(rank, sizeof(long));
    status_n  = H5Sget_simple_extent_dims(dataspace, datasetDimensions, NULL);

    sds->rank = rank;
    sds->type = datatype;
    sds->typeh5 = datatype;
    sds->datasize = size;
    sds->nattr = 0;
   
    return true;
}

bool GetSDSDimInfoV(int32 sds_id, Myhdf_dim_t *dim, int irank) {
    hid_t datasetHandle = (hid_t) sds_id;
    dataspace = H5Dget_space(datasetHandle);
    int32 rank      = H5Sget_simple_extent_ndims(dataspace);
    long *datasetDimensions = (long *) calloc(rank, sizeof(long));
    status_n  = H5Sget_simple_extent_dims(dataspace, datasetDimensions, NULL);

    dim->nval = datasetDimensions[irank];
    dim->name = DupString("dim name unspecified");

    return true;    
}

bool GetSingleValueAttrAsDouble(hid_t dataset, char *attrName, double *dblVal) {
    hid_t attr, attrType, nativeType;
    herr_t status;
    
    unsigned int uintValue;
    int intValue;
    unsigned short ushortValue;
    short shortValue;
    float floatValue;
    double doubleValue;
    signed char sbyteValue;
    unsigned char ubyteValue;
    char msg[256];
    
    if (!(H5Aexists(dataset, attrName))) {
        sprintf(msg, "no attribute %s \n", attrName);
	LogInfomsg(msg);
        return false;
    }
    
    attr = H5Aopen(dataset, attrName, H5P_DEFAULT);
    attrType = H5Aget_type (attr);
    nativeType = getNativeType(attrType);
    if (nativeType == H5I_INVALID_HID) {
        LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
    }
    
    if (nativeType == H5T_NATIVE_UINT) {
        status = H5Aread(attr, nativeType, &uintValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        dblVal[0] = (double) uintValue;
    }
    else if (nativeType == H5T_NATIVE_INT) {
        status = H5Aread(attr, nativeType, &intValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        dblVal[0] = (double) intValue;        
    }
    else if (nativeType == H5T_NATIVE_USHORT) {
        status = H5Aread(attr, nativeType, &ushortValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        dblVal[0] = (double) ushortValue;
    }
    else if (nativeType == H5T_NATIVE_SHORT) {
        status = H5Aread(attr, nativeType, &shortValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        dblVal[0] = (double) shortValue;
    }    
    else if (nativeType == H5T_NATIVE_FLOAT) {
        status = H5Aread(attr, nativeType, &floatValue);
        if (status < 0) {
             LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);           
        }
        dblVal[0] = (double) floatValue;
    }
    else if (nativeType == H5T_NATIVE_DOUBLE) {
        status = H5Aread(attr, nativeType, &doubleValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);            
        }
        dblVal[0] = doubleValue;
    }
    else if (nativeType == H5T_NATIVE_SCHAR) {
        status = H5Aread(attr, nativeType, &sbyteValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);            
        }
        dblVal[0] = (double) sbyteValue;
    }
    else if (nativeType == H5T_NATIVE_UCHAR) {
        status = H5Aread(attr, nativeType, &ubyteValue);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);            
        }
        dblVal[0] = (double) ubyteValue;
    }
                

    return true;
}

bool GetAttrArrayAsDouble(hid_t dataset, char *attrName, int *arrayLen, double **dblArray)
{
    hid_t attr, attrType, nativeType;
    herr_t status;
    H5A_info_t ainfo;
    
    int i;
    char msg[256];
    size_t size;
    double *dblVal;
    
    void *attrValues;
    
    if (!(H5Aexists(dataset, attrName))) {
        snprintf(msg, 256, "no attribute %s \n", attrName);
	LogInfomsg(msg);
        return false;
    }
       
    attr = H5Aopen(dataset, attrName, H5P_DEFAULT);
    
    status = H5Aget_info(attr, &ainfo);
    if (status < 0) {
        LOG_RETURN_ERROR(" ", "GetArrayAttrAsDouble", false);
    }
    attrValues = malloc(ainfo.data_size);
    
    attrType = H5Aget_type (attr);
    size = H5Tget_size(attrType);
    int numVals = ainfo.data_size/size;
    arrayLen[0] = numVals;
    dblVal = (double *) malloc(ainfo.data_size);
    
    nativeType = getNativeType(attrType);
    
    if (nativeType == H5I_INVALID_HID) {
        LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
    }
    
    if (nativeType == H5T_NATIVE_UINT) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((unsigned int *)attrValues)[i];
        }
    }
    else if (nativeType == H5T_NATIVE_INT) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((int *)attrValues)[i];
        }        
    }
    else if (nativeType == H5T_NATIVE_USHORT) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((unsigned short *)attrValues)[i];
        }
    }
    else if (nativeType == H5T_NATIVE_SHORT) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((short *)attrValues)[i];
        }        
    }    
    else if (nativeType == H5T_NATIVE_FLOAT) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
             LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);           
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((float *)attrValues)[i];
        }        
    }
    else if (nativeType == H5T_NATIVE_DOUBLE) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);            
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((double *)attrValues)[i];
        }        
    }
    else if (nativeType == H5T_NATIVE_SCHAR) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);            
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((signed char *)attrValues)[i];
        }        
    }
    else if (nativeType == H5T_NATIVE_UCHAR) {
        status = H5Aread(attr, nativeType, attrValues);
        if (status < 0) {
            LOG_RETURN_ERROR(" ", "GetSingleValueAttrAsDouble", false);            
        }
        for (i=0; i<numVals; i++) {
            dblVal[i] = (double) ((unsigned char *)attrValues)[i];
        }        
    }
    
    *dblArray = dblVal;
      
    free(attrValues);

    return true;    
}

bool DetermineResolution(Myhdf_sds_t *sds, Img_coord_int_t *ls_dim, int *ires)
/*
!C*****************************************************************************
!Description: 'DetermineResolution' determines the resolution of the SDS,
 based on the nominal 1KM frame size.

!Input Parameters:
 sds            SDS structure. The dimension values are used to determine
                the number of lines in the SDS.
 ls_dim         Dimension locations for the line and sample dimensions.

!Output Parameters:
 ires           Resolution value of the SDS
 (returns)      Status:
                  'true' = okay
                  'false' = error reading the bounding coords

!Team Unique Header:

! Design Notes:

!END*****************************************************************************/
{
  Img_coord_int_t size;                /* input file size */

  /* Check the line and sample dimensions */
  size.s = sds->dim[ls_dim->s].nval;
  *ires = -1;
  *ires = (int)((size.s / (double)NFRAME_1KM_MODIS) + 0.5);
  
  /* TDR, VIIRS always = 1, hardcode for now*/
  *ires = 1;

  /* Verify the resolution. If not 250m, 500m, or 1km product, then this
     is an error. */
  if (*ires != 1 && *ires != 2 && *ires != 4) {
    LOG_RETURN_ERROR("invalid resolution", "DetermineResolution", false);
  }

  return true;
}


bool DeterminePixelSizeV(char *geoloc_file_name, int num_input_sds, 
  char *geoProductName, int out_proj_num,
  double output_pixel_size[MAX_SDS_DIMS])
/*
!C*****************************************************************************
!Description: 'DeterminePixelSize' uses the ires values (calculated by
 DetermineResolution) to determine the output pixel size.

!Input Parameters:
 Param_t        Input parameter structure. The output_pixel_size values
                are updated for each SDS specified.

!Output Parameters:
 (returns)      Status:
                  'true' = okay
                  'false' = error reading the bounding coords

!Team Unique Header:

! Design Notes:
  1. If the output projection is Geographic, then the pixel size will need
     to be computed in degrees.  The algorithm used in this routine uses the
     nominal pixel size to determine the pixel size in meters.  For Geographic,
     the difference between two lat/long locations will be used to determine
     the pixel size in degrees.  That value will be modified depending on
     whether the resolution is 250m, 500m, or 1km (no modification necessary).

!END*****************************************************************************/
{
  int i;
  Geoloc_t *geoloc = NULL;             /* geolocation file */
  int center_loc;                      /* center location in the scan */
  int midscan;                         /* middle scan location */
  int32 start[MYHDF_MAX_RANK];         /* start reading at this location */
  int32 nval[MYHDF_MAX_RANK];          /* read this many values */
  double center = 0.0;                 /* longitude value at the center */
  double centerp1 = 0.0;               /* longitude value at the center+1 */

  /* If dealing with an output projection of Geographic, then read the
     center pixel values from the Geolocation file for use later in the
     pixel size calculation. */
  if (out_proj_num == PROJ_GEO)
  {
      /* Open geoloc file */
      geoloc = OpenGeolocSwath(geoloc_file_name);
      if (geoloc == (Geoloc_t *)NULL)
        LOG_RETURN_ERROR("bad geolocation file", "DeterminePixelSize",
                              false);

      /* Grab the line representing the center of the swath (nscan/2) */
      midscan = geoloc->nscan / 2;
      start[0] = midscan * geoloc->scan_size.l;
      start[1] = 0;
      nval[0] = 1;
      nval[1] = geoloc->scan_size.s;

      if (readData(geoloc->sds_lon.id, start, nval, geoloc->lon_buf) == FALSE) {
         LOG_RETURN_ERROR("reading longitude", "DeterminePixelSize", false);          
      }

      /* Get the center of the scan and the center of the scan plus one */
      center_loc = geoloc->scan_size.s / 2;
      center = geoloc->lon_buf[center_loc];
      centerp1 = geoloc->lon_buf[center_loc + 1];
      
      /* Close geolocation file */
      if (!CloseGeoloc(geoloc))
        LOG_RETURN_ERROR("closing geolocation file", "DeterminePixelSize",
                              false);
  }
  
  double pixelSize = getVIIRSpixelResolutionFromGeoProductName(geoProductName);
  if (pixelSize == 0) {
        LOG_RETURN_ERROR(strcat("can't determine pixel size for product: ", geoProductName), 
                "DeterminePixelSize", false);
  }

  /* Loop through all the SDSs to be processed */
  for (i = 0; i < num_input_sds; i++)
  {
    /* Determine the output pixel size */
    /* If the output projection is not Geographic, then the output pixel
       size is in meters */
    if (out_proj_num != PROJ_GEO)
    {
        output_pixel_size[i] = (double) pixelSize;
    }
    /* If the output projection is Geographic, then we need to read the
       Geolocation file */
    else
    {
       double posDiff = fabs(centerp1 - center);
       if (posDiff > 180) {
           output_pixel_size[i] = 360 - posDiff;
       }
       else {
           output_pixel_size[i] = posDiff;
       }
    }
  }

  if (out_proj_num == PROJ_GEO)
  {
    /* Free geolocation structure */
    if (!FreeGeoloc(geoloc))
      LOG_RETURN_ERROR("freeing geoloc file struct", "DeterminePixelSize",
      false);
  }

  return true;
}

/*
!C*****************************************************************************
!Description: 'readData' returns a hyperslab from dataset integer coordinates
 in a 1D array.

!Input Parameters:
 hid_t dataset identifier.
 int *start, int *count: data subset
 void *buf: buffer allocated by the caller to hold subset values.

!Output Parameters:
 (returns)      Status:
                  'true' = okay
                  'false' = error reading data

! Design Notes:
  1. Only handles 2D datasets. TODO: generalize.

!END*****************************************************************************/
bool readData(hid_t dataset, int *start, int *count, void *buf) {
    hid_t       datatype, outtype, dataspace;
    hid_t       memspace;
    int         status_n;
    int         rank;
    hsize_t     dims2D[2];
    hsize_t     hyperSlabStart2D[2];
    hsize_t     hyperSlabCount2D[2];
    
    
    char msg[M_MSG_LEN+1];
    
    if (dataset < 0) {
        sprintf(msg, "problem retrieving datatype from dataset %d \n", dataset);
        LOG_RETURN_ERROR(msg, "readData", false);        
    }
 
    
    datatype  = H5Dget_type(dataset);
    if (datatype < 0) {
        sprintf(msg, "problem retrieving datatype from dataset %d \n", dataset);
        LOG_RETURN_ERROR(msg, "readData", false);
        
    }
    
    outtype = getNativeType(datatype);

    dataspace = H5Dget_space(dataset);
    if (dataspace < 0) {
        sprintf(msg, "problem retrieving dataspace from dataset %d \n", dataset);
        LOG_RETURN_ERROR(msg, "readData", false);
    }
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims2D, NULL);   
    
    hyperSlabStart2D[0] = start[0];
    hyperSlabStart2D[1] = start[1];
    
    hyperSlabCount2D[0] = count[0];
    hyperSlabCount2D[1] = count[1];
    
    
    status_n = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hyperSlabStart2D, NULL, hyperSlabCount2D, NULL);
    if (status_n < 0) {
        sprintf(msg, "problem retrieving hyperslab from dataset %d \n", dataset);
        LOG_RETURN_ERROR(msg, "readData", false);               
    }
    if (!H5Sselect_valid(dataspace)) {
        sprintf(msg, "hyperslab selection not valid \n");
        LOG_RETURN_ERROR(msg, "readData", false);
    }
    
    
    dims2D[0] = count[0];
    dims2D[1] = count[1];
    memspace = H5Screate_simple(rank, dims2D, NULL);
    
    status_n = H5Dread(dataset, outtype, memspace, dataspace, H5P_DEFAULT, buf);
    if (status_n < 0) {
        sprintf(msg, "problem reading from dataset %d \n", dataset);
        LOG_RETURN_ERROR(msg, "readData", false);   
    }
    
    /*
     * Close/release resources.
     */
    H5Tclose(datatype);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    
    return true;
}

/*
!C*****************************************************************************
!Description: 'getNativeType' maps dataset datatype to a native datatype

!Input Parameters:
 hid_t datatype: dataset datatype

!Output Parameters:
 (returns)  hid_t native datatype

!END*****************************************************************************/
hid_t getNativeType(hid_t datatype) {
    hid_t outtype;
    H5T_class_t t_class;
    size_t      size;
    H5T_sign_t  sign;
    
    char msg[M_MSG_LEN+1];
    
    t_class   = H5Tget_class(datatype);
    size      = H5Tget_size(datatype);
    sign      = H5Tget_sign(datatype);

    if (t_class == H5T_INTEGER) {
        if (size == 1) {
            outtype = H5T_NATIVE_SCHAR;
            if (sign == H5T_SGN_NONE) {
                outtype = H5T_NATIVE_UCHAR;
            }
        }
        else if (size == 2) {
            outtype = H5T_NATIVE_SHORT;
            if (sign == H5T_SGN_NONE) {
                outtype = H5T_NATIVE_USHORT;
            }
        }
        else if (size == 4) {
            outtype = H5T_NATIVE_INT;
            if (sign == H5T_SGN_NONE) {
                outtype = H5T_NATIVE_UINT;
            }            
        }
        else if (size == 8) {
            outtype = H5T_NATIVE_LONG;
            if (sign == H5T_SGN_NONE) {
                outtype = H5T_NATIVE_ULONG;
            }            
        }
        else {
            sprintf(msg, "can only handle INTEGER sizes of 1, 2, 4, 8, not: %d \n", (int) size);
            LogInfomsg(msg);
            return H5I_INVALID_HID;            
        }
    }
    else if (t_class == H5T_FLOAT) {
        if (size == 4) {
            outtype = H5T_NATIVE_FLOAT;
        }
        else if (size == 8) {
            outtype = H5T_NATIVE_DOUBLE;
        }
        else {
            sprintf(msg, "can't handle dataset FLOAT type with size: %d \n", (int) size);
            LogInfomsg(msg);
            return H5I_INVALID_HID;
        }
    }
    else {
        sprintf(msg, "can only handle numeric (FLOAT, INTEGER) datatypes, not: %d \n", (int) t_class);
        LogInfomsg(msg);
        return H5I_INVALID_HID;
    }
    
    return outtype;
}