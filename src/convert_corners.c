/* 
!C******************************************************************************

!Description: 'ConvertCorners' uses the output UL and LR lat/long corners to
 determine the UL corner in output space and the number of lines and samples
 in output space (if output spatial subset type is LAT_LONG). If the output
 spatial subset type is PROJ_COORDS, then the number of lines and samples
 are calculated from the UL and LR proj coords. If the output spatial subset
 type is LINE_SAMPLE, then the UL and LR corners are calculated given the
 UL and LR line/sample values in input space.
 
!Input Parameters:
 param             list of user parameters; the following fields are input:
                   output_spatial_subset_type

 output_space_def  output grid space definition; the following fields are
                   input: pixel_size, ul_corner.lat, ul_corner.lon,
                   lr_corner.lat, lr_corner.lon, proj_num, zone,
                   sphere, proj_param[*]
                   the following fields are output: ul_corner.x, ul_corner.y,
                   lr_corner.x, lr_corner.y, ul_corner_geo.x, ul_corner_geo.y,
                   lr_corner_geo.x, lr_corner_geo.y, img_size.l, img_size.s

!Output Parameters:
 (returns)      an integer value (true, false) specifying no error or an error

!Developers:
 Gail Schmidt
 SAIC / USGS EROS Data Center
 Rapid City, SD 57701
 gschmidt@usgs.gov
    
!Notes:

!END****************************************************************************
*/

#if !defined(__CYGWIN__) && !defined(__APPLE__) && !defined(WIN32)
#include <values.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "param.h"

#include "geoloc.h"
#include "parser.h"
#include "usage.h"

#include "myproj.h"

#include "myerror.h"
#include "mydtype.h"
#include "deg2dms.h"
#include "cproj.h"
#include "const.h"
#include <float.h>
#include "space.h"

#define DEBUG

/* Functions */

int ConvertCorners(Param_t *param)
{
  int i;
  int ul_line, ul_samp, lr_line, lr_samp;  /* UL/LR line/sample values */
  Map_coord_t ul, lr;
  double minx = FLT_MAX,
         maxx = -FLT_MAX,
         miny = FLT_MAX,
         maxy = -FLT_MAX;
  Space_t *out_space = NULL;
  Space_def_t *output_space_def = NULL;
  Geo_coord_t geo_coord_p;
  Map_coord_t output_coord_p;
  Geoloc_t *geoloc = NULL;             /* geolocation file */
  
  int32 start[2];
  int32 nval[2];
  
  
  Map_coord_t output_coord_UL;
  Map_coord_t output_coord_LR;
  Map_coord_t output_coord_MP;
  
  
  float *lon_side_a;
  float *lat_side_a;
  float *lon_side_b;
  float *lat_side_b;
  float *lon_side_c;
  float *lat_side_c;
  float *lon_side_d;
  float *lat_side_d;
  int numSidesSpan = 0;
  
  int userDefinedProjParams = FALSE;
  
  /* Make the pointers cleaner */
  output_space_def = &(param->output_space_def);
  output_space_def->straddlesDateline = false;
  output_space_def->containsPole = false;
  
  for (i=0; i<NPROJ_PARAM; i++) {
      if (output_space_def->proj_param[i] != -999.0) {
          userDefinedProjParams = TRUE;
          break;
      }
  }

  /* If the output spatial subset type is LINE_SAMPLE, then the UL and LR
     corner points have been provided as line/sample (in input space). First
     determine the lat/long UL/LR, then treat this as a LAT_LONG spatial
     subset type and determine the minimum bounding rectangle. */
  if (param->output_spatial_subset_type == LINE_SAMPLE)
  {
    ul_samp = (int)output_space_def->ul_corner.x;
    ul_line = (int)output_space_def->ul_corner.y;
    lr_samp = (int)output_space_def->lr_corner.x;
    lr_line = (int)output_space_def->lr_corner.y;
    
    int del_line = lr_line - ul_line;
    int del_samp = lr_samp - ul_samp;
    
    geoloc = OpenGeolocSwath(param, (Input_t *)NULL);
    
    if (geoloc == (Geoloc_t *)NULL)
      LOG_RETURN_ERROR("bad geolocation file", "ConvertCorners", false);

    
    /* Grab the UL to UR longitude (0-based start values) */
    start[0] = ul_line;
    start[1] = ul_samp;
    nval[0] = 1;
    nval[1] = del_samp + 1;
    lon_side_a = calloc(nval[1], sizeof(float32));
    if (!readData(geoloc->sds_lon.id, start, nval, lon_side_a)) {
        LOG_RETURN_ERROR("reading longitude", "ConvertCorners", false);          
    }    

    /* Grab the UL to UR latitude (0-based start values) */
    start[0] = ul_line;
    start[1] = ul_samp;
    nval[0] = 1;
    nval[1] = del_samp + 1;
    lat_side_a = calloc(nval[1], sizeof(float32));
    if (!readData(geoloc->sds_lat.id, start, nval, lat_side_a)) {
       LOG_RETURN_ERROR("reading latitude", "ConvertCorners", false);          
    }

    /* Grab the UR to LR longitude (0-based start values) */
    start[0] = ul_line;
    start[1] = ul_samp + del_samp;
    nval[0] = del_line + 1;
    nval[1] = 1;
    lon_side_b = calloc(nval[0], sizeof(float32));
    if (!readData(geoloc->sds_lon.id, start, nval, lon_side_b)) {
       LOG_RETURN_ERROR("reading longitude", "DeterminePixelSize", false);          
    }   

    /* Grab the UR to LR latitude (0-based start values) */
    start[0] = ul_line;
    start[1] = ul_samp + del_samp;
    nval[0] = del_line + 1;
    nval[1] = 1;
    lat_side_b = calloc(nval[0], sizeof(float32));
    if (!readData(geoloc->sds_lat.id, start, nval, lat_side_b)) {
       LOG_RETURN_ERROR("reading latitude", "ConvertCorners", false);          
    }    
     
    start[0] = ul_line + del_line;
    start[1] = ul_samp;
    nval[0] = 1;
    nval[1] = del_samp + 1;
    lon_side_c = calloc(nval[1], sizeof(float32));
    if (!readData(geoloc->sds_lon.id, start, nval, lon_side_c)) {
       LOG_RETURN_ERROR("reading longitude", "ConvertCorners", false);          
    }

    start[0] = ul_line + del_line;
    start[1] = ul_samp;
    nval[0] = 1;
    nval[1] = del_samp + 1;
    lat_side_c = calloc(nval[1], sizeof(float32));
    if (!readData(geoloc->sds_lat.id, start, nval, lat_side_c)) {
       LOG_RETURN_ERROR("reading latitude", "ConvertCorners", false);          
    }
    
    start[0] = ul_line;
    start[1] = ul_samp;
    nval[0] = del_line + 1;
    nval[1] = 1;
    lon_side_d = calloc(nval[0], sizeof(float32));
    if (!readData(geoloc->sds_lon.id, start, nval, lon_side_d)) {
       LOG_RETURN_ERROR("reading longitude", "ConvertCorners", false);          
    }

    start[0] = ul_line;
    start[1] = ul_samp;
    nval[0] = del_line + 1;
    nval[1] = 1;
    lat_side_d = calloc(nval[0], sizeof(float32));
    if (!readData(geoloc->sds_lat.id, start, nval, lat_side_d)) {
       LOG_RETURN_ERROR("reading latitude", "ConvertCorners", false);          
    }
    
    float32 cntrLon[1];
    float32 cntrLat[1];
    if (userDefinedProjParams == TRUE) {
    }
    else {
        start[0] = ul_line + del_line/2;
        start[1] = ul_samp + del_samp/2;
        nval[0] = 1;
        nval[1] = 1;
        if (!readData(geoloc->sds_lon.id, start, nval, cntrLon)) {
            LOG_RETURN_ERROR("reading longitude", "ConvertCorners", false);          
        }  
        if (!readData(geoloc->sds_lat.id, start, nval, cntrLat)) {
            LOG_RETURN_ERROR("reading latitude", "ConvertCorners", false);          
        }
        
        if (output_space_def->proj_num == 9)  { // Tranverse Mercator
            output_space_def->proj_param[2] = 1.0;
        }
        else {
            output_space_def->proj_param[2] = cntrLat[0];
        }
        output_space_def->proj_param[3] = cntrLat[0];
        output_space_def->proj_param[4] = cntrLon[0];
        output_space_def->proj_param[5] = cntrLat[0];
        output_space_def->proj_param[6] = 0.0;
        output_space_def->proj_param[7] = 0.0;
        if (output_space_def->proj_num == 31) { // Integerized Sinusoidal
            output_space_def->proj_param[8] = 2.0;
            output_space_def->proj_param[10] = 0.0;
        }
    }
    

    /* Copy the projection parameters to orig_proj_param to use the decimal
       degree values later (GeoTiff output) */
    for (i = 0; i < NPROJ_PARAM; i++) {
         param->output_space_def.orig_proj_param[i] =
         param->output_space_def.proj_param[i];
    }

    /* Convert the output projection parameter lat/long values from decimal
       degrees to DMS */
    if (!Deg2DMS (param->output_space_def.proj_num,
                 param->output_space_def.proj_param)) {
        FreeParam(param);
        LOG_RETURN_ERROR("error converting projection parameters from decimal degrees to DMS",
                         "ConvertCorners", false);
    }
    
    output_space_def->pixel_size = param->output_pixel_size[0];
    output_space_def->img_size.s = 1;
    output_space_def->img_size.l = 1;
    
    /* Get the forward and reverse transformation functions */
    out_space = SetupSpace(output_space_def);
    if (out_space == (Space_t *)NULL) {
         LOG_RETURN_ERROR("setting up output space", "ConvertCorners", false); 
    }
    
    minx = FLT_MAX,
    maxx = -FLT_MAX,
    miny = FLT_MAX,
    maxy = -FLT_MAX;
    
    bool spansDateline = false;
    float lon, lat;
    bool pos;
    
    //int numSidesSpan = 0;
    
    if (output_space_def->proj_num == PROJ_GEO)
    {
        for (i=0; i<del_samp-1; i++) {
            if (fabs(lon_side_a[i+1] - lon_side_a[i]) > 180) {
                numSidesSpan++;
                break;
            }
        }

        for (i=0; i<del_line-1; i++) {
            if (fabs(lon_side_b[i+1] - lon_side_b[i]) > 180) {
                numSidesSpan++;
                break;
            }            
        }

        for (i=0; i<del_samp-1; i++) {
            if (fabs(lon_side_c[i+1] - lon_side_c[i]) > 180) {
                numSidesSpan++;
                break;
            }
        }

        for (i=0; i<del_line-1; i++) {
            if (fabs(lon_side_d[i+1] - lon_side_d[i]) > 180) {
                numSidesSpan++;
                break;
            }
        }
        
        minx = 180;
        maxx = -180;
        miny = 90;
        maxy = -90;
        
        for (i=0; i<del_samp; i++) {
            lat = lat_side_a[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
            
            lon = lon_side_a[i];
            if (lon < minx) minx = lon;
            if (lon > maxx) maxx = lon;            
        }
        
        for (i=0; i<del_line; i++) {
            lat = lat_side_b[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
            
            lon = lon_side_b[i];
            if (lon < minx) minx = lon;
            if (lon > maxx) maxx = lon;                        
        }  
        
        for (i=0; i<del_samp; i++) {
            lat = lat_side_c[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
            
            lon = lon_side_c[i];
            if (lon < minx) minx = lon;
            if (lon > maxx) maxx = lon;                        
        }  
        
        for (i=0; i<del_line; i++) {
            lat = lat_side_d[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
            
            lon = lon_side_d[i];
            if (lon < minx) minx = lon;
            if (lon > maxx) maxx = lon;                        
        }
    
    if (numSidesSpan == 2) {
        output_space_def->straddlesDateline = true;
        float lon_west = 180;
        float lon_east = -180;
        
        for (i=0; i<del_samp; i++) {
            lon = lon_side_a[i];
            if (lon > 0) {
                if (lon < lon_west) lon_west = lon;
            }
            else {
                if (lon > lon_east) lon_east = lon;
            }
        }
        
        for (i=0; i<del_line; i++) {
            lon = lon_side_b[i];
            if (lon > 0) {
                if (lon < lon_west) lon_west = lon;
            }
            else {
                if (lon > lon_east) lon_east = lon;
            }
            
        }  
        
        for (i=0; i<del_samp; i++) {
            lon = lon_side_c[i];
            if (lon > 0) {
                if (lon < lon_west) lon_west = lon;
            }
            else {
                if (lon > lon_east) lon_east = lon;
            }            
        }  
        
        for (i=0; i<del_line; i++) {
            lon = lon_side_d[i];
            if (lon > 0) {
                if (lon < lon_west) lon_west = lon;
            }
            else {
                if (lon > lon_east) lon_east = lon;
            }            
        }
        minx = lon_west;
        maxx = lon_east;
    }
    else if (numSidesSpan == 1) { // must contain north or south pole
        output_space_def->containsPole = true;
        minx = -180;
        maxx = 180;
        miny = 90;
        maxy = -90;
        
        for (i=0; i<del_samp; i++) {
            lat = lat_side_a[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
        }
        
        for (i=0; i<del_line; i++) {
            lat = lat_side_b[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
            
        }  
        
        for (i=0; i<del_samp; i++) {
            lat = lat_side_c[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
        }  
        
        for (i=0; i<del_line; i++) {
            lat = lat_side_d[i];
            if (lat < miny) miny = lat;
            if (lat > maxy) maxy = lat;
        }
    }
        minx *= RAD;
        miny *= RAD;
        maxx *= RAD;
        maxy *= RAD;
    }
    else { // Not GEO, everything else.
        for (i=0; i<del_samp; i++) {
            geo_coord_p.lon = lon_side_a[i]*RAD;
            geo_coord_p.lat = lat_side_a[i]*RAD;
            geo_coord_p.is_fill = false;
            if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
                 &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
                 FreeSpace(out_space);
                 LOG_RETURN_ERROR("converting UL to output map coordinates",
                              "ConvertCorners", false);
            }
            if (output_coord_p.x < minx) minx = output_coord_p.x;
            if (output_coord_p.x > maxx) maxx = output_coord_p.x;
            if (output_coord_p.y < miny) miny = output_coord_p.y;
            if (output_coord_p.y > maxy) maxy = output_coord_p.y;        
        }
        for (i=0; i<del_line; i++) {
            geo_coord_p.lon = lon_side_b[i]*RAD;
            geo_coord_p.lat = lat_side_b[i]*RAD;
            geo_coord_p.is_fill = false;
            if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
                 &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
                 FreeSpace(out_space);
                 LOG_RETURN_ERROR("converting UL to output map coordinates",
                              "ConvertCorners", false);
            }
            if (output_coord_p.x < minx) minx = output_coord_p.x;
            if (output_coord_p.x > maxx) maxx = output_coord_p.x;
            if (output_coord_p.y < miny) miny = output_coord_p.y;
            if (output_coord_p.y > maxy) maxy = output_coord_p.y;        
        }
        for (i=0; i<del_samp; i++) {
            geo_coord_p.lon = lon_side_c[i]*RAD;
            geo_coord_p.lat = lat_side_c[i]*RAD;
            geo_coord_p.is_fill = false;
            if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
                 &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
                 FreeSpace(out_space);
                 LOG_RETURN_ERROR("converting UL to output map coordinates",
                              "ConvertCorners", false);
            }
            if (output_coord_p.x < minx) minx = output_coord_p.x;
            if (output_coord_p.x > maxx) maxx = output_coord_p.x;
            if (output_coord_p.y < miny) miny = output_coord_p.y;
            if (output_coord_p.y > maxy) maxy = output_coord_p.y;        
        }
        for (i=0; i<del_line; i++) {
            geo_coord_p.lon = lon_side_d[i]*RAD;
            geo_coord_p.lat = lat_side_d[i]*RAD;
            geo_coord_p.is_fill = false;
            if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
                 &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
                 FreeSpace(out_space);
                 LOG_RETURN_ERROR("converting UL to output map coordinates",
                              "ConvertCorners", false);
            }
            if (output_coord_p.x < minx) minx = output_coord_p.x;
            if (output_coord_p.x > maxx) maxx = output_coord_p.x;
            if (output_coord_p.y < miny) miny = output_coord_p.y;
            if (output_coord_p.y > maxy) maxy = output_coord_p.y;        
        }
    }

     for (i = 0; i < param->num_input_sds; i++)
     {
       if (output_space_def->proj_num == PROJ_GEO)
       {
         /* Convert Geographic pixel size from degrees to radians */
           
           if (output_space_def->containsPole) {
               param->output_pixel_size[i] = (maxx - minx)/(lr_line - ul_line);
           }
           else {
               param->output_pixel_size[i] *= RAD;
           }
       }

       /* Calculate the number of output lines and samples */
       param->output_img_size[i].s =
         (int) ((maxx - minx) / param->output_pixel_size[i] + 0.5);
       param->output_img_size[i].l =
         (int) ((maxy - miny) / param->output_pixel_size[i] + 0.5);
       
       if (output_space_def->proj_num == PROJ_GEO) {
         
           float del;
           
           if (minx < 0 && maxx < 0) {
               del = abs(minx*DEG - maxx*DEG);
           }
           else if (minx > 0 && maxx > 0) {
               del = abs(minx*DEG - maxx*DEG);
           }
           else if (minx > 0 && maxx < 0) {
               del = (180 - minx*DEG) - (-180 - maxx*DEG);
           }
           else if (minx < 0 && maxx > 0) {
               del = maxx*DEG - minx*DEG;
           }

           param->output_img_size[i].s =
             (int) ((del*RAD) / param->output_pixel_size[i] + 0.5);           

           param->output_img_size[i].l =
             (int) ((maxy - miny) / param->output_pixel_size[i] + 0.5);           
           }

     }

     /* Use the first pixel size and number of lines/samples for the
        rest of the calculations */
     output_space_def->pixel_size = param->output_pixel_size[0];
     output_space_def->img_size.l = param->output_img_size[0].l;
     output_space_def->img_size.s = param->output_img_size[0].s;

     /* Redefine the output coords to use the minimum bounding box in
        output projection coords. Also make the LR corner a factor of the
        number of lines and samples. */
     output_space_def->ul_corner.x = minx;
     output_space_def->ul_corner.y = maxy;
     
     if (output_space_def->proj_num == PROJ_GEO) {
         output_space_def->lr_corner.x = maxx;
         output_space_def->lr_corner.y = miny;
     }
     else {
         output_space_def->lr_corner.x = output_space_def->ul_corner.x +
             output_space_def->img_size.s * output_space_def->pixel_size;
     
         output_space_def->lr_corner.y = output_space_def->ul_corner.y -
             output_space_def->img_size.l * output_space_def->pixel_size;
     }
     
     
    /* Close geolocation file */
    if (!CloseGeoloc(geoloc)) {
      FreeGeoloc(geoloc);
      LOG_RETURN_ERROR("closing geolocation file", "ConvertCorners", false);
    }

    /* Free geolocation structure */
    if (!FreeGeoloc(geoloc)) {
      LOG_RETURN_ERROR("freeing geoloc file struct", "ConvertCorners", false);
    } 
    
    free(lon_side_a);
    free(lat_side_a);
    free(lon_side_b);
    free(lat_side_b);
    free(lon_side_c);
    free(lat_side_c);
    free(lon_side_d);
    free(lat_side_d);
  } /* if LINE_SAMPLE */

  /* If the output spatial subset type is LAT_LONG, use the lat/long corner
     points to get the UL corner in output space and the number of
     lines/samples in the output image. If the output spatial subset type
     is LINE_SAMPLE, then the corner points have been converted to lat/long.
     So handle them as LAT_LONG. */
  else if (param->output_spatial_subset_type == LAT_LONG)
  {
      if (userDefinedProjParams == FALSE) {
          if (output_space_def->proj_num == 9)  { // Tranverse Mercator
             output_space_def->proj_param[2] = 1.0;
          }
          else {
             output_space_def->proj_param[2] = output_space_def->ul_corner.y;
          }          
          param->output_space_def.proj_param[3] = output_space_def->ul_corner.y;
          param->output_space_def.proj_param[4] = output_space_def->ul_corner.x;
          param->output_space_def.proj_param[5] = output_space_def->ul_corner.y;
          param->output_space_def.proj_param[6] = 0.0;
          param->output_space_def.proj_param[7] = 0.0;
          if (output_space_def->proj_num == 31) { // Integerized Sinusoidal
              param->output_space_def.proj_param[8] = 2.0;
              param->output_space_def.proj_param[10] = 0.0;
          }
      }

      /* Convert the output projection parameter lat/long values from decimal
         degrees to DMS */
      if (!Deg2DMS (param->output_space_def.proj_num,
                 param->output_space_def.proj_param)) {
        FreeParam(param);
        LOG_RETURN_ERROR("error converting projection parameters from decimal degrees to DMS",
                         "ConvertCorners", false);
      }
      
      /* Copy the projection parameters to orig_proj_param to use the decimal
         degree values later (GeoTiff output) */
      for (i = 0; i < NPROJ_PARAM; i++) {
          param->output_space_def.orig_proj_param[i] =
          param->output_space_def.proj_param[i];
      }      
      
     /* Get the UR and LL lat/longs from the UL and LR lat/longs.
        Convert the corner points to RADIANS for the GCTP call. */
     output_space_def->ul_corner.x *= RAD;
     output_space_def->ul_corner.y *= RAD;
     output_space_def->lr_corner.x *= RAD;
     output_space_def->lr_corner.y *= RAD;
     ul.x = output_space_def->ul_corner.x;
     ul.y = output_space_def->ul_corner.y;
     lr.x = output_space_def->lr_corner.x;
     lr.y = output_space_def->lr_corner.y;

     /* Use the first pixel size for these calculations */
     if (output_space_def->proj_num == PROJ_GEO)
     {
       /* Convert Geographic pixel size from degrees to radians */
       output_space_def->pixel_size = param->output_pixel_size[0] * RAD;

       /* If the output projection is Geographic, then the pixel size needs
          to be in degrees. Verify that the pixel size is less than 1.0. */
       if (param->output_pixel_size[0] > 1.0)
       {
         LOG_RETURN_ERROR("for output to geographic the pixel size needs to be "
		      "in degrees", "ConvertCorners", false);
       }
     }
     else
     {
       output_space_def->pixel_size = param->output_pixel_size[0];

       /* If the output projection is non-Geographic, then the pixel size
          needs to be in meters. Verify that the pixel size is larger
          than 1.0. */
       if (param->output_pixel_size[0] < 1.0)
       {
         LOG_RETURN_ERROR("for output to non-geographic projections the pixel "
                      "size needs to be in meters", "ConvertCorners", false);
       }
     }

     /* Initialize the number of lines and samples to 1,1 just so SetupSpace
        won't complain.  These lines and samples are not used in for_trans. */
     output_space_def->img_size.s = 1;
     output_space_def->img_size.l = 1;

     /* Get the forward and reverse transformation functions */
     out_space = SetupSpace(output_space_def);
     if (out_space == (Space_t *)NULL)
         LOG_RETURN_ERROR("setting up output space", "ConvertCorners", false);

     /* UL --------------------------------------------------*/
     #ifdef DEBUG
        printf ("Input UL lat/long (DEG): %f %f\n", ul.y*DEG, ul.x*DEG);
        printf ("Input UL lat/long (RAD): %f %f\n", ul.y, ul.x);
     #endif
     geo_coord_p.lon = ul.x;
     geo_coord_p.lat = ul.y;
     geo_coord_p.is_fill = false;
     if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
         &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
         FreeSpace(out_space);
         LOG_RETURN_ERROR("converting UL to output map coordinates",
                      "ConvertCorners", false);
     }
     #ifdef DEBUG
        printf ("Output UL projection coords x/y: %f %f\n", output_coord_p.x, output_coord_p.y);
     #endif
     output_coord_UL.x = output_coord_p.x;
     output_coord_UL.y = output_coord_p.y;
     if (output_coord_p.x < minx) minx = output_coord_p.x;
     if (output_coord_p.x > maxx) maxx = output_coord_p.x;
     if (output_coord_p.y < miny) miny = output_coord_p.y;
     if (output_coord_p.y > maxy) maxy = output_coord_p.y;

     /* LR -------------------------------------------------------*/
     #ifdef DEBUG
        printf ("Input LR lat/long (DEG): %f %f\n", lr.y*DEG, lr.x*DEG);
        printf ("Input LR lat/long (RAD): %f %f\n", lr.y, lr.x);
     #endif
     geo_coord_p.lon = lr.x;
     geo_coord_p.lat = lr.y;
     geo_coord_p.is_fill = false;
     if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
         &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
         FreeSpace(out_space);
         LOG_RETURN_ERROR("converting LR to output map coordinates",
                      "ConvertCorners", false);
     }
     #ifdef DEBUG
        printf ("Output LR projection coords x/y: %f %f\n", output_coord_p.x, output_coord_p.y);
     #endif
     output_coord_LR.x = output_coord_p.x;
     output_coord_LR.y = output_coord_p.y;
   
     if (userDefinedProjParams == TRUE) {
         if (output_coord_p.x < minx) minx = output_coord_p.x;
         if (output_coord_p.x > maxx) maxx = output_coord_p.x;
         if (output_coord_p.y < miny) miny = output_coord_p.y;
         if (output_coord_p.y > maxy) maxy = output_coord_p.y;        
     }
     else {
         output_coord_MP.x = output_coord_UL.x + (output_coord_LR.x - output_coord_UL.x)/2;
         output_coord_MP.y = output_coord_UL.y + (output_coord_LR.y - output_coord_UL.y)/2;

         if (out_space->inv_trans(output_coord_MP.x, output_coord_MP.y, &geo_coord_p.lon, &geo_coord_p.lat) != GCTP_OK) {
            FreeSpace(out_space);
            LOG_RETURN_ERROR("converting LR to output lat/long coordinates", "ConvertCorners", false);
         }   
         #ifdef DEBUG
            printf("Computed MidPoint lon,lat: %f, %f \n", geo_coord_p.lon*DEG, geo_coord_p.lat*DEG);
         #endif

         for (i = 0; i < NPROJ_PARAM; i++) {
             param->output_space_def.proj_param[i] = param->output_space_def.orig_proj_param[i];
         }
         param->output_space_def.proj_param[2] = geo_coord_p.lat*DEG;
         param->output_space_def.proj_param[3] = geo_coord_p.lat*DEG;         
         param->output_space_def.proj_param[4] = geo_coord_p.lon*DEG;
         param->output_space_def.proj_param[5] = geo_coord_p.lat*DEG;

          /* Copy the projection parameters to orig_proj_param to use the decimal
             degree values later (GeoTiff output) */
          for (i = 0; i < NPROJ_PARAM; i++) {
             param->output_space_def.orig_proj_param[i] = param->output_space_def.proj_param[i];
          }

          /* Convert the output projection parameter lat/long values from decimal
             degrees to DMS */
          if (!Deg2DMS (param->output_space_def.proj_num,
                     param->output_space_def.proj_param)) {
            FreeParam(param);
            LOG_RETURN_ERROR("error converting projection parameters from decimal degrees to DMS",
                             "ConvertCorners", false);
          }         

         /* reset projection center to MP. Recalculate map coord UL, LR*/
         out_space = SetupSpace(output_space_def);
         if (out_space == (Space_t *)NULL)
             LOG_RETURN_ERROR("setting up output space", "ConvertCorners", false);


         minx = FLT_MAX,
         maxx = -FLT_MAX,
         miny = FLT_MAX,
         maxy = -FLT_MAX;

         /* UL --------------------------------------------------*/
         #ifdef DEBUG
            printf ("Input UL lat/long (DEG): %f %f\n", ul.y*DEG, ul.x*DEG);
            printf ("Input UL lat/long (RAD): %f %f\n", ul.y, ul.x);
         #endif
         geo_coord_p.lon = ul.x;
         geo_coord_p.lat = ul.y;
         geo_coord_p.is_fill = false;
         if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
             &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
             FreeSpace(out_space);
             LOG_RETURN_ERROR("converting UL to output map coordinates",
                          "ConvertCorners", false);
         }
         #ifdef DEBUG
            printf ("Output UL projection coords x/y: %f %f\n", output_coord_p.x, output_coord_p.y);
         #endif
         output_coord_UL.x = output_coord_p.x;
         output_coord_UL.y = output_coord_p.y;
         if (output_coord_p.x < minx) minx = output_coord_p.x;
         if (output_coord_p.x > maxx) maxx = output_coord_p.x;
         if (output_coord_p.y < miny) miny = output_coord_p.y;
         if (output_coord_p.y > maxy) maxy = output_coord_p.y;     

         /* LR -------------------------------------------------------*/
         #ifdef DEBUG
            printf ("Input LR lat/long (DEG): %f %f\n", lr.y*DEG, lr.x*DEG);
            printf ("Input LR lat/long (RAD): %f %f\n", lr.y, lr.x);
         #endif
         geo_coord_p.lon = lr.x;
         geo_coord_p.lat = lr.y;
         geo_coord_p.is_fill = false;
         if (out_space->for_trans(geo_coord_p.lon, geo_coord_p.lat,
             &output_coord_p.x, &output_coord_p.y) != GCTP_OK) {
             FreeSpace(out_space);
             LOG_RETURN_ERROR("converting LR to output map coordinates",
                          "ConvertCorners", false);
         }
         #ifdef DEBUG
            printf ("Output LR projection coords x/y: %f %f\n", output_coord_p.x, output_coord_p.y);
         #endif
         output_coord_LR.x = output_coord_p.x;
         output_coord_LR.y = output_coord_p.y;     

         if (output_coord_p.x < minx) minx = output_coord_p.x;
         if (output_coord_p.x > maxx) maxx = output_coord_p.x;
         if (output_coord_p.y < miny) miny = output_coord_p.y;
         if (output_coord_p.y > maxy) maxy = output_coord_p.y;

     } /* if LAT_LONG */

     
     /* The pixel size and number of lines/samples is specified for each
        SDS.  The overall UL and LR corners will be the same for each SDS,
        however the pixel size and number of lines/samples might be
        different. */
     for (i = 0; i < param->num_input_sds; i++)
     {
       if (output_space_def->proj_num == PROJ_GEO)
       {
         /* Convert Geographic pixel size from degrees to radians */
         param->output_pixel_size[i] *= RAD;
       }

       /* Calculate the number of output lines and samples */
       param->output_img_size[i].s =
         (int) ((maxx - minx) / param->output_pixel_size[i] + 0.5);
       param->output_img_size[i].l =
         (int) ((maxy - miny) / param->output_pixel_size[i] + 0.5);

       #ifdef DEBUG
           if (output_space_def->proj_num == PROJ_GEO)
             printf ("Pixel size is %f\n", param->output_pixel_size[i]*DEG);
           else
             printf ("Pixel size is %f\n", param->output_pixel_size[i]);
           printf ("Output number of lines is %d\n",
             param->output_img_size[i].l);
           printf ("Output number of samples is %d\n",
             param->output_img_size[i].s);
       #endif
     }

     /* Use the first pixel size and number of lines/samples for the
        rest of the calculations */
     output_space_def->pixel_size = param->output_pixel_size[0];
     output_space_def->img_size.l = param->output_img_size[0].l;
     output_space_def->img_size.s = param->output_img_size[0].s;

     /* Redefine the output coords to use the minimum bounding box in
        output projection coords. Also make the LR corner a factor of the
        number of lines and samples. */
     output_space_def->ul_corner.x = output_coord_UL.x;
     output_space_def->ul_corner.y = output_coord_UL.y;  
     
     output_space_def->lr_corner.x = output_space_def->ul_corner.x +
         output_space_def->img_size.s * output_space_def->pixel_size;
     
     output_space_def->lr_corner.y = output_space_def->ul_corner.y -
         output_space_def->img_size.l * output_space_def->pixel_size;

     #ifdef DEBUG
         if (output_space_def->proj_num == PROJ_GEO)
         {
           printf ("Output UL projection coords: %f %f\n",
             output_space_def->ul_corner.x*DEG, output_space_def->ul_corner.y*DEG);
           printf ("Output LR projection coords: %f %f\n",
             output_space_def->lr_corner.x*DEG, output_space_def->lr_corner.y*DEG);
         }
         else
         {
           printf ("Output UL projection coords: %f %f\n",
             output_space_def->ul_corner.x, output_space_def->ul_corner.y);
           printf ("Output LR projection coords: %f %f\n",
             output_space_def->lr_corner.x, output_space_def->lr_corner.y);
         }
     #endif
  }

  /* If the output spatial subset type is PROJ_COORDS, use the UL and LR
     projection coords to determine the number of lines and samples in the
     output image. Also determine the lat/long values. */
  else if (param->output_spatial_subset_type == PROJ_COORDS)
  {
     /* If the input was in projection coords and the projection is geographic,
        then the corners are in degrees.  We need to make sure they are in
        radians for the rest of the processing. */
     if (output_space_def->proj_num == PROJ_GEO)
     {
        output_space_def->ul_corner.x *= RAD;
        output_space_def->ul_corner.y *= RAD;
        output_space_def->lr_corner.x *= RAD;
        output_space_def->lr_corner.y *= RAD;
     }

     /* The pixel size and number of lines/samples is specified for each
        SDS.  The overall UL and LR corners will be the same for each SDS,
        however the pixel size and number of lines/samples might be
        different. */
     for (i = 0; i < param->num_input_sds; i++)
     {
       /* If the output projection is Geographic, then the pixel size needs
          to be in degrees. Verify that the pixel size is larger than 1.0. */
       if (output_space_def->proj_num == PROJ_GEO &&
           param->output_pixel_size[i] > 1.0)
       {
         FreeSpace(out_space);
         LOG_RETURN_ERROR("for output to geographic the pixel size needs to be "
                      "in degrees", "ConvertCorners", false);
       }
       else if (output_space_def->proj_num == PROJ_GEO)
       {
         /* Convert Geographic pixel size from degrees to radians */
         param->output_pixel_size[i] *= RAD;
       }

       /* Determine the number of lines and samples */
       param->output_img_size[i].l =
          (int) (((fabs (output_space_def->lr_corner.y -
          output_space_def->ul_corner.y)) / param->output_pixel_size[i]) +
          0.5);
       param->output_img_size[i].s =
          (int) (((fabs (output_space_def->lr_corner.x -
          output_space_def->ul_corner.x)) / param->output_pixel_size[i]) +
          0.5);
     }

     /* Use the first pixel size and number of lines/samples for the
        rest of the calculations */
     output_space_def->pixel_size = param->output_pixel_size[0];
     output_space_def->img_size.l = param->output_img_size[0].l;
     output_space_def->img_size.s = param->output_img_size[0].s;

     /* Recalculate the LR corner based on the number of lines and samples,
        since the projection coords might not have been an exact number of
        lines and samples. */
     output_space_def->lr_corner.x = output_space_def->ul_corner.x +
         output_space_def->img_size.s * output_space_def->pixel_size;
     output_space_def->lr_corner.y = output_space_def->ul_corner.y -
         output_space_def->img_size.l * output_space_def->pixel_size;

     /* Get the forward and reverse transformation functions */
     out_space = SetupSpace(output_space_def);
     if (out_space == (Space_t *)NULL)
         LOG_RETURN_ERROR("setting up output space", "ConvertCorners", false);
  }

  /* Determine/Redetermine the lat/long of the UL corner */
  ul.x = output_space_def->ul_corner.x;
  ul.y = output_space_def->ul_corner.y;
  if (out_space->inv_trans(ul.x, ul.y, &geo_coord_p.lon, &geo_coord_p.lat)
      != GCTP_OK) {
      FreeSpace(out_space);
      LOG_RETURN_ERROR("converting UL to output lat/long coordinates",
                   "ConvertCorners", false);
  }
  output_space_def->ul_corner_geo.lat = geo_coord_p.lat;
  output_space_def->ul_corner_geo.lon = geo_coord_p.lon;

  /* Determine the lat/long of the LR corner */
  lr.x = output_space_def->lr_corner.x;
  lr.y = output_space_def->lr_corner.y;
  if (out_space->inv_trans(lr.x, lr.y, &geo_coord_p.lon, &geo_coord_p.lat)
      != GCTP_OK) {
      FreeSpace(out_space);
      LOG_RETURN_ERROR("converting LR to output lat/long coordinates",
                   "ConvertCorners", false);
  }
  output_space_def->lr_corner_geo.lat = geo_coord_p.lat;
  if (output_space_def->proj_num == PROJ_GEO && output_space_def->straddlesDateline) {
      if (geo_coord_p.lon < 0) {
          output_space_def->lr_corner_geo.lon = 2*3.14159267 + geo_coord_p.lon;
          output_space_def->lr_corner.x += 2*3.14159267;
      }
  }
  
  /* Free the output space */
  FreeSpace(out_space);

  return true;
}
