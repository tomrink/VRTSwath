
CC = gcc
CFLAGS = -O3 -Wall -W -Wno-switch -g
#Debug: CFLAGS = -O0 -Wall -W -Wno-switch -g
LDFLAGS = $(HDF5LIB) $(NETCDF4LIB) $(GEOLIB) $(TIFFLIB) -g -lsz -lz -ldl -lm -lcurl -Wl,-rpath,/home/rink/hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared/lib
MV = mv
CP = cp
RM = rm -f
INCS  = -I../include -I/opt/netcdf4/4.4.1.1-gcc-4.9.2/include -I../../hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared/include


#--------------------------
# Define the include files:
#--------------------------
INC	= resamp.h param.h geoloc.h input.h scan.h output.h space.h kernel.h \
	 patches.h myhdf.h mystring.h myerror.h bool.h range.h parser.h \
	 myproj.h myproj_const.h usage.h const.h deg2dms.h GeoS2G.h \
	 geowrpr.h myendian.h

#-----------------------------------------------------------
# Define the source code and object files for the resampler:
#-----------------------------------------------------------
SRC	= param.c geoloc.c input.c scan.c output.c space.c kernel.c \
	  patches.c myhdf.c mystring.c parser.c myerror.c \
	  InitGeoTiff.c deg2dms.c degdms.c convert_corners.c \
	  geo_trans.c write_hdr.c write_rb.c logh.c geowrpr.c \
	  filegeo.c myendian.c

OBJ = $(SRC:.c=.o)

#----------------------------------
# Define the resampler executables:
#----------------------------------
EXE = swath2grid

#-----------------------------
# Define the object libraries:
#-----------------------------
GEOLIB = ../lib/libgctp.a
HDF5LIB = -L../../hdf5-1.8.18-linux-centos7-x86_64-gcc485-shared/lib -lhdf5_hl -lhdf5
TIFFLIB = -L/usr/lib64 -ljpeg /usr/lib64/libtiff.so.5.2.0 ../lib/libgeotiff.a
NETCDF4LIB = -L/opt/netcdf4/4.4.1.1-gcc-4.9.2/lib -lnetcdf

#-------------------------------------------------------------------
# Targets for each executable.
# Note: stripping executable with -s flag on SGI may cause problems.
#-------------------------------------------------------------------
all: $(EXE)

swath2grid: resamp.o $(OBJ) $(INC)
	$(CC) -o swath2grid resamp.o $(OBJ) $(LDFLAGS)

install:
	$(MV) swath2grid dumpmeta hdf2hdr ../bin

clean:
	$(RM) *.o swath2grid dumpmeta hdf2hdr *~

copy-makefile:
	@if [ ! -f Makefile.orig ]; then \
	  $(CP) Makefile Makefile.orig; fi
	$(CP) Makefile.$(CPMAKEFILEEXT) Makefile

#-----------------------------------------------------------------------------
# Rules for compiling the object files.
#-----------------------------------------------------------------------------

"$(OBJ)": $(INC)

.c.o:
	$(CC) $(CFLAGS) $(INCS) -c $<


