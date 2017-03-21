#-----------------------------------------------------------------------
# Makefile for GCTP library.
#
# Note: Makefile.lnx can be used to compile on the MacOS X architecture.
#-----------------------------------------------------------------------

CC      = gcc
CFLAGS  = -O3 -Wall -W -Wno-switch
CP      = cp
RM      = rm -f
MV      = mv

#-----------------------------------------------------------------------------
# Rules for compiling the object files.  If you have multiple groups of object
# files, add one dependency of the form "$(OBJS): $(INC)" for each group of
# object files.
#-----------------------------------------------------------------------------

all: spload nad27sp.lut nad83sp.lut

# only build if the source file is present to allow for the check-in process
spload:
	if [ -f spload.c ]; then \
            $(CC) $(CFLAGS) spload.c -o spload; \
        fi

# only run if the spload is present to allow for the check-in process
nad27sp.lut: spload
	if [ -f spload ]; then \
            ./spload; fi

# only run if the spload is present to allow for the check-in process
nad83sp.lut: spload
	@if [ -f spload ]; then \
            ./spload; fi

# need to copy the nad27 and nad83 binary tables to the MRT_DATA_DIR
install:
	$(CP) nad27sp.lut ../data/nad27sp
	$(CP) nad83sp.lut ../data/nad83sp

clean:
	$(RM) spload nad27sp.lut nad83sp.lut 

copy-makefile:
	@if [ ! -f Makefile.orig ]; then \
	  $(CP) Makefile Makefile.orig; fi
	$(CP) Makefile.$(CPMAKEFILEEXT) Makefile


