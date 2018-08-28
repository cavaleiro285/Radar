#############################################################################

.SUFFIXES:
.SUFFIXES: .c .o $(SUFFIXES)
.SUFFIXES: .f .o $(SUFFIXES)
.SUFFIXES: .c .a $(SUFFIXES)
.SUFFIXES: .c .a $(SUFFIXES)
.SUFFIXES: .f .a $(SUFFIXES)
.SUFFIXES: .so .sl $(SUFFIXES)

# Options 
F77_LIBDIR =/usr/local/lang/SC0.0   #SUNOS ONLY
CDEBUG =
F_OPTIONS = -g
# Platform dependent C_OPTIONS 
#ifdef __ALPHA
C_OPTIONS = -D__ALPHA
#endif	
#C_OPTIONS = -DSUNOS 
#C_OPTIONS = -DAIX 
#C_OPTIONS = -DHPUX

# Dec Alpha/OSF Flags
#ifdef __ALPHA
 CFLAGS_alpha = -c
 FFLAGS_alpha = -c
 LDFLAGS_alpha = -shared
 F_LIBS_alpha = -lUfor -lfor -lFutil -lm -lots -lc 
#endif

# HP/UX Flags
##ifdef HPUX
# CFLAGS_hpux = +z -c
# FFLAGS_hpux = +z -c
# LDFLAGS_hpux = -b
##endif

# SGI Flags
##ifdef SGI
# CFLAGS_sgi = -c -KPIC
# FFLAGS_sgi = -c -KPIC
# LDFLAGS_sgi = -shared
# F_LIBS = -lftn
##endif

# SunOs Flags
##ifdef SUNOS
# CFLAGS_sunos = -c -pic -fsingle
# FFLAGS_sunos = -c -pic
# F_LIB_sunos = $(F77_LIBDIR)/libF77.a $(F77_LIBDIR)/libm.a 
##endif

# AIX Flags
##ifdef AIX
# F77_aix = xlf
# FFLAGS_aix = -c
# CFLAGS_aix = -c
# LDFLAGS_aix = -bM:SRE 
# F_LIBS_aix = -e _nostart -lxlf
##endif

# Source Code
FSOURCE = cg_alt_dip.f mlt.f solar_loc.f mlt1.f sfc_convert_geo_coord.f \
          rylm.f eqn_of_time.f cnvcoord.f cnv_sec_mdhms.f L_mlat_mlt.f \
          coodecc.f 
CSOURCE = cnv_mdhms_sec_idl.c cnv_mdhms_sec.c cnvcoord_idl.c mlt_idl.c \
           L_mlat_mlt_idl.c coodecc_idl.c

F_OBJS = cg_alt_dip.o mlt.o solar_loc.o mlt1.o sfc_convert_geo_coord.o \
         rylm.o eqn_of_time.o cnvcoord.o cnv_sec_mdhms.o L_mlat_mlt.o \
         coodecc.o


# Alpha specific
#F_OBJS = fsource.o

C_OBJS = cnv_mdhms_sec.o cnv_mdhms_sec_idl.o cnvcoord_idl.o mlt_idl.o \
         L_mlat_mlt_idl.o coodecc_idl.o

HOMER = /usr/users/johnson/idl/darn/Radar
TARGET = $(HOMER)/LIB_PGM.so
#TARGET = LIB_PGM.so

# Alpha specific make
#ifdef __ALPHA
$(TARGET): $(C_OBJS) $(F_OBJS)
	$(LD) -o $@ $(LDFLAGS_alpha) $(C_OBJS) $(F_OBJS) $(F_LIBS_alpha)
	$(RM) *.o

$(F_OBJS): $(FSOURCE)
	$(FC) $(FFLAGS_alpha) $(F_OPTIONS) $(FSOURCE)

# Alpha specfic
#$(F_OBJS): $(FSOURCE)
#        $(FC) $(FFLAGS_alpha) -o $(F_OBJS) -O4 $(FSOURCE)

$(C_OBJS): $(CSOURCE) 
	$(CC) $(CFLAGS_alpha) $(C_OPTIONS) $(CSOURCE)
#endif

# SUNOS specific make
##ifdef SUNOS
#$(TARGET): $(C_OBJS) $(F_OBJS)
#	$(LD) -o $@ $(C_OBJS) $(F_OBJS) $(F_LIB_sunos)
#$(F_OBJS): $(FSOURCE)
#	$(FC) $(FFLAGS_sunos) $(F_OPTIONS) $(FSOURCE)
#$(C_OBJS): $(CSOURCE)
#	$(CC) $(CFLAGS_sunos) $(C_OPTIONS) $(CSOURCE)
##endif

# AIX specifc make
##ifdef AIX
#$(TARGET): $(C_OBJS) $(F_OBJS)
#	$(LD) -o $@ $(LDFLAGS_aix) $(C_OBJS) $(F_OBJS) -BE:example.exp \
#        $(F_LIBS_aix) 
#$(F_OBJS): $(FSOURCE)
#	$(F77_aix) $(FFLAGS_aix) $(F_OPTIONS) $(FSOURCE)
#$(C_OBJS): $(CSOURCE)
#	$(CC) $(CFLAGS_aix) $(C_OPTIONS) $(CSOURCE)
##endif
