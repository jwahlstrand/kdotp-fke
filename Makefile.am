INCLUDES = \
	-I$(top_srcdir) \
	-Wall \
        -O3 \
        -ffast-math \
    $(GLIB_CFLAGS) \
    $(GSL_CFLAGS) \
    $(FFTW_CFLAGS) 

c14bandlibdir = $(pkglibdir)
noinst_LTLIBRARIES = libkdotp.la

savebands_SOURCES = savebands.c
savebands_LDADD = libkdotp.la $(GLIB_LIBS) $(GSL_LIBS) $(FFTW_LIBS)
savebands_LDFLAGS = -static  -lm -lc -lgcc -lgcc_eh

saveV_SOURCES = saveV.c
saveV_LDADD = libkdotp.la $(GLIB_LIBS) $(GSL_LIBS) $(FFTW_LIBS)
saveV_LDFLAGS = -static  -lm -lc -lgcc -lgcc_eh

bin_PROGRAMS = savebands saveV

libkdotp_la_SOURCES = \
    three-vector.c \
    k-dot-p.c \
    zincblende2pba.c \
    zincblende14.c \
    zincblende14nr.c \
    zincblende14ns.c \
    diagonalize.c \
    calc-gamma.c \
    gamma-spline.c \
    matrix-element-pert.c \
    zener.c \
    pert-ode.c \
    coeffs-pert.c \
    coeffs-pert-other.c
    
EXTRA_DIST = px8.txt \
             py8.txt \
             pz8.txt \
             px14.txt \
             py14.txt \
             pz14.txt \
             k-dot-p.h 

MAINTAINERCLEANFILES = Makefile.in
