# MakeDiskGalaxy -- initial conditions generator for isolated disk galaxies
#
#   make            build ./MakeDiskGalaxy
#   make clean      remove objects and binary
#   make CC=icc     build with a different compiler

EXEC     = MakeDiskGalaxy
SRCDIR   = src
BUILDDIR = build

CC       ?= gcc
OPTIMIZE ?= -O2 -g

# Compile-time options (uncomment / add as needed):
#   -DREDSHIFT_SCALING        scale halo properties with redshift
#   -DV_SCALING               set halo via V200 instead of M200
#   -DMY_STRUCTURE            use structure_my.c
#   -DOLD_STRUCTURE           use structure.c
#   -DADIABATIC_CONTRACTION   use structure_withcontraction.c
#   -DPURE_HALO               dark matter halo only
OPT      = -DREDSHIFT_SCALING -DV_SCALING -DMY_STRUCTURE

CFLAGS   = $(OPTIMIZE) $(OPT) -I$(SRCDIR)
LIBS     = -lm

SRCS = main.c cooling.c read_param.c effmodel.c globvars.c \
       halo.c haloset.c escape.c \
       disk.c diskset.c \
       gashalofns.c gashaloset.c gashalo.c \
       bulge.c bulgeset.c \
       gas.c gasset.c \
       save.c init.c structure.c structure_withcontraction.c structure_my.c \
       misc.c toomre.c force.c forcetree.c \
       nrsrc/bessi0.c nrsrc/bessi1.c nrsrc/bessj0.c nrsrc/bessj1.c \
       nrsrc/bessk0.c nrsrc/bessk1.c nrsrc/erff.c nrsrc/gammln.c \
       nrsrc/gammp.c nrsrc/gasdev.c nrsrc/gcf.c nrsrc/gser.c \
       nrsrc/indexx.c nrsrc/nrutil.c nrsrc/polint.c nrsrc/qromb.c \
       nrsrc/ran1.c nrsrc/spline.c nrsrc/splint.c nrsrc/trapzd.c \
       nrsrc/zriddr.c

OBJS = $(addprefix $(BUILDDIR)/,$(SRCS:.c=.o))
INCL = $(wildcard $(SRCDIR)/*.h $(SRCDIR)/nrsrc/*.h)

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(INCL) Makefile
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) $(EXEC)
