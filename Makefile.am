
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS = unrotated_square \
                  rotated_square \
                  rotated_square_upgraded \
                  circular_disc \
                  star_trek \
				  pull_strip \
				  free_swelling_square

#-KS sources---------------------------------------------------------------
unrotated_square_SOURCES = \
 unrotated_square.cc
rotated_square_SOURCES = \
 rotated_square.cc
rotated_square_upgraded_SOURCES = \
 rotated_square_upgraded.cc
circular_disc_SOURCES = \
 circular_disc.cc
star_trek_SOURCES = \
 star_trek.cc
pull_strip_SOURCES = \
 pull_strip.cc
free_swelling_square_SOURCES = \
 free_swelling_square.cc

#-KS libraries-------------------------------------------------------------
unrotated_square_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
rotated_square_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
rotated_square_upgraded_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
circular_disc_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
star_trek_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
pull_strip_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
free_swelling_square_LDADD = -L@libdir@ -lc1_koiter_steigmann \
 -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#circular_sector_LDADD = -L@libdir@ -lc1_koiter_steigmann \
# -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#corner_modes_LDADD = -L@libdir@ -lc1_koiter_steigmann \
# -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#square_corner_modes_LDADD = -L@libdir@ -lc1_koiter_steigmann \
# -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#triangular_plate_LDADD = -L@libdir@ -lc1_koiter_steigmann \
# -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#element_test_LDADD = -L@libdir@ -lc1_koiter_steigmann \
# -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
#jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
# jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.

AM_CPPFLAGS += -I@includedir@ -Wall -ggdb3 -O0
#EXTRA_DIST = Circle1.1.ele Circle1.1.node Circle1.1.poly
