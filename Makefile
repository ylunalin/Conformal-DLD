include config.mk
lflags=-L. `gsl-config --libs`
iflags=`gsl-config --cflags`

# Lists of files to be built
objs=common.o kdtree.o worklist.o dld_sim.o dld_msgs.o dld_data_procs.o sim_params.o
# Replace objs suffix with .cc
src=$(patsubst %.o,%.cc,$(objs))

execs=radial flat
# recursive make command, stored in $(MAKE) variable
all:
#	$(MAKE) -C ../../shared
	$(MAKE) executables

executables: $(execs)

include Makefile.dep

depend:
	$(cxx) $(cflags) $(iflags) -MM $(src) >Makefile.dep

clean:
	rm -rf $(objs) $(execs) *.dSYM

# -c flag compiles but not link
%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

# first remove existing libbac.a, then create a new one from the all the object files
libdld.a: $(objs)
	rm -f libdld.a
	ar rs libdld.a $^
# -o flag compile to the file named specified by the string follws it

# while -L looks in directory for library files
# -l(library name) links a library file, e.g. -ldld links static library file libdld.a

# $@ is the automatic variable for target name
# $< is the name of the first prerequisite.

radial: radial.cc libdld.a sim_manager.hh
	$(cxx) $(cflags) $(iflags) -o $@ $< -ldld $(lflags)

flat: flat.cc libdld.a sim_manager.hh
	$(cxx) $(cflags) $(iflags) -o $@ $< -ldld $(lflags)

# Since c_mapping only needs common.o, not the entire libdld.a
# we set prerequiistes to be the entire list
c_mapping: c_mapping.cc common.o
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend executables libdld.a
