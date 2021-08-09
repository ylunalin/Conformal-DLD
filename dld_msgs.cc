#include "dld_sim.hh"

// A file that contains all the common error messages.
// So that the core code is cleaner

const char * dld_sim::DIV_LINE = "##################################################";
const char * dld_sim::LOC_NAN = "Random walker location is NAN.";
const char * dld_sim::TREE_NO_PBC = "PBC isn't implemented in quadtree search yet.";
const char * dld_sim::MINDIST_NO_PBC = "PBC isn't implemented in min_distance(...) search yet.";
const char * dld_sim::NO_PBC = "PBC isn't implemented in this part yet.";
const char * dld_sim::EQUI_DIST = "1st and 2nd nearest neighbor are about equi-distance. Use uniform step to an empty circle instead of pacman.";
const char * dld_sim::FAR= "Distance too far, use simple uniform step to an empty circle.";
const char * dld_sim::TOO_FAR= "Distance WAY too far, remap using first passage + conformal mapping.";
const char * dld_sim::NNB_SEARCH_FAIL = "Sweep didn't find 2nd nearest or either neighbors within maximum sweep radius.\nTaking a walk to the edge of an empty circle.";
const char * dld_sim::XSTEP_OUT_OF_BOUND = "Stepx out of bound. Decouple x,y. The row we are searching exists.";
const char * dld_sim::SPECIAL_SEARCH = "Random walker on the right most column, which should contains also the first particles in the left most column. Searching there as well.";
const char * dld_sim::TOO_FEW = "Total number of particles < 10. Almost no particles left. Should this happen?";
const char * dld_sim::BIN_MAT_ERR = "Binary Matrix is not written for hexagonal lattice.";
const char * dld_sim::INSIDE_ERR = "Random walker walked inside a particle. Check f.*_fail file.";
const char * dld_sim::INIT_FILE_ERR = "Warning, initialization file isn't opened, starting from seed at origin.";
const char * dld_sim::HEX_ERR = "Initialization of cluster failed. Periodic boundary condition is only implemented for flat geometry on square lattice.";
const char * dld_sim::SURVEY_ERR = "Height survey returns an array larger than the one provided to record data. Provide a larger container!";
const char * dld_sim::PARAM_NOT_SET = "In sim_params object, some geometry parameters are not set.";
