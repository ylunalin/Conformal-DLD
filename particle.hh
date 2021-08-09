#ifndef PARTICLE_HH
#define PARTICLE_HH
struct particle {
	/** The particle's ID. */
//	int id;
	/** The x position of the particle. */
	double x;
	/** The y position of the particle. */
	double y;
	particle() {}
	particle(double x_,double y_) :
		 x(x_), y(y_) {}
};
#endif
