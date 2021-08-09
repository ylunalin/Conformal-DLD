#include "common.hh"

/** \brief Opens a file and checks the operation was successful.
 *
 * Opens a file, and checks the return value to ensure that the operation
 * was successful.
 * \param[in] filename the file to open.
 * \param[in] mode the cstdio fopen mode to use.
 * \return The file handle. */
FILE* safe_fopen(const char* filename,const char* mode) {
	FILE *temp=fopen(filename,mode);
	if(temp==NULL) fprintf(stderr,"rmap_sim: error opening file \"%s\"",filename);
	return temp;
}

/** \brief Function for printing fatal error messages and exiting.
 *
 * Function for printing fatal error messages and exiting.
 * \param[in] p a pointer to the message to print.
 * \param[in] status the status code to return with. */
void fatal_error(const char *p,int code) {
	fprintf(stderr,"dld_sim: %s\n",p);
	exit(code);
}

/** \brief Calculates the argument of two-dimensional position vector.
 *
 * Calculates the argument of the two-dimensional position vector.
 * \param[in] (x,y) the coordinates of the vector.
 * \return The argument. */
double argument(double x,double y) {
	const double pi=3.1415926535897932384626433832795;
	return x+y>0?(x>y?atan(y/x):0.5*pi-atan(x/y)):(x>y?-atan(x/y)-0.5*pi:atan(y/x)+(y>0?pi:-pi));
}
