/** A File that contains all the data collection, processing utilities.
 These tools include using square voronoi cells to represent mass density
 of each particle, and create density bump functions according to lattice
 tilt; surveying the whole simulation to collect CDF data; etc. */
#include "dld_sim.hh"

/** Gives a bump function by scanning across a particle (considered a tilted square) vertically.
 * Used in flat interface width computation.
 * \param[in,out] cdf the pointer to the starting position of the CDF to be filled by this bump.
 * \param[in] y the y coordinate of current particle.
 * \param[in] h the starting y coordinate.
 * \param[in] dh the height increment.
 * \param[in] avail_n integer number of available length of CDF array from starting position. */
void dld_sim::bump(double * cdf, double y, double h, double dh, int avail_n){
	double eps = epsilon/height;
	if (true_tilt!=0) {
		double costheta = cos(true_tilt);
		double sintheta = sin(true_tilt);
		// the top level of the bump
		double level = 1.0 / costheta;
		// the slope of bump function
		double slope = level/sintheta;
		// distance between center and lower vertex
		double h_dis = sqrt(2)*cos(pi/4 - true_tilt)/2;
		// the total distance (vertical) between top and bottom vertices
		double h_dis2 = h_dis*2;
		// the height of the bottom vertex
		double h_low = y - h_dis;

		double dis = h-h_low;
		int count = 0;
		while(dis < h_dis2 && count<avail_n){
			double val=0;
			// dist from h to h_low
			if (dis>=0 && dis<sintheta){
				val=slope*dis;
			}
			else if(dis>=sintheta && dis<=costheta){
				val = level;
			}
			else if(dis>costheta && dis<=h_dis2){
				val=(h_dis2 - dis)*slope;
			}
			else{
				val=eps;
			}
			*(cdf+count) = *(cdf+count) + val;
			dis+=dh;
			count++;
		}
	}
	else{
		// the total distance (vertical) between top and bottom vertices
		double h_dis = 0.5;
		double h_dis2 = 1;
		// the height of the bottom vertex
		double h_low = y - h_dis;

		double dis = h-h_low;
		int count = 0;
		while(dis < h_dis2 && count<avail_n){
			double val=0;
			// dist from h to h_low
			if(dis>=0 && dis<=h_dis2) {val=1.0;}
			else{
				val=eps;
			}
			*(cdf+count) = *(cdf+count) + val;
			dis+=dh;
			count++;
		}


	}
}

/** Starting at the lowest y on the interface, surveys the density of particles as function of height,
 * which is equivalent to one instance of a discrete CDF of particles at the interface.
 * \param[in] dh height increment.
 * \param[in,out] n the length of CDF array (returned as pointer), set by ylow. */
double * dld_sim::survey(double dh, int &n, bool init){
	double *cdf;
	double *hi_arr;

	// This is the length from the particle true position, to the lowest vertex if the mass were to be smeared out to a square of side = 1.
	const double Y_DIS = 0.5*sqrt(2)*cos(pi/4- true_tilt);
	double local_ylow = ylow;
	int local_blclow = blclow;
	if (init) {
		local_ylow = -rad;
		local_blclow = 0;
	}
	const double BOTTOM_LINE = double (std::floor(local_ylow - Y_DIS));
	if (BOTTOM_LINE < double (std::floor (-rad-Y_DIS))) fatal_error("survey()", "Height survey starts lower than the bottom of the simulation");

	// Add the lowest vertex y coordinate of the lowest dissolved particle as the bottom of the interface we care about
	const int NELEM = int ((mr-BOTTOM_LINE)/dh);
	if (NELEM==0) return NULL;
	else if (NELEM<0) fatal_error("survey()", "NELEM < 0, check mr and ylow.");
	cdf = new double [NELEM];
	hi_arr = new double [NELEM];

	// Add a point in height every dh
	for (int i =0; i<NELEM; i++){
		hi_arr[i] = BOTTOM_LINE+i*dh;
		cdf[i] = 0.0;
	}

	// It's more natural to scan over mn boxes one by one
	for (int i =0; i<mn; i++){

		// If it's way below the lowest box, it can be ignored
		// only care about particles above the lowest particles, since all particles would add to unity below
		// we still need to get normalization factor
		if (int (i/nh) >= local_blclow-1) {
			for (int j=0; j<co[i];j++){

				// if the particle's top vertex is below the lowest dissolved particle, it can also be ignored
				double y = ba[i][j].y;
				// if the bottom line go through the particle
				// we feed the pointer cdf into bump function
				if (y + Y_DIS  > BOTTOM_LINE  &&  y- Y_DIS < BOTTOM_LINE){
					bump(cdf, y, BOTTOM_LINE, dh, NELEM);
				}
				// if the lower vertex of the particle square is above bottom line
				// we need to figure out where the bump should be added to the cdf array
				else if (y-Y_DIS > BOTTOM_LINE) {
					int ind = int ( (y-Y_DIS - BOTTOM_LINE)/dh);
					bump(cdf+ind, y, hi_arr[ind], dh, NELEM-ind);
				}
			}
		}
	}
	delete [] hi_arr;
	double * re_cdf=0;
	// nonzero angle has little risk of having a h going through the boundary between rows of particles
	int nonone=0;
	int nonzero=0;
	double eps = epsilon;
	double cdffull = cdf[0];
	if (true_tilt==0) eps=1.0/twidth;


	if (!init) {
		for (int i =0; i<NELEM; i++){
			cdf[i] = cdf[i]/cdffull;
			// find the position before cdf becomes smaller than 1
			if(cdf[i]<=1+eps&&cdf[i]>=1-eps){
				nonone=i;
			}
			// find the position when cdf first becomes 0
			if(cdf[i]==0) {
				nonzero = i+1;
				break;
			}
		}
	}

	else{
		// If it's the initial survey, cdf[0] would not always be nonzero, but then we know it should roughly be twidth
		cdffull = twidth;
		for (int i =0; i<NELEM; i++) {
			cdf[i] = cdf[i]/cdffull;
			if(cdf[i]>=1-eps && cdf[i]<=1+eps) nonone = i;
			// Since there could be some nonzero entries at the beginning of CDF, we requrie nonzero>nonone
			if(cdf[i]==0 && i>nonone) {
				nonzero = i+1;
				break;
			}
		}
	}

	n = nonzero-nonone;
	if(n>0) {
		re_cdf = new double[n];
		for (int i=0; i<n; i++) {
			*(re_cdf+i) = *(cdf+nonone+i);
		}

	}
	delete [] cdf;

	return re_cdf;
}

/** Write pdf at kth frame to a 2D array f.
* \param[in] f the 2D array to be written.
* \param[in] k the integer representing frame number. */
void dld_sim::set_pdf_field(double **f, int k){
	// Grid dimensions for PDF
	const double ax=rad*(-1),ay=-1*rad;

	// Horizontal and vertical grid size
	int mm=(int) 2*rad+1;

	//Increment 1/nr per particle to the 2D grid f.
	for(int s=0;s<mn;s++) {
		for(int q=0;q<co[s];q++) {
			particle &b=ba[s][q];
			int grid_id = (int) mm*(b.y-ay)+(b.x-ax)+1;
			f[k][grid_id] += tot_nrsp;
		}
	}
}

/** Compute radius of gyration of a cluster of particles. */
double dld_sim::rgyration(){
	double rgsq=0;
	double xmean=0, ymean=0;
	double normal =(double) total_particles();
	normal = 1/normal;
	for(int i=0;i<mn;i++){
		for(int j=0;j<co[i];j++){
			xmean+=ba[i][j].x;
			ymean+=ba[i][j].y;
		}
	}
	xmean*=normal;
	ymean*=normal;
	for(int i=0;i<mn;i++){
		for(int j=0;j<co[i];j++){
			rgsq += ((ba[i][j].x-xmean)*(ba[i][j].x-xmean)+(ba[i][j].y-ymean)*(ba[i][j].y-ymean));
		}
	}
	rgsq*=normal;
	double rg = sqrt(rgsq);
	return rg;
}

/** Record data for P(r, N)dr.
 * We need to keep track of the distance from Nth particle to the center of mass (radial cluster)
 * or distance to the average interface height, defined wrt an arbitrary horizon.
 * Store this information in some array managed by the sim manager class.
 * \param[in,out] nN the index of the element to write in the active zone array, gets incremented by 1.
 * \param[in] p the pointer to the particle that is being dissolved or aggregated next. */
void dld_sim::active_zone_dist(int &nN, particle *p){
	// quit if we exceeding the length of data collection
	if (nN >= spars_p -> N_len) return;
	double * rp = spars_p -> r_nr;
	double *r0p = spars_p-> r0_nr;
	double * xp = spars_p -> x_nr;
	double * yp = spars_p -> y_nr;
	int N = spars_p -> N[nN];

	// 1. store (x,y) of the particle
	double myx = p->x, myy = p->y;
	xp[nN] = myx; yp[nN] = myy;

	if(n_dissd == N-1){
		if(!is_flat){
			// 2. compute center of mass or average height, at N-1 particle dissolved.
			double xmean=0, ymean=0;
			double normal =1.0/ (double) total_particles();
			for(int i=0;i<mn;i++){
				for(int j=0;j<co[i];j++){
					xmean+=ba[i][j].x;
					ymean+=ba[i][j].y;
				}
			}
			xmean*=normal;
			ymean*=normal;
			// 3. compute distance from last dissolved particle (meaning call this function when N-1 particle has been dissolved.
			xmean = (myx - xmean)*(myx - xmean);
			ymean = (myy - ymean)*(myy - ymean);
			// 4. store in the member data
			rp[nN] = sqrt(xmean+ymean);
			// 5. compute distance wrt origin and store
			xmean = myx*myx;
			ymean = myy*myy;
			r0p[nN] = sqrt(xmean+ymean);
		}
		else{
			// Store the distance from the particle to the reference ceiling height yref
			//printf("n_dissed %d nN %d, p.y %f, yref, %f, delta y %f\n", n_dissd, nN, p->y, yref, p->y - yref);
			rp[nN] = fabs(p->y-yref);
		}
		nN++;
	}
}
