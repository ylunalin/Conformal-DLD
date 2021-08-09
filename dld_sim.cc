#include "dld_sim.hh"

typedef std::complex<double> Comp;
dld_sim::dld_sim():
	n_part(0), n_dissd(0), max_swp_thr(5), mrsq(1.), mr(1.),abmr(1.),
	co(NULL), mem(NULL), wlpos(NULL),
	rando(NULL), buf(NULL), shape(NULL), ba(NULL)
	{}
	
/*  #################### Constructor and Setup functions ################### */
/** Constructs a dld_sim class, allocating a block structure to store the
 * particles positions and opening a diagnostic output file.
 * Good for all geometry, parameters controlled by the sim_params object.
 * \param[in] spars the pointer to a sim_params object that controls the type of simulation to initialize. */
dld_sim::dld_sim(sim_params *spars_p_):
	n_part(0), n_dissd(0), max_swp_thr(5), mrsq(1.), mr(1.),abmr(1.),
	rando(gsl_rng_alloc(gsl_rng_taus)), buf(new char[256]), spars_p(spars_p_)
	{
		get_sim_params();

		co=new int[mn];
		mem=new int[mn];
		ba=new particle*[mn];
		for(int i=0;i<mn;i++) {
			co[i]=0;
			if(!(spars_p->aggre)){
				ba[i]=new particle[init_region_memory];
				mem[i]=init_region_memory;
			}
			else{
				ba[i] = NULL;
				mem[i]=0;
			}
		}
		dspcell = double (nocell)/dx;
		int min[] = {0,0};
		int max[] = {128,128};
		wlpos = new int[nocell*nocell+1];
		wl_gen my_list(min, max, max_swp_thr, nocell);
		wl=my_list.worklist;
		// Store the starting points of each worklist in wlpos[]
		wl_pos();

		//hard coded threshold value, change this later
		/*
		Grid g=extract_grid();
		mytree=KDTree(g);
		thr.push_back(1);
		thr.push_back(2);
		thr.push_back(4);
		thr.push_back(8);
		thr.push_back(10);
		thr.push_back(13);
		thr.push_back(18);
		*/
	}


/** Extracts paramters from spars object and call set_rnd() and set_blc_size(...).
 * For flat geometry, it also sets some parameters in sim_params object, such as f_rate, target, burnin, etc.*/
void dld_sim::get_sim_params(){
	sim_params &spars = *spars_p;
	tot_nr = spars.nr;
	seed = spars.randseed;
	nocell = spars.nocell;
	rad = spars.rad;
	redge = rad;
	height = 0;
	tot_nrsp = 1.0 / (double) tot_nr;
	cell_normalizer = 1/double (nocell*nocell);
	filename = spars.sim_dir;

	set_rnd(seed);

	// adld stuff
	if(spars.simtype == spars.ADLD){
		if(!spars.adld_set) fatal_error("Initialization", PARAM_NOT_SET);
		set_blc_size(spars.params, spars.terms, spars.length, rad, spars.dtheta);
	}
	else if(spars.simtype == spars.CIRCLE){
		shape = new Comp[1];
		set_blc_size();
	}
	else if(spars.simtype == spars.FLAT1){
		if(!spars.aflat_set) fatal_error("Initialization", PARAM_NOT_SET);
		rad /= 2;
		redge = rad;
		height = spars.height;
		ul = height;
		cftol = spars.cftol;
		mr = -rad;
		abmr = 5.0;
		shape = new Comp[1];
		set_blc_size(spars.tilt);
		spars.set_output_params(twidth, nv, dx, true_tilt);
	}
	else if(spars.simtype == spars.FLAT2){
		if(!spars.nflat_set) fatal_error("Initialization", PARAM_NOT_SET);
		rad /= 2;
		redge = rad;
		height = spars.height;
		ul = height;
		mr = -rad;
		abmr = 5.0;
		shape = new Comp[1];
		set_blc_size(spars.siden);
		spars.set_output_params(twidth, nv, dx, true_tilt);
	}
}

/** Set vital parameters in the simulation: dx, nv, nh, epsilon, xsp. For circles.*/
void dld_sim::set_blc_size(){

	double tmp_r = rad;
	int i=0;
	while (tmp_r>min_blc_size){
		tmp_r /= 2.;
		i++;
	}
	nh = std::pow(2,i+1);
	if (height!=0) nv = int(ceil(height/tmp_r));
	else nv = nh;
	mn = nv*nh;
	dx = tmp_r;
	epsilon = 0.0001;
	xsp = 1/tmp_r;
	twidth = rad+redge;
	is_flat = false;
}

/** Set vital parameters in the simulation: dx, nv, nh, epsilon, xsp. For adlds.
 * \param[in] params the array of coefficients for Laurent series that defines the shape.
 * \param[in] terms the array of powers in the Laurent series.
 * \param[in] length the first dimension of params array.
 * \param[in] ur the size of the cluster.
 * \param[in] dtheta sampling rate around the unit circle. */
void dld_sim::set_blc_size(double ** params, int* terms, int length, double ur, double dtheta){
	int i=0;
	adld_len = int(2*pi/dtheta);
	//printf("adld_len = %d\n", length);
	shape = new Comp[adld_len];
	for(i=0;i<adld_len;i++) shape[i] = Comp(0,0);

	adld_outline(shape, params, terms, length, ur, dtheta);

	double rmax = -(rad*rad*10000);
	for(i=0; i<adld_len; i++){
		double shape_rsq = real(shape[i])*real(shape[i]) + imag(shape[i])* imag(shape[i]);
		if(shape_rsq > rmax) rmax = shape_rsq;
	}
	rad= std::ceil(sqrt(rmax));
	double tmp_r = rad;
	i=0;
	while (tmp_r>min_blc_size){
		tmp_r /= 2.;
		i++;
	}

	//i+1 because x-length is double the radius
	nh = std::pow(2,i+1);
	if (height!=0) nv = int(ceil(height/tmp_r));
	else nv = nh;
	mn = nv*nh;
	dx = tmp_r;
	epsilon = 0.0001;
	xsp = 1/tmp_r;
	twidth = rad+redge;
	is_flat = false;
}

/** Set vital parameters in the simulation: par_lim, dx, nv, nh, epsilon, xsp
 * true_tilt, redge, height, blclow, ylow.
 * For flat geometry, with vertical number of particle specified.
 * \param[in] n user provided particle numerator (in tan(theta) calculation).*/
void dld_sim::set_blc_size(int n){
	int d = int (sqrt(rad*rad*4 - n*n));
	ul = sqrt(d*d+n*n);
	redge = ( (int) (rad*2/ul)) * ul  - rad;
	height = int (std::ceil(height/ul)) *ul;
	if(height==0) height=ul;
	true_tilt = atan2(double(n), double(d));
	par_n = n;
	par_d = d;
	// Set block size
	double tmp_r = rad + redge;
	int i=0;
	while (tmp_r>min_blc_size){
		tmp_r /= 2.;
		i++;
	}

	//i because x-length is the radius
	nh = std::pow(2,i);

	// nv is multiplied by two to allow shifting the whole array of particle up when it's dissolved away partially.
	if (height!=0) nv = 2*int(ceil(height/tmp_r));
	else nv = nh;
	mn = nv*nh;
	dx = tmp_r;
	xsp = 1/tmp_r;
	epsilon =0.0001;

	// This puts ylow just outside of the maximum height, such at ylow-Y_DIS is at the max height
	blclow=nv-1;
	// This might lead to survey(dh, nelem) breaking, since mr is bounding the block of particles, and ylow is bounding the simulation domain
	// mr - ylow, which should give the total high difference, would be negative before any particle was dissolved
	ylow=nv*dx-rad +sfac*sqrt(2)*cos(pi/4- true_tilt); // this is the highest possible y (accounting the particle size)
	yref= round(nv*dx - rad + 0.5);
	twidth = rad+redge;
	is_flat = true;
}

/** Set vital parameters in the simulation: par_lim, dx, nv, nh, epsilon, xsp
 * true_tilt, redge, height, blclow, ylow.
 *  For flat geometry, with tilting angle specified.
 * \param[in] tilt user provided tilting angle.*/
void dld_sim::set_blc_size(double tilt){
	int d=1, n=0, par_lim=20;
	bool success = 0;
	if (tilt !=0){
		success = cont_frac(tilt, n, d, par_lim);
		while(!success){
			par_lim +=10;
			success = cont_frac(tilt, n, d, par_lim);
		}
		ul = sqrt(d*d+n*n);
		par_n = n;
		par_d = d;
		redge = ( (int) (rad*2/ul)) * ul  - rad;
		height = int (std::ceil(height/ul)) *ul;
		if(height==0) height=ul;
		if(redge<= -rad) {
			sprintf(buf,"Width is too small to accomodate tilt = %.0f estimate tolerance (%f) and unitlength %f\n", tilt, cftol, ul);
			fatal_error("set_blc_size()", buf);
		}
		true_tilt = atan2(double(n), double(d));
	}
	else{
		true_tilt = 0;
		ul = rad;

	}
	// Set block size
	double tmp_r = rad + redge;
	int i=0;
	while (tmp_r>min_blc_size){
		tmp_r /= 2.;
		i++;
	}

	//i because x-length is the radius
	nh = std::pow(2,i);

	// nv is multiplied by two to allow shifting the whole array of particle up when it's dissolved away partially.
	if (height!=0) nv = 2*int(ceil(height/tmp_r));
	else nv = nh;
	mn = nv*nh;
	dx = tmp_r;
	xsp = 1/tmp_r;
	epsilon =0.0001;

	// This puts ylow just outside of the maximum height, such at ylow-Y_DIS is at the max height
	blclow=nv-1;
	// This might lead to survey(dh, nelem) breaking, since mr is bounding the block of particles, and ylow is bounding the simulation domain
	// mr - ylow, which should give the total high difference, would be negative before any particle was dissolved
	ylow=nv*dx-rad +sfac*sqrt(2)*cos(pi/4- true_tilt);
	yref= round(nv*dx - rad + 0.5);
	twidth = rad+redge;
	is_flat = true;
}

/** The class destructor closes the diagnostic file and frees the dynamically allocated memory. */
dld_sim::~dld_sim() {
	if(wlpos!=NULL) delete [] wlpos;
	if(ba!=NULL){
		for(int i=mn-1;i>=0;i--) delete [] ba[i];
		delete [] ba;
	}
	if(mem!=NULL) delete [] mem;
	if(co!=NULL) delete [] co;
	if(shape!=NULL) delete [] shape;
	if(buf!=NULL) delete [] buf;
	if(rando!=NULL) gsl_rng_free(rando);

}

/* ####################  House Keeping ####################*/

/** Reset all the counters, and the storage array.
 * Used when multiple simulations are using the same dld_sim object.
 * \param[in] dla a flat to tell how much initial region memory to allocate
 * if true then initial memory is zero. Default value is false. */
void dld_sim::clean(bool dla){

	n_part=0;
	n_dissd=0;
	yref= round(nv*dx - rad + 0.5);

	// Reinitialize the simulation blocks
	for(int i=0; i<mn;i++){
		co[i]=0;
		if(!dla){
			mem[i]=init_region_memory;
			particle* nba=new particle[init_region_memory];
			delete [] ba[i];
			ba[i]=nba;
		}
		else{
			mem[i]=0;
			if(ba[i]!=NULL) delete [] ba[i];
			ba[i] = NULL;
		}

	}
	mrsq=1.0;
	if (is_flat) mr = -rad;
	else mr = 1;
}

/** Checks if there are duplicate particles in the system */
void dld_sim::check_dup(){

	for (int i=0; i<mn;i++){
		for (int j =0;j<co[i]; j++){
			particle * dummy = &ba[i][j];
			for (int k=0; k<co[i]; k++){
				if (k!=j){
					if(ba[i][k].x<dummy->x+epsilon && ba[i][k].x>dummy->x - epsilon && ba[i][k].y <dummy->y+epsilon && ba[i][k].y>dummy->y- epsilon) {
						sprintf(buf, "Found duplication particles at block %d, (%f %f)", i, ba[i][k].x, ba[i][k].y);
						fatal_error("check_dup()", buf);
					}
				}
			}
		}
        }
}

/** Take a particle from the simulation.
 * \param[in] doom the pointer to the particle to be taken. */
void dld_sim::take(particle* doom) {
	double x=doom->x;
	double y=doom->y;
	double ax = (-1)*rad, ay = (-1)*rad;
	// Compute which block the particle is within
	int nx=int((x-ax)*xsp);if(nx<0) nx=0;if(nx>=nh) nx=nh-1;
	int ny=int((y-ay)*xsp);if(ny<0) ny=0;if(ny>=nv) ny=nv-1;
	int s=nx+nh*ny;
	*doom = ba[s][--co[s]];
	n_part --;
}

/** Adds a particle to the simulation.
 * \param[in] x the x coordinate of the particle.
 * \param[in] y the y coordinate of the particle. */
particle * dld_sim::put(double x,double y) {
	double ax = (-1)*rad, ay = (-1)*rad;
	// Compute which block the particle is within
	int nx=int((x-ax)*xsp);if(nx<0) nx=0;if(nx>=nh) nx=nh-1;
	int ny=int((y-ay)*xsp);if(ny<0) ny=0;if(ny>=nv) ny=nv-1;
	int s=nx+nh*ny;
	//printf("put: nh = %d, nv=%d, nx=%d, ny=%d\n", nh, nv, nx, ny);
	if (s>=mn) {
		sprintf(buf,"Putting particle at x=%.2f, y=%.2f, nx=%d, ny=%d, s=%d.\nBlock index out of bound. Probably from shiftup().\n", x, y, nx, ny, s);
		fatal_error("put():", buf);
	}

	// Check there's enough memory for the particle,
	// and if not, extend the memory
	if(co[s]==mem[s]) add_region_memory(s);
	ba[s][co[s]++]=particle(x,y);
	n_part++;
	if(!is_flat){
		double rsq=sqrt(x*x+y*y)+1;
		// Update the maximum radius of the cluster
		if(rsq>mr) {
			mr=rsq;
			mrsq=mr*mr;
		}
	}
	else{
		double R = y+abmr;
		if (R>mr){
			mr = R;
			// hypotenuse squared
			mrsq = mr*mr + 4*rad*rad;
		}
	}
	return ba[s]+co[s]-1;
}

/** Initialize the cluster from file, return if the file was read successfully
 * \param[in] init_filename the full path and name to the file for initialization */
bool dld_sim::init_from_file(char* init_filename){
	std::fstream fn(init_filename, std::ios_base::in);
	double a,x=0,y=0;
	int counter=0;
	while(fn>>a){
		if(counter%2==0) x=a;
		else {
			y=a;
			put(x,y);
		}
		counter++;
	}
	return counter;
}

/** Consider a walker at (x,y)=(r,0). This samples the angle at which the walker first
 * touches the unit circle.
 * \param[in] r the x position of the walker. */
double dld_sim::bias_angle(double r) {
	return 2*std::atan((r+1)/(r-1)*std::tan((rnd()-0.5)*pi));
}

/** Consider a walker at (x,y). This maps the walker to the line y=mr
 * the x position is determined by uniform sampling from 0 to 2PI, then conformal map to real line
 * \param[in,out] x the x position of the walker.
 * \param[in,out] y the y position of the walker. */
void dld_sim::pbcremap(double &x, double &y){
	if(y>(mr+abmr)){
		double b = y-mr;
		double theta =  2*pi*rnd();
		Comp z_denom = Comp(1.0+cos(theta), sin(theta));
		Comp z_num   = Comp(b*sin(theta), b*(1.0-cos(theta)));
		Comp w= z_num / z_denom;
		double new_x = real(w)+x;

		//printf("pbcremap(): (%f, %f)\n", new_x, new_y);
		y = mr;
		if(new_x< -rad) new_x=fmod(new_x+rad, twidth)+redge;
		else if(new_x >=redge) new_x=fmod(new_x-redge, twidth) -rad;
		x=new_x;
        }
}

/** Calculates the total number of particles in the simulation. */
int dld_sim::total_particles() {
	int tot=*co;
	for(int i=1;i<mn;i++) tot+=co[i];
	return tot;
}

/** Check whether north, south, east, west neighbors exist
 * 0, 1 = north, south, 2, 3 = west, east.
 * \param[in] x x-coordinate of the current particle.
 * \param[in] y y-coordinate of the current particle.
 * \param[in] dir a integer representation of the inquired neighbor.*/
bool dld_sim::ex_neighbor(double x, double y, short dir){
	double sintheta = sin(true_tilt);
	double costheta = cos(true_tilt);
	double tol = 0.1, ngx, ngy;
	//north
	if(dir==0){
		ngx = x-sintheta;
		ngy = y+costheta;
		if(is_flat && ngx<-rad) ngx +=twidth;
	}
	//south
	else if (dir==1){
		ngx = x+sintheta;
		ngy = y-costheta;
		if(is_flat && ngx>=redge) ngx -=twidth;
	}
	//west
	else if (dir==2){
		ngx = x-costheta;
		ngy = y-sintheta;
		if(is_flat && ngx<-rad) ngx +=twidth;
	}
	//east
	else if (dir ==3){
		ngx = x+costheta;
		ngy = y+sintheta;
		if(is_flat && ngx>=redge) ngx -=twidth;
	}
	else fatal_error("ex_neighbor()", "dir must be 0, 1, 2, 3 for north, south, east, west");

	bool exist = 0;

	// highest possible is nv*dx -rad since the particles can be shifted up
	if (ngy>(-1)*rad && ngy<(nv*dx-rad)){
		double ax = (-1)*rad;
		double ay = ax;
		int nx=int((ngx-ax)*xsp);if(nx<0) nx=0;if(nx>=nh) nx=nh-1;
		int ny=int((ngy-ay)*xsp);if(ny<0) ny=0;if(ny>=nv) ny=nv-1;
		int s = nx+nh*ny;
		for (int i =0; i<co[s]; i++){
			double px = ba[s][i].x;
			double py = ba[s][i].y;
			if (px<ngx+tol && px>ngx-tol && py <ngy+tol && py>ngy-tol){
				exist=1;
				break;
			}
		}
	}
	return exist;
}

/** Sets the lowest y coordinate in the system and the corresponding lowest block.
 * \param[in] pp the pointer to the particle that is just dissolved. */
void dld_sim::set_lowest(particle *pp){
	double ay = (-1)*rad;
	double y = pp->y, x =pp->x;
	int ny=int((y-ay)*xsp);if(ny<0) ny=0;if(ny>=nv) ny=nv-1;
	if (ny<blclow) {
		blclow = ny;
	}
	if (y<ylow){
		ylow=y;
		xlow=x;
	}
}

/** Determines whether the cluster has dissolved enough to be shifted up by a unit height. */
bool dld_sim::diss_enough(){
	if(blclow>1) return false;
	else{
		double ay = (-1)*rad;
		double y;
		int ny;
		for(int i=mn-1;i>=0;i--){
			for (int j=0; j<co[i]; j++){
				y = ba[i][j].y;
				ny=int((y-ay)*xsp);if(ny<0) ny=0;if(ny>=nv) ny=nv-1;
				if((ny+1) > (int) nv/2) {
					return false;
				}
			}
		}
		return true;
	}
}

/** Shifts the entire cluster up by member data, height, a unit length determined by cont_frac().
 * \param[in] dupcheck a switch to check if there are duplicate particles.
 * \param[in] verbose a switch to turn on self-identifying and debug info. */
void dld_sim::shiftup(bool dupcheck, bool verbose){
	bool enough = diss_enough();
	if (!enough) {
		if(verbose) msg("shiftup()", "Have not dissolved enough.");
		return;
	}
	else{
		if (blclow < 1 ) fatal_error("shiftup()", "Might have dissolved too much, or blclow artifially set too low.");

		for (int i=mn-1; i>=0; i--){
			for (int j=0; j<co[i]; j++){
				double xx = ba[i][j].x;
				double yy = ba[i][j].y;
				put(xx, yy+height);
			}
			if (co[i]!=0){
				n_part -= co[i];
				mem[i]  = init_region_memory;
				co[i] = 0;
			}
		}
		ylow += height;
		yref += height;
		int ny = (int) ((ylow+rad)*xsp); if(ny<0) ny=0; if(ny>=nv) ny= nv-1;
		blclow=ny;

		flat(verbose);
		if(dupcheck) check_dup();
	}
}

/** Doubles the memory allocation in a simulation block.
 * \param[in] s the block to double. */
void dld_sim::add_region_memory(int s) {
	if(mem[s]==0) mem[s]=1;
	else mem[s]<<=1;
	if(mem[s]>=max_region_memory) {
		sprintf(buf, "Memory allocation exceeded in region %d, number of particles %d. But not fatal.\n",s, co[s]);
		msg("add_region_memory():",buf);
	}
	particle* nba=new particle[mem[s]];
	if (ba[s]!=NULL) {
		for(int i=0;i<co[s];i++) nba[i]=ba[s][i];
		delete [] ba[s];
	}
	ba[s]=nba;
}

/** Halves the memory allocation in a simulation block.
 * \param[in] s the block to halve. */
void dld_sim::halve_region_memory(int s) {
	// >1 because we want non-zero memory, since add_region_memory
	// does bit shifting, shifting zero still gives zero.
	if(mem[s]>1){
		mem[s]>>=1;
		particle* nba=new particle[mem[s]];
		for(int i=0;i<co[s];i++) nba[i]=ba[s][i];
		delete [] ba[s];
		ba[s]=nba;
	}
}

/** Use sim_params info to initialize a cluster, gives back the total number of particles;
 * cleans up all the particles afterward. This is used in output filename of active zone files.*/
void dld_sim::calculate_cluster_size(){
	sim_params &spars = *spars_p;
	bool hexa = spars.hexa, adld = (spars.simtype==1);
	if(hexa) printf("Nothing has been done. Hexagonal grid not implemented.\n");
	if(!adld) square(rad);
	else adld_cluster();
	int tb = total_particles();
	spars.set_output_params(tb);
	clean();
}

/** Saves the particle positions as a text file containing integer IDs,
 * positions, lengths, and rotations.
 * \param[in] suffix the string suffix to add to the filename. */
void dld_sim::write(const char* k) {
	// create the output filename and open the output file
	sprintf(buf,"%s/f.%s",filename,k);
	FILE *fp=safe_fopen(buf,"w");

	//output the particle positions to the file
	for(int s=0;s<mn;s++) {
		for(int q=0;q<co[s];q++) {
			particle &b=ba[s][q];
			fprintf(fp,"%g %g\n",b.x,b.y);
		}
	}

	// close the file
	fclose(fp);
}
/** Saves the particle positions as a text file containing integer IDs,
 * positions, lengths, and rotations.
 * \param[in] k the integer suffix to add to the filename. */
void dld_sim::write(int k) {
	// create the output filename and open the output file
	sprintf(buf,"%s/f.%d",filename,k);
	FILE *fp=safe_fopen(buf,"w");

	//output the particle positions to the file
	for(int s=0;s<mn;s++) {
		for(int q=0;q<co[s];q++) {
			particle &b=ba[s][q];
			fprintf(fp,"%g %g\n",b.x,b.y);
		}
	}

	// close the file
	fclose(fp);
}

/** Scans all the blocks and give the location of the last particle. */
particle dld_sim::last() {
	if(total_particles()==1){
		particle loner;
		for(int i=0; i<mn; i++){
			if(co[i]==1) {
				loner=ba[i][0];
				// take away the particle
				co[i]--;
				break;
			}
		}
		return loner;
	}
	else {
		fatal_error("last()", "More than just one particle left, but other parts of the code claim only 1 particle left. Check your code.");
		exit(0);
	}
}

/** Scan all the blocks and shrink the max radius if applicable. */
void dld_sim::shrink_r() {
	double rsq=0;

	if(is_flat) {
		rsq = -rad;
		bool update = false;
		for(int i=0;i<mn;i++){
			for(int j=0;j<co[i];j++){
				double y = ba[i][j].y;
				if(y>=rsq) {
					rsq = y;
					update = true;
				}
			}
		}
		if(update==true){
			mr = rsq+abmr;
			mrsq = mr*mr + 4*rad*rad;
		}
	}
	else{
		for(int i=0;i<mn;i++){
			// This bit frees up unused memory.
			if(mem[i]>1){
				int gold_fish=mem[i];
				while(co[i]<=(gold_fish>>=1)){
					halve_region_memory(i);
					if(gold_fish == 1) break;
				}
			}

			// This bit shrinks down the radius.
			for(int j=0;j<co[i];j++){
				double my_rsq;
				my_rsq=ba[i][j].x*ba[i][j].x+ba[i][j].y*ba[i][j].y;
				if(my_rsq>=rsq) rsq=my_rsq;

			}
		}
		mr=sqrt(rsq)+2.*sfac;
		mrsq=mr*mr;
	}
}


/* ####################  Creating the grids ################### */

/** Create a circular cluster of radius r on hexagonal grid. Particle radius = 0.5.
 * \param[in] r cluster radius. */
void dld_sim::hex(double r){
	int rmax = (int) r +1;
	// The side of the hexagonal cells is 1
	// h is the height from the center of the cell to the side
	double h=sqrt(3)*0.5;
	// Radius of the cluster in the units of h.
	int rh= (int)(r/h)+1;

	// j is the index of the rows
	for(int j=(-1)*rmax; j<=rmax; j++){
		// i is the index of the columns.
		// In this particular hexagonal lattice, all the points line up vertically,
		// spaced by h. While horizontally, spacing alternates.
		for(int i=(-1)*rh; i<=rh; i++){
			if(i%2==0){
				double x=(double) (i*h);
				double y=(double) j;
				if(x*x+y*y<=r*r) put(x,y);
			}
			else{
				double x=(double) (i*h);
				double y=(double) (j+0.5);
				if(x*x+y*y<=r*r) put(x,y);
			}
		}
	}
}


/** Create a circular cluster of radius r, on cartesian grid. Particle radius=0.5.
 * \param[in] r the radius of the cluster. */
void dld_sim::square(double r) {
	// Start a circular cluster by putting a seed at (0,0)
	int rmax = (int) r+1;
	for(int i=(-1)*rmax; i<=rmax;i++){
		for (int j=(-1)*rmax; j<=rmax;j++){
			double x = (double) i;
			double y = (double)j;
			if (x*x+y*y<=r*r) put(x, y);
		}
	}

}

/** Generates the outline of shapes in ADLD.
 * \param[out] result the complex number that represent the outline of the shape.
 * \param[in] params a 2D double pointers array specifying the conformal mapping coefficient of the shape.
 * \param[in] terms the pointer to an array specifying the powers in the Laurent series of the shape.
 * \param[in] length length of the 2D double points, row wise, i.e. number of terms in Laurent series.
 * \param[in] ur the size of the cluster.
 * \param[in] dtheta sampling rate around the unit circle. */
void dld_sim::adld_outline(Comp *&result, double ** params, int * terms, int length, double ur, double dtheta){
	Comp *mypars = new Comp[length];
	for(int i =0; i<length; i++){
		mypars[i] = Comp(params[i][0], params[i][1]);
		//printf("%f %f\n", real(mypars[i]), imag(mypars[i]));
	}
	int domain = int (2*pi/dtheta);
	for (int i =0; i<domain; i++){
		double theta = i*dtheta;
		Comp z = Comp(cos(theta), sin(theta));
		for (int j= 0; j<length; j++){
			result[i] += ur*mypars[j]* pow(z, terms[j]);
		}
	}
	delete [] mypars;
}

/** Creates cluster according to ADLD shape outline. */
void dld_sim::adld_cluster(){
	Comp * outline = shape;
	int outline_len = adld_len;

	// grad the x and y coordinate in 2D
	double * xcoord = new double[outline_len];
	double * ycoord = new double[outline_len];

	// find xmax and xmin
	double xlargest = -rad;
	double xsmallest = rad;
	double ylargest = xlargest, ysmallest = xsmallest;
	for (int i=0; i<outline_len; i++){
		xcoord[i] = real(outline[i]);
		ycoord[i] = imag(outline[i]);
		if(xcoord[i]>xlargest) {
			xlargest = xcoord[i];
		}
		if(xcoord[i]<xsmallest){
			xsmallest = xcoord[i];
		}
		if(ycoord[i]>ylargest) {
			ylargest = ycoord[i];
		}
		if(ycoord[i]<ysmallest){
			ysmallest = ycoord[i];
		}
		// check whether xcoord[i], ycoord[i] are integers, to be continued here
	}
	int xmax = int(std::ceil(xlargest));
	int xmin = int(std::floor(xsmallest));
	int ymax = int(std::ceil(ylargest));
	int ymin = int(std::floor(ysmallest));

	// create a vector to get store pair of verteces where ray crosses the segment
	std::vector<vertex_pair> crosses;

	// value of side =1 if vertex is to the left of ray; otherwise 0
	bool side = 0;
	// side value for next vertex
	bool next_side = 0;
	// true if current vertex is on the ray
	bool on_vertex = 0;
	// true if next vertex is on the ray
	bool next_on_vertex = 0;
	for (int x=xmax; x>=xmin; x--){
		//printf("x =%d\n", x);
		int vp_co = 0;
		// compare first point on the polygon to last point
		if(xcoord[0]<x) {
			side=1;
			on_vertex = 0;
		} // if to the left of the ray
		else if (xcoord[0]>x) {
			side = 0;
			on_vertex = 0;
		} // if to the right
		else {
			vertex_pair single_vp(xcoord[0], ycoord[0]);
			on_vertex =1 ;
			// check whether the ray is grazing this single vertex
			bool graze_check1 = 0;
			bool graze_check2 = 0;
			if (xcoord[1] < xcoord[0]) graze_check1 = 1;
			else if (xcoord[1] > xcoord[0]) graze_check1 =0;
			else fatal_error("adld_cluste()", "Encounters vertical lines\n");

			if (xcoord[outline_len-1] < xcoord[0]) graze_check2=1;
			else if( xcoord[outline_len-1]>xcoord[0]) graze_check2=0;
			else fatal_error("adld_cluster()", "Encounters vertical lines\n");

			// if two neighbors are on the same side, then the ray grazes current vertex
			if(graze_check1 == graze_check2) {
				single_vp.graze=1;
			}
			// if it doesn't graze it, then add it to the ylim list
			else{ vp_co++ ; }
			crosses.push_back(single_vp);
		} // if neither, ray go through vertex

		// since while loop only look at current vertex starting from 1, compared to the one before,
		// the loop will miss 0 comparing to the one before
		// here we test vertex 0 and vertex -1 (end)
		if(xcoord[outline_len-1] < x) {
			next_side=1;
			next_on_vertex =0 ;
		}
		else if (xcoord[outline_len -1] > x) {
			next_side =0;
			next_on_vertex = 0;
		}
		else {
			//crosses.push_back( vertex_pair(xcoord[outline_len-1], ycoord[outline_len-1]) );
			next_on_vertex = 1;
		}

		// on the in the case of neither point is on the ray, one is left and one is right
		if(side!=next_side && !(on_vertex | next_on_vertex) ){
			crosses.push_back( vertex_pair( xcoord[0], ycoord[0], xcoord[outline_len-1], ycoord[outline_len-1] ));
			vp_co++;
		}

		// starting at second vertex, compare to the one before
		int j = 1;
		while(j < outline_len){
			if (xcoord[j] < x ) {
				next_side = 1;
				next_on_vertex = 0;
			}
			else if (xcoord[j] > x) {
				next_side =0;
				next_on_vertex = 0;
			}
			else {
				next_on_vertex = 1;
				vertex_pair single_vp(xcoord[j], ycoord[j]);
				// check whether the ray is grazing this single vertex
				bool graze_check1 = 0;
				bool graze_check2 = 0;
				if (xcoord[j+1] < xcoord[j]) graze_check1 = 1;
				else if (xcoord[j+1] > xcoord[j]) graze_check1 =0;
				else fatal_error("adld_cluste()", "Encounters vertical lines");

				if (xcoord[j-1] < xcoord[j]) graze_check2=1;
				else if( xcoord[j-1]>xcoord[j]) graze_check2=0;
				else fatal_error("adld_cluste()", "Encounters vertical lines");

				// if two neighbors are on the same side, then the ray grazes current vertex
				if(graze_check1 == graze_check2) {
					single_vp.graze=1;
				}
				else{
					vp_co++;
				}
				crosses.push_back(single_vp);
			}


			if(next_side!=side && !(on_vertex | next_on_vertex)) {
				crosses.push_back( vertex_pair(xcoord[j], ycoord[j], xcoord[j-1], ycoord[j-1]) );
				vp_co++;
			}

			// set the current vertex
			side = next_side;
			on_vertex = next_on_vertex;
			j++;
			// if on the same side don't do anything
		}
		std::vector<vertex_pair>::iterator it;
		// if there are crossings
		if (vp_co!=0) {

			double * ylim = new double[vp_co];
			int dummy_index = 0;
			double dummy_y;
			for (it=crosses.begin(); it!=crosses.end(); it++){
				if ( (!(*it).single && !(*it).graze) || ((*it).single && !(*it).graze) ) {
					double lx1 = (*it).v1[0], ly1= (*it).v1[1], lx2= (*it).v2[0], ly2=(*it).v2[1];
					if ((*it).single) {
						 dummy_y = ly1;
					}
					else{
						double t = (double(x) - lx1)/(lx2 - lx1);
						dummy_y = ly1*(1-t)+ly2*t;
					}
					ylim[dummy_index] = dummy_y;
					dummy_index++;
				}
				else if ( (*it).single && (*it).graze ){ put( (*it).v1[0], (*it).v1[1] );  }
				else fatal_error("adld_cluste()", "A situation with a pair of verteces grazing the ray hasn't been considered, but detected anyway.");
			}
			// Gotta sort the ylim first!
			if (dummy_index >2 ){
				for (int p = 0; p<dummy_index; p++){
					double lips = ylim[p];
					int lips_pos = p;
					for (int q = p+1; q<dummy_index; q++){
						if(ylim[q] < lips){
							lips = ylim[q];
							lips_pos  = q;
						}
					}
					if(lips_pos!=p){
						double lips_tmp = ylim[p];
						ylim[p] = ylim[lips_pos];
						ylim[lips_pos] = lips_tmp;
					}
				}
				//for(int p=0; p<dummy_index; p++) printf("ylim[%d]= %f\n", p, ylim[p]);
			}

			dummy_index = 0;

			while(dummy_index<vp_co/2){
				double y_hi = ylim[dummy_index*2];
				double y_low = ylim[dummy_index*2+1];
				if(y_hi<y_low){
					double tmp = y_low;
					y_low = y_hi;
					y_hi = tmp;
				}
				//printf("yhi=%f, ylow=%f\n", yhi, ylow);
				double y = double(ymax);
				while( y>=ymin ){
					if( y>=y_low && y<= y_hi) put(double(x), y);
					y-=1.0;
				}
				dummy_index++;

			}
			delete [] ylim;
		}
		// if not, either ray goes through some verteces, or vector is empty
		else{
			for (it=crosses.begin(); it!=crosses.end(); it++)
				if((*it).single) put( (*it).v1[0], (*it).v1[1] );
		}
		// clear the vector
		crosses.clear();

	}
	delete [] xcoord;
	delete [] ycoord;
}


/** Use continued fraction to approximate tangent of a given angle
 * by fraction of integers within certain limit.
 * \param[in] deg the angle whose tangent is to be approximated.
 * \param[in,out] n numerator in the tangent fraction.
 * \param[in,out] d denomenator in the tangent fraction.
 * \param[in] par_lim limit on the number of n and p. */
bool dld_sim::cont_frac(double deg, int& n, int & d, int par_lim){
	bool success =0;
	if (deg > 45 || deg <=0){
		fatal_error("cont_frac()", "Provide an angle in (0, 45]");
	}
	// Numerical issue with using (int) to round tan(pi/4) to integer (round to 0, instead of 1)
	// so use this to treat angle=pi/4 as special case
	else if (deg == 45){
		n = par_lim;
		d = par_lim;
		success =1;
	}
	// All other values between (0,45) degs
	else{
		double t = tan(deg/180.0*pi);
		bool nxt = true;
		std::vector<int> frac;
		int first = (int) t;
		frac.push_back(first);
		t -= first;
		while(nxt && t!= 0.0) {
			// invert the remainder, and add the integer part to the list of integers
			double inv  = 1.0 / t;
			int tmp = (int) inv;
			t = inv - tmp;
			frac.push_back(tmp);
			// iterate through the list, add up the fractions, check that the numerator and the denomenator are both within limit
			std::vector<int>::iterator it;
			int num = 0, denom = 1;
			for(it=frac.end()-1;it!=frac.begin()-1;--it){
				num += (*it)*denom ;
				tmp = num;
				num = denom;
				denom = tmp;
				if (num > par_lim || denom > par_lim){
					nxt = false;
				}
			}
			// only update the numerator and denomenator when they are within limit
			// there's one extra inversion for the first element in the vector, so invert again
			if(nxt == true){ d = num; n = denom;}
		}
		if(n!=0){
			int mmin = par_lim / d;
			int moth  = par_lim / n;
			if (mmin>moth) mmin = moth;
			n = mmin*n;
			d = mmin*d;
			//printf("n=%d, d = %d\n", n, d);
			double ad_angle = fabs(atan2((double) n, (double) d)/pi*180 - deg);
			//printf("ad_angle = %f\n", ad_angle);
			if (ad_angle < cftol) {
				success= true;
			}
		}

	}
		return success;
}

/** Create a flat interface of npart particles per row. Particle radius =0.5.
 *  \param[in] verbose a switch to turn on debug info. */
void dld_sim::flat(bool verbose){
	// Start a flat interface, with certian tilt angle
	// a small number for inequality comparison
	double stop_y =rad;

	// everything starts at the lower left corner
	double llc = (-1) *rad;

	// a preliminary y value to stop adding particles
	if (height!=0) {
		stop_y = height-rad;

	}

	// if there's no tilt, just string particles on y=-rad, -rad+1, ...
	if(true_tilt==0){
		if(verbose) printf("flat(): true tilt = %f, width = %f, height =%f\n\n", true_tilt/pi*180, rad+redge, height);
		double y=llc;
		while(y<stop_y){
			double x=llc;

			//strictly small, because the real particle lives on the left boundary,
			//the image of it lives on the right boundary
			while (x<rad){
				put(x,y);
				x+=1.0;
			}
			y+=1.0;
		}
	}

	// If tan(theta) isn't zero, i.e. if there's a non-zero intersection angle
	// populate the left boundary first , then lower boundary
	else {
		double x = llc;
		double y = llc;
		// updating the coordinate of the right edge
		double sintheta = par_n/ul;
		double costheta = par_d/ul;

		// subtract the right length for it to match, see my notebook for illustration
		stop_y = height - rad;
		double stop_x = redge;
		if(verbose) printf("flat(): Adj tilt = %f, unitlengt= %f, adj width = %f, adj height =%f\n\n", true_tilt/pi*180, ul, rad+redge, height);

		// stringing particles along the gridlines intersecting left boundary
		while(y<=stop_y){

			// find the first particle whose x coord is larger than -rad on the next grid line above current line
			// shifting to the next gridline above to do so:
			double tmpx = x-sintheta;
			double tmpy = y+costheta;

			// find first particle whose x > -rad
			while(tmpx>llc+costheta){
				tmpx-=costheta;
				tmpy-=sintheta;
			}
			while(tmpx<llc-epsilon){
				tmpx+=costheta;
				tmpy+=sintheta;
			}

			// put the rest of the particles on current line
			// since we add the value first then put down particle, the upper bound is one step size smaller than the bounds
			// strictly small because last particle needs to live on the left edge
			while(x<stop_x-epsilon &&  y <stop_y-epsilon){
				put(x,y);
				x+=costheta;
				y+=sintheta;
			}

			x = tmpx; y = tmpy;
		}

		// reset to lower left corner llc
		x =llc; y=llc;

		// stringing particle along the gridlines intersecting lower boundary
		// main different from above is that here, I start with previous starting point, not current one
		// so the same code is actually finding current starting point.
		while(x<stop_x){
			double tmpx = x+sintheta;
			double tmpy = y-costheta;
			while(tmpy>llc+sintheta){
				tmpx -=costheta;
				tmpy -=sintheta;
			}
			while(tmpy<llc-epsilon){
				tmpx+=costheta;
				tmpy+=sintheta;
			}
			x = tmpx; y = tmpy;

			// put the rest of the particles on current line
			while(x<stop_x-epsilon && y<stop_y-epsilon){
				put(x,y);
				x+=costheta;
				y+=sintheta;
			}
			x = tmpx; y = tmpy;
		}


	}
}

/* ####################  Dissolution & Aggregation ################### */

/** Carries out a disrete dissolution simulation for a number of frames. */
void dld_sim::dissolve(){
	const char * FUNC_NAME = "dissolve()";
	sim_params &spars = *spars_p;
	// Get parameters from sim_params object
	bool adld=(spars.simtype==1), hexa = spars.hexa, w = spars.write, verbose = spars.verbose, dupcheck = spars.dupcheck;
	int f_rate = spars.f_rate, total_target = spars.target, targ=0, k=0, act_n=0;
	double tot_time=0, com_time;

	// Initialize the cluster
	if(!adld){
		if(hexa && !is_flat) hex(rad);
		else if (!hexa && is_flat) flat(verbose);
		else if (!hexa && !is_flat) square(rad);
		else fatal_error(FUNC_NAME, HEX_ERR);
	}
	else adld_cluster();

	// Write initial frame
	if(w) write(0);

	if (!is_flat){
		int tb_static = total_particles(), tb = tb_static;
		printf("%s: #Frame 0 [%d particles, max rad=%.6g]\n",FUNC_NAME, tb,mr);
		tot_time = wtime();
		while(tb>1) {
			k++;
			if(k>1) shrink_r();
			targ = tb-f_rate;
			if (targ <=0 ) targ =1;

			com_time=wtime();
			while(tb>targ) {
				particle * p = minus_particle();
				active_zone_dist(act_n, p);
				take(p);
				tb--;
				n_dissd++;
			}
			com_time = wtime()-com_time;

			// Output the fields
			if(w) write(k);

			// Print diagnostic information
			printf("%s: Frame %d [%d particles, max rad=%.6g, com=%.6g s]\n",FUNC_NAME,k,tb,mr,com_time);
		}
		tot_time = wtime() - tot_time;
		printf("%s: Use total time: %g s.\n",FUNC_NAME,tot_time);
	}// End of (!is_flat)

	else {
		targ = (int) twidth;
		// we massage the parameter a little to hasten the while loop
		// i.e. we don't want to do mod every particle
		int rate_ratio, eaten_sofar=0;
		if(targ>f_rate) {
			f_rate = targ;
			rate_ratio = 1;
		}
		else{
			rate_ratio = f_rate / targ;
		}
		// smallrate is what we will use in inner while loop
		printf("Flat %s: #Frame %d [%d particles, max height=%.6g, min height=%.6g, adjusted output rate is once per %d particles]\n", FUNC_NAME, k, n_dissd, mr, ylow, rate_ratio*targ);
		while(n_dissd<total_target){
			k++;
			// Since we have require that f_rate >= targ
			// we can loop through every targ # particles have been dissolved
			// and shrink and shift accordingly.
			// After rate_ratio * targ have been done, we save output
			for(int i=0;i<rate_ratio;i++){

				eaten_sofar = 0;
				com_time=wtime();
				while(eaten_sofar<targ){
					particle * p = minus_particle();
					active_zone_dist(act_n, p);
					take(p);
					n_dissd++;
					eaten_sofar++;
				}
				com_time = wtime()-com_time;
				tot_time+=com_time;

				shrink_r();
				shiftup(dupcheck,verbose);
			}
			printf("FLAT %s: Frame %d [%d particles, max height=%.6g, min height = %.6g, tot compute time=%.6g s]\n",FUNC_NAME,k,n_dissd,mr,ylow,tot_time);
			if(w) write(k);

		}
		printf("FLAT %s: Dissolved pars: %d\nUse total time: %.2g s.\n",FUNC_NAME, n_dissd, tot_time);
	}// End of is_flat
}

/** Fast dissolve() for many realizations on a flat interface with some tilting angle.
 * \param[in,out] x the x coordinate of the collapse point.
 * \param[in,out] y the y coordinate of the collapse point. */
void dld_sim::rapid_dissolve(double &x, double &y){
	sim_params &spars = *spars_p;
	// If we want to get asymptotic interface width on flat geometry
	if(is_flat){
		const char * FUNC_NAME = "Flat rapid_dissolve()";

		if(!spars.aflat_set && !spars.nflat_set) fatal_error(FUNC_NAME, PARAM_NOT_SET);

		bool verbose = spars.verbose, dupcheck = spars.dupcheck, lv=spars.local_verbose;
		int NELEM = spars.f_len, num_frame = spars.num_frame, npar = spars.target, burnin = spars.burnin, f_rate = spars.f_rate;
		double dh = spars.dh, **f = spars.f;

		int targ=(int) twidth, random_amount, checkpt;
		int k = 0, nelem;
		int rate_ratio=1;
		if(f_rate>targ) rate_ratio = f_rate/targ;
		int large_targ = targ*rate_ratio;
		double time=0, burn_time=0, ratiosp = 1.0 / (double) rate_ratio;

		// initialize the array
		flat(verbose);

		// internal data storagee
		double * int_data = new double[NELEM];
		for (int i=0;i<NELEM;i++) int_data[i]=0;

		// run initial survey to get initial width of the interface
		double * tmpf_init = survey(dh, nelem, true);
		if(tmpf_init!=NULL){
			if (nelem > NELEM) fatal_error(FUNC_NAME, SURVEY_ERR);
			for(int i=0; i<nelem; i++) {
				f[0][i] += tmpf_init[i];
			}
			delete [] tmpf_init;
		}

		checkpt = int( burnin/f_rate/10);
		if (checkpt<1) checkpt=1;
		// First dissolve a whole bunch to reach steady state
		k=0;
		while(n_dissd<burnin){
			k++;
			for(int j=0; j<rate_ratio; j++){
				int inc_t = wtime();
				for(int i=0; i<targ;i++) {
					particle *p = minus_particle();
					take(p);
				}
				inc_t = wtime() - inc_t;
				burn_time += inc_t;
				shrink_r();
				shiftup(dupcheck,verbose);
			}
			n_dissd += large_targ;
		}
		if(lv) printf("%s: Done burnin %d particles.\n", FUNC_NAME, burnin);

		// reusing the variable burn_checkpt is now number of particle to dissolve
		// burn_diss is a local counter
		int act_n = 0;
		checkpt = (int) ((npar-burnin)/f_rate/10);
		// keep dissolving until reaching target number of dissolved particles
		k=0;
		while(n_dissd < npar) {
			k++;
			if(k>=num_frame) fatal_error(FUNC_NAME, "output frame out of bound.");
			// Fixing nearest neighbor finding to be sweep(), not using KD tree
			// Periodic boundary condition is also fixed to be True

			for(int j=0;j<rate_ratio;j++){
				// Dissolve a random fraction of a layer
				random_amount = (int) rnd()*targ;
				double inc_t = wtime();
				for(int i=0;i<random_amount; i++) {
					particle * p = minus_particle();
					active_zone_dist(act_n, p);
					take(p);
					n_dissd++;
				}
				inc_t = wtime()-inc_t;
				time += inc_t;

				// Survey the height, record data to the 2D data array every (output_rate) particle,
				// including the first frame, i.e. when no particle has been dissolved
				double * tmpf = survey(dh, nelem);
				if(tmpf!=NULL){
					if (nelem > NELEM) fatal_error(FUNC_NAME, SURVEY_ERR);
					for(int i=0; i<nelem; i++) {
						int_data[i] += tmpf[i];
					}
					delete [] tmpf;
				}

				// continue to dissolve until one layer is dissolved away
				random_amount = targ-random_amount;
				inc_t = wtime();
				for(int i=0; i<random_amount; i++) {
					particle *p = minus_particle();
					active_zone_dist(act_n, p);
					take(p);
					n_dissd++;
				}
				inc_t = wtime() - inc_t;
				time += inc_t;

				// shrink down the bounding line, where random walkers are generated
				// and prevent the dissovling front breaching through to the bottom by shifting things up
				shrink_r();
				shiftup(dupcheck,verbose);
			}

			// Move the randomly collected cdf data to f
			// Assuming that data had reached steady state.
			for(int i=0;i<NELEM;i++){
				// each frame we have the averaged stats from after each particle is dissolved.
				f[k][i] += int_data[i]*ratiosp;
				int_data[i] = 0;
			}

			if (lv && k%checkpt==0) printf("%s: Frame %d [Dissolved pars %d out of %d]\n", FUNC_NAME, k, n_dissd, npar);

			if(blclow==0){
				int tb = total_particles();
				if (tb<=10) fatal_error(FUNC_NAME, TOO_FEW);
			}
		}

		printf("\n%s\n%s: Done for tilt=%.2f. Speed = %.1f p/s (Auxiliary not counted).\nTotal dissolved particles = %d, total time = %.1f secs.\nBurn particles = %d, burn time = %.1f secs.\n%s\n", \
			DIV_LINE, FUNC_NAME, true_tilt/pi*180, double(n_dissd)/time, n_dissd, time, burnin, burn_time, DIV_LINE);
		delete [] int_data;
	}// End of is_flat
	// If we are dissolving radial type geometry and want to know the collapse point.
	else{
		bool hexa = spars.hexa, adld = (spars.simtype==1), local_verbose = spars.local_verbose;
		bool output_set = spars.output_set;
		const char * FUNC_NAME = "Circular rapid_diss()";
		if(hexa) fatal_error(FUNC_NAME, BIN_MAT_ERR);
		if(!adld) square(rad);
		else adld_cluster();

		int tb_static = total_particles(),tb=total_particles(),targ=1, grid_edge = (int) rad;
		int grid_dim = 2*grid_edge +1, k=0, act_n=0;
		if(local_verbose) printf("%s: Cluster initially has %d particles.\n",FUNC_NAME, tb);
		double time = 0;

		// Do the first frame.
		targ = (int) (tb_static*0.8);
		while(tb>targ){
			particle *p = minus_particle();
			active_zone_dist(act_n, p);
			take(p);
			tb--;
			n_dissd++;
		}

		k=1;
		while(true) {
			k++;

			if(output_set && k<=5) {
				for(int s=0; s<mn; s++){
					for(int p=0;p<co[s];p++){
						int ind_x =(int) ba[s][p].x;
						int ind_y =(int) ba[s][p].y;
						int grid_pos = grid_dim*(grid_edge-ind_y)+(grid_edge-ind_x);
						//why only record from frame 2 to 5?
						spars.f[k-2][grid_pos]+=tot_nrsp;
					}
				}
			}

			if(tb==1) {
				// When there's only one particle left, record the position
				// and realization id of the collapse point.
				particle loner = last();
				// last() function takes the particle away, not recordable
				x=loner.x;
				y=loner.y;
				if(local_verbose) printf("%s: # %d frames to reach collapse point (%f, %f).\n", FUNC_NAME, k, x, y);
				break;
			}

			// At culster size 80%, 60%, 40%, 20%, shrink the max radius to speed up random walk
			// At these frames, increment by 1/nr at each lattice point that has a particle
			if(k<5) targ=(int) tb_static*(1-0.2*k);
			else targ = (int) tb*0.8;
			shrink_r();

			double inc_t = wtime();
			while(tb>targ) {
				particle *p =minus_particle();
				active_zone_dist(act_n, p);
				take(p);
				tb--;
				n_dissd++;
			}
			inc_t = wtime() - inc_t;
			time+=inc_t;
		}
		if(local_verbose) printf("\n%s\n%s: Done for cluster of %d particles. Speed = %.1f p/s (Auxiliary not counted).\nTotal dissolved particles = %d, total time = %.1f secs.\n%s\n", \
			DIV_LINE, FUNC_NAME, tb_static, double(n_dissd)/time, n_dissd, time, DIV_LINE);
	}
}

/** Simulates a process where dla and dld are present simultaneous,
 * but the ratio of probabilities of dissolution to aggregation can be adjusted.
 * dl2a stands for diffusion limited aggregation and annihilation. */
int dld_sim::dl2a(){
	sim_params &spars = *spars_p;
	const char * FUNC_NAME = "dl2a()";
	// TODO how to make this more general, i.e. aggregation on arbitrary seed geometry
	msg(FUNC_NAME, "Radius of Gyration calculation only makes sense for radial geometry!");
	if(!spars.dl2a_set) fatal_error(FUNC_NAME, PARAM_NOT_SET);
	char* init_filename = spars.init_filename;
	int start_fr = spars.start_fr, end_fr = spars.end_fr, npar = spars.target;
	double p_aggre = spars.p_aggre, *data = spars.rg_data;
	// local_verbose is true means the function is not used in streamline fashion
	bool w = spars.write, lv=spars.local_verbose;

	double rg;
	// If there isn't anything file to initize from, start from scratch
	if(init_filename==NULL) start_fr=0;

	// initialize the cluster by a seed at origin
	if(start_fr == 0) put(0,0);

	// or we can initialize from previous simulation
	else{
		bool success = init_from_file(init_filename);
		if(!success){
			msg(FUNC_NAME, INIT_FILE_ERR);
			start_fr=0;
			put(0,0);
		}
	}

	// write the initial file
	if(w) write(start_fr);

	rg=rgyration();
	int tb=total_particles(),targ,k, checkpt=int(std::ceil( (double) (end_fr-start_fr)/10));

	*(data+2*start_fr)=rg;
	*(data+2*start_fr+1)=tb;

	printf("%s: # Initial frame %d [%d particles, max rad=%.6g, radius of gyration=%.6g.]\n", FUNC_NAME, start_fr, tb,mr, rg);

	for(k=start_fr+1;k<end_fr;k++) {
		targ=int(pow(k+1,1.71));
		// if particle limit isn't provided, then set it to be target number of particles
		while(mr<=rad && tb<targ && tb>0 && tb<npar){
			double prob=rnd();
			if(prob>p_aggre){
				// TODO Active zone to be implemented here
				particle * p = minus_particle();
				take(p);
				tb--;
			}

			else{
				add_particle();
				tb++;
			}
		}

		// print diagnostic and write to file
		if(lv || k%checkpt==0 || k== end_fr-1)
			printf("%s: # Frame %d [%d particles, max rad=%.6g, radius of gyration=%.6g.]\n",FUNC_NAME, k,tb,mr, rg);
		if(w) write(k);
		rg = rgyration();
		*(data+2*k)=rg;
		*(data+2*k+1)=tb;

		// break if all is dissolved
		if(tb==0) {
			msg(FUNC_NAME, "All particle annihilated.\n");
			return k+1;
		}

		// break out the for loop if maximum radius or maximum number of particles are reached
		if (tb>=npar || mr>rad) {
			printf("%s: One of the two max limits reached. # par=%d, mr=%.2f\n", FUNC_NAME, tb, mr);
			// if the file hasn't been written along the way, and if the function is being used in a slow way
			// write the last file out
			if(!w && lv) write(k);
			// return the total number of frames
			return k+1;
		}
	}
	return k;
}

/** A test function for kdtree search. */
void dld_sim::test_search(bool swp){
	double largest = 5.76*rad*rad;
	double tree_threshold = 25*dx*dx;
	double rsq=0, sec_d=0., dsq=0., rot_angle=0.;
	double theta=2*pi*rnd(),x=(rad)*cos(theta),y=(rad)*sin(theta);
	particle *target = 0;
	while (target==NULL){
		if(swp){
			sec_d=dsq=largest;

			// When pbc is implemented, don't forget to add pbc in the argument list
			target =sweep(x, y, dsq, sec_d, rot_angle, false);
		}
		else{
			sec_d=dsq=tree_threshold;
			target = treesearch(x, y, mytree, dsq, sec_d, rot_angle);
		}
		if(target==NULL){
			if (dsq < tree_threshold) dsq=sqrt(dsq)-sfac*3.;
			else dsq=double (max_swp_thr) - sfac*3.;
			theta=2*pi*rnd();
			x+=dsq*cos(theta);
			y+=dsq*sin(theta);
			rsq=x*x+y*y;
		}

		if(rsq>1.44*(rad*rad)) {
			rsq=sqrt(rsq)/rad;
			theta=arg(x,y);
			theta+=bias_angle(rsq);
			x=(rad)*cos(theta);
			y=(rad)*sin(theta);
		}
	}
}

/** Simulates a random walker being introduced far away from the cluster, and
 * come into contact with it, and annihilate a particle in the cluster.
 * There are 3 ways this can be accomplished:
 * 1) sweeping using a worklist;
 * 2) by searching the quatree built on the grid of the simulation;
 * 3) by the most basic way to executing uniform random walk to a circle. */
particle * dld_sim::minus_particle(){
	sim_params &spars = *spars_p;
	// some variables needed to do this
	const char* FUNC_NAME = "minus_particle()";
	bool swp = spars.sweep, tree = spars.tree, verbose = spars.verbose, eaten=false;
	double largest = 5.76*mrsq, tree_threshold = 25*dx*dx, rsq=0, sec_d=0., dsq=0., rot_angle=0., loc_time=0;
	Comp pos_walker = Comp(0,0);
	particle *target = 0;
	// Initialize the walker position to be uniformly distributed on a
	// circle surrounding the cluster
	double theta=2*pi*rnd(),x=mr*cos(theta),y=mr*sin(theta);
	// or if if we are dealing with flat geometry, on a line
	if (is_flat) {
		x = rnd()*twidth - rad;
		y = mr;
	}

	// walk the random walker until it touches something
	while(!eaten) {
		// both the tree and the regular method guarantee second nearest neighbor
		sec_d = largest;
		dsq = largest;
		if(verbose) printf("%s: Random walker at (%g, %g).\n",FUNC_NAME, x,y);

		// if the time it takes to do the previous round is large
		// shrink down the bounding circle
		if (loc_time > 0.01)  shrink_r();
		loc_time = wtime();

		if(std::isnan(x) or std::isnan(y)) fatal_error(FUNC_NAME, LOC_NAN);


		// Search for nearest neighbor
		// sweep takes priority
		if(swp) target =sweep(x, y, dsq, sec_d, rot_angle, verbose);
		// then tree search
		else if(tree){
			if(is_flat) fatal_error(FUNC_NAME, TREE_NO_PBC);
			sec_d = dsq = tree_threshold;
			target = treesearch(x, y, mytree, dsq, sec_d, rot_angle);
		}
		else if (!swp && !tree){
			if(is_flat) fatal_error(FUNC_NAME, MINDIST_NO_PBC);
			target = min_distance(x,y,dsq,sec_d,rot_angle,verbose);
		}

		if(target!=NULL){
			bool inside, on_edge, outside, far_away;
			pos_of_rw(dsq, inside, on_edge, outside, far_away);

			// If the distance from random walker (rw) to
			// nearest neighbor is < 1, that means rw is inside the particle.
			if(inside) {
				printf("%s: Random walker at (%g, %g), nearest neighbor at (%g, %g), sec_d=%g, rot_angle=%g\n", FUNC_NAME, x, y, target->x, target->y, sec_d, rot_angle);
				write("minus_particle_fail");
				fatal_error(FUNC_NAME, INSIDE_ERR);
			}

			// If the distance from rw to nearest neighbor is within epsilon, exactly 1, annihilate!
			if (on_edge){
				eaten=true;
				if(verbose) printf("%s: Annihilate particle at (%g, %g).\n", FUNC_NAME, target->x, target->y);
				set_lowest(target);
			}

			// If the distance is near enough and the second nearest neighbor is far enough
			// we use pacman to take a walk
			if (outside){

				if(sec_d+1 > sqrt(dsq) + epsilon){
					if(verbose) printf("%s: Use pacman with sqrt(dsq)=%g, sec_d=%g, rot_ang=%g.\n", FUNC_NAME, sqrt(dsq), sec_d, rot_angle);
					// pacman will take care of pbc
					pos_walker = pacman(dsq, sec_d, rot_angle, target, eaten);
					x = real(pos_walker);
					y = imag(pos_walker);

					if(eaten){
						if(verbose) printf("%s: Annihilate particle, eaten by pacman, at (%g, %g).\n", FUNC_NAME, target->x, target->y);
						set_lowest(target);
					}

					else {

						if (is_flat){
							// first calculate distance between walker and either left or right edge
							// then fmod to wrap around because of periodic boundary condition
							if(x<-rad) x=fmod(x+rad, twidth)+redge;
							else if(x>=redge) x=fmod(x-redge, twidth)-rad;
							// in flat geometry, the distance we care for remapping is the y distance
							rsq = y*y;
						}

						else rsq = x*x+y*y;
						if(verbose) printf("%s: Random walker take a walk to (%g, %g).\n", FUNC_NAME, x,y);
					}
				}// end of if(sec_d+1 > sqrt(dsq) + epsilon)

				else {
					if(verbose) msg(FUNC_NAME, EQUI_DIST);
					// Take a step to the edge of the empty circle
					dsq=sqrt(dsq)-2*sfac-epsilon;
					theta=2*pi*rnd();
					x+=dsq*cos(theta);
					y+=dsq*sin(theta);
					rsq=x*x+y*y;
					if (is_flat) {
						rsq = y*y;
						if(x<-rad) x=fmod(x+rad, twidth)+redge;
						else if(x>=redge) x=fmod(x-redge, twidth)-rad;
					}
				}
			}// end of pacman

			// If the distance is quite far, i.e. dsq>9. we use jump to a circle
			if(far_away){
				if(verbose) msg(FUNC_NAME, FAR);
				// Take a step to the edge of the empty circle
				dsq=sqrt(dsq)-3.*sfac;
				theta=2*pi*rnd();
				x+=dsq*cos(theta);
				y+=dsq*sin(theta);
				rsq=x*x+y*y;
				if (is_flat) {
					if(x<-rad) x=fmod(x+rad, twidth)+redge;
					else if(x>=redge) x=fmod(x-redge, twidth)-rad;
					rsq = y*y;
				}
			}
		} //end of if(target!=NULL)

		// If (target == NULL), then sweep didn't find either the 2nd neighbor, or both neighbors,  within threshold, just take a random jump up to the threshold
		else {

			if(verbose) msg(FUNC_NAME, NNB_SEARCH_FAIL);
			//Need to distinguish between the case when 1 neighbor is found, and when neither is found
				//double dummydsq =sqrt(double(*(thr.end()-1)))*dx- 3.*sfac;
				double dummydsq = double (max_swp_thr) * dx - 3. *sfac;
				dsq = sqrt(dsq)-3.*sfac;
				if(dummydsq <dsq) dsq = dummydsq;

				// Take a step to the edge of the empty circle
				theta=2*pi*rnd();
				x+=dsq*cos(theta);
				y+=dsq*sin(theta);
				rsq=x*x+y*y;
				if (is_flat) {
					rsq = y*y ;
					if(x<-rad) x=fmod(x+rad, twidth)+redge;
					else if(x>=redge) x=fmod(x-redge, twidth)-rad;
				}
		} // end of discussion on target pointer

		// If the walker is outside the bounding circle of the cluster,
		// then remap it back onto the circle
		if (!is_flat){

			if(rsq>1.44*(mrsq)) {
				if(verbose) {
					msg(FUNC_NAME, TOO_FAR);
					printf("From (%g %g) to ", x, y);
				}
				rsq=sqrt(rsq)/mr;
				theta=arg(x,y);
				theta+=bias_angle(rsq);
				x=(mr)*cos(theta);
				y=(mr)*sin(theta);
				if(verbose) printf(" (%g, %g)\n",x,y);
			}
		}
		else pbcremap(x, y);

		loc_time = wtime() - loc_time;
	}// end of while loop

	return target;
}

/** Simulates a random walker being introduced far away from the cluster, and
 * come into contact with it, and adds a particle to the cluster.
 * It uses the sweeping routine to search for nearest neighbor. */
particle * dld_sim::add_particle(){
	const char * FUNC_NAME = "add_particle()";
	bool verbose = spars_p -> verbose;
	double largest=5.76*mrsq;
	double rsq=0., sec_d=0., dsq=0., rot_angle = 0.;
	bool stick=false;
	Comp pos_walker = Comp(0,0);
	particle *target = 0;
	// initialize random walker
	double theta = 2*pi*rnd(), x = mr*cos(theta), y = mr*sin(theta);
	if (is_flat) {
		x = rnd()*twidth-rad;
		y = mr;
	}
	while(!stick){
		// check random walker position, for debug purpose
		if(std::isnan(x) || std::isnan(y)) fatal_error(FUNC_NAME, LOC_NAN);

		// use sweep to find nearest and second nearest neighbor
		sec_d = largest;
		dsq = largest;
		target = sweep(x, y, dsq, sec_d, rot_angle, verbose);
		// check sweep result
		if(target!=NULL){
			bool inside, on_edge, outside, far_away;
			pos_of_rw(dsq, inside, on_edge, outside, far_away);
			// if random walker is inside a particle, something went wrong
			if(inside){
				write("add_particle_fail");
				printf("%s: Random walker at (%g, %g), nearest neighbor at (%g, %g), sec_d=%g, rot_angle=%g\n", FUNC_NAME, x, y, target->x, target->y, sec_d, rot_angle);
				fatal_error(FUNC_NAME, INSIDE_ERR);
			}
			// if random walker is distance 1 away from a particle, it's tangent with that particle
			if(on_edge){
				stick=true;
				if(verbose) printf("%s: Aggregated to particle at (%g, %g).\n", FUNC_NAME, target->x, target->y);
			}
			// if random walker is 3 particle away, use pacman
			if(outside){
				if(sec_d+1 > sqrt(dsq)+epsilon) {
					if(verbose) printf("%s: Using pacman with dsq= %g, sec_d=%g, rot_ang=%g.\n", FUNC_NAME, sqrt(dsq),sec_d,rot_angle);
					pos_walker = pacman(dsq, sec_d, rot_angle, target, stick);
					x = real(pos_walker);
					y = imag(pos_walker);
					if(stick){
						if(verbose) printf("%s: Aggregated particle at (%g, %g).\n", FUNC_NAME, target->x, target->y);
					}
					else {
						if (is_flat) {
							if(x<-rad) x=fmod(x+rad, twidth)+redge;
							else if(x>=redge) x=fmod(x-redge, twidth)-rad;
							rsq = y*y;
						}
						else rsq = x*x+y*y;
						if(verbose) printf("%s: Random walker took a walk to (%g, %g).\n",FUNC_NAME,x,y);
					}
				}
				else {
					if(verbose) msg(FUNC_NAME, EQUI_DIST);
					dsq=sqrt(dsq)-2*sfac-epsilon;
					theta=2*pi*rnd();x+=dsq*cos(theta);y+=dsq*sin(theta);rsq=x*x+y*y;
					if (is_flat) {
						rsq = y*y;
							if(x<-rad) x=fmod(x+rad, twidth)+redge;
							else if(x>=redge) x=fmod(x-redge, twidth)-rad;
						}

				}
			}
			// if random walker is too far, just use circular random walk
			if(far_away){
				if(verbose) msg(FUNC_NAME, FAR);
				// Take a step to the edge of the empty circle
				dsq=sqrt(dsq)-3.*sfac;
				theta=2*pi*rnd();
				x+=dsq*cos(theta);
				y+=dsq*sin(theta);
				rsq=x*x+y*y;
			}
		}// if target is not NULL
		else{
			if(verbose) msg(FUNC_NAME, NNB_SEARCH_FAIL);
			double dummydsq = double (max_swp_thr) * dx - 3. *sfac;
			dsq = sqrt(dsq)-3.*sfac;
			//if the distance to nearest neighbor is larger than maximum sweep distance,
			//that means sweep didn't find nearest neighbor
			if(dummydsq <dsq) dsq = dummydsq;
			theta=2*pi*rnd();
			x+=dsq*cos(theta);
			y+=dsq*sin(theta);
			rsq=x*x+y*y;
			if (is_flat) {
				rsq = y*y;
				if(x<-rad) x=fmod(x+rad, twidth)+redge;
				else if(x>=redge) x=fmod(x-redge, twidth)-rad;
			}

		}
		// If the walker is outside the bounding circle of the cluster,
		// then remap it back onto the circle
		if (!is_flat){
			if(rsq>1.44*(mrsq)) {

				if(verbose) {
					msg(FUNC_NAME, TOO_FAR);
					printf("From (%g %g) to ", x, y);
				}
				rsq=sqrt(rsq)/mr;
				theta=arg(x,y);
				theta+=bias_angle(rsq);
				x=(mr)*cos(theta);
				y=(mr)*sin(theta);
				if(verbose) printf(" (%g, %g)\n",x,y);
			}
		}
		else{
			pbcremap(x, y);
		}
	}

	particle *p = put(x,y);
	return p;
}

/*  #################### #################### ##################### #################### ################## */

/*  #################### #################### Pacman Essentials #################### ################## */
/* Test whether the given particle is isolated,
 * meaning it does not have any of the 8 immediate neighbors.
 * \param[in] target the pointer to the particle to test on. */
bool dld_sim::test_isolation(particle *target){
	if(is_flat) fatal_error("test_isolation()", NO_PBC);

	//radius needed to search for 8 nearest neighbors
	bool isolate=true;
	double r = 4.*sfac;
	int ci, cj, cij, li, ui, lj, uj;

	// Find the extent of the boxes that are within the threshold distance
	// of the walker
	int_box(target->x,target->y,r,li,ui,lj,uj);

	// Loop over the boxes, looking for any immediate neighbors,
	// If found, then particle isn't isolated
	for(cj=lj;cj<=uj;cj++)  for(ci=li; ci<=ui;ci++){
		cij=ci+nh*cj;
		for(int qq=0; qq<co[cij];qq++) {
			particle &bb=ba[cij][qq];
			double rsq = dis_sq(bb.x, bb.y, target->x, target->y);

			if(rsq< 4. and rsq!=0){ isolate=false;}
		}
	}
	return isolate;
}

/* Returns the angle by which the pacman shape needs to
 * rotate so that its branch pts are on real axis.
 * \param[in] (x,y) the position of the random walker.
 * \param[in] target the pointer to the particle to be eliminated, */
double dld_sim::pacman_tilt(double x, double y, particle *target){
	double a = x-target->x;
	double b = y-target->y;
	double rot_angle = std::atan2(b,a) - pi/2.;
	return rot_angle;
}

/* Returns a random location on the pacman, where the random walker first reach the boundary.
 * reference for the equations and rationales used in this function, see pacman.tex.
 * \param[in] dsq the square of distance from walker to the center of nearest particle.
 * \param[in] sec_d the radius of pacman, i.e. distance from random walker to the edge of 2nd nearest neighbor.
 * \param[in] rot_angle the rotation angle to recover original pacman from the post mapping pacman.
 * \param[in] target the pointer to the nearest particle, we need the coordinate to shift to original pacman.
 * \param[in,out] eaten bool to indicate whether pacman has eaten a particle. */
Comp dld_sim::pacman(double dsq, double sec_d, double rot_angle, particle * target, bool &eaten){

	//inner arc radius
	double lilR=1.0;

	//distance between centers of outer and inner arc
	double dist = sqrt(dsq);

	//radius of the outer arc
	sec_d = sec_d - epsilon;

	//cos and sin of theta, theta is half of chord angle in the inner arc
	double cosT = (lilR*lilR + dsq - sec_d*sec_d)/(2*lilR*dist);

	//printf("cosTheta=%g\n", cosT);
	double sinT = sqrt(1.-cosT*cosT);
	double theta = acos(cosT);

	// half of chord angle in the outer arc
	double phi = acos((sec_d*sec_d + dsq - lilR*lilR)/(2.*dist*sec_d));

	// alpha * pi is the angle at branch point
	double al = (pi - theta -phi)/pi;

	// 1 over alpha, used in inverse map
	double alsp = 1.0 /al;

	//printf("alpha,theta,phi =%g, %g, %g\n", al, theta, phi);
	// conformal mapping coefficients
	Comp A = Comp(lilR*sinT*(cosT+1.), lilR*sinT*sinT);
	Comp B = Comp(-1.*lilR*sinT*(cosT+1.), 1*lilR*sinT*sinT);
	Comp C = Comp(1.+cosT, sinT);
	Comp D = Comp(1.+cosT, -1*sinT);

	// map pacman's center to upper plane, i.e. apply f^-1(w)
	// z plane: upper complex plane, w plane: pacman plane
	// pmw: pacman center in w plane, pmz: pacman center mapped to z plane
	/*****Suppose we've already
		1) moved inner arc center to origin
		2) rotate pacman so that the chord is parallel to real axis
		3) shift pacman so that the chord is on the real axis*/

	Comp pmw =Comp(0., dist-lilR*cosT);

	//printf("%g, %g\n",real(pmw), imag(pmw));
	// jj: dummy variable to calculate part of the inverse mapping
	Comp jj = pow((B-D*pmw)/(pmw*C-A),alsp);

	// Apply inverse mapping to find where pacman center pmw maps to in the upper complex plane
	Comp pmz = (jj+1.)/(1.-jj);

	// generate a random number between -pi/2 to pi/2, then take tangent of that.
	// then multiply by the y0, the imaginary part of pmz. (real part of pmz is 0)
	double probSample = imag(pmz) * tan((rnd()-0.5)*pi);
	if (probSample >= 1. or probSample <=-1.) eaten=true;

	// This is the location on real axis that random walker will first come into contact with
	Comp z = Comp(probSample, 0.);

	//Mapping from real axis to the adjusted pacman:
	Comp p1, p2;
	p1=pow((z-1.),al);
	p2=pow((z+1.),al);
	Comp w = (A*p1+B*p2)/(C*p1+D*p2);

	//Need to reverse all the shift and rotation to recover original pacman
	//Shift vertically, so that inner arc center is back at origin
	w += Comp(0., lilR*cosT);

	//imag(w) +=lilR*cosT;
	//Rotate by the rotate angle we found
	w = w*exp(Comp(0.,1.)*rot_angle);

	//Shift so that the center of the inner arc matches target particle location
	w += Comp(target->x, target->y);

	//real(w)+=target->x;
	//imag(w)+=target->y;
	//printf("%g, %g\n", real(w), imag(w));
	return w;
}


/* Find the second nearest neighbor, which determines the pacman big R.
 * Suitable for radial geometry simulations.
 * \param[in] (x,y) the position of the random walker.
 * \param[in,out] dsq the shortest distance from walker to cluster.
 * \param[in,out] sec_d the second shortest distance from walker to cluster, minus two times particles radius.
 * \param[in,out] rot_angle the angle to rotate to go from original mapping to the one drawn around walker.
 * \param[in] verbose if true debug messages are printed.*/
particle * dld_sim::min_distance(double x, double y, double &dsq,  double &sec_d, double &rot_angle, bool verbose){
	if(is_flat) fatal_error("min_distance()", NO_PBC);

	int ci,cj,cij,li,ui,lj,uj;
	particle *nearest=0;

	// Find the extent of the boxes that are within the threshold distance of the walker
	int_box(x,y,sqrt(dsq),li,ui,lj,uj);
	// Loop over the boxes, looking for any close particles
	for(cj=lj;cj<=uj;cj++) for(ci=li;ci<=ui;ci++) {
		cij=ci+nh*cj;
		for(int qq=0;qq<co[cij];qq++) {
			particle &bb=ba[cij][qq];
			double rsq=dis_sq(x,y,bb.x,bb.y);
			if(rsq<dsq){
				sec_d=dsq;
				dsq=rsq;
				nearest = &ba[cij][qq];
			}
			else {
				if(rsq<sec_d) sec_d=rsq;
			}

		}
	}

	// sec_d is really the distance from random walker to the edge of the 2nd nearest neighbor
	// but the particles in pacman has double the radius.
	double R = sqrt(dsq)+2*sfac;
	sec_d = sqrt(sec_d)-2*sfac;
	if (sec_d >= R) sec_d = R - 0.001;
	if(nearest!=NULL){
		rot_angle = pacman_tilt(x, y, nearest);
	}
	return nearest;

}


/*  #################### ##################### ################### #################### ################## */

/*  #################### #################### Sweep Essentials #################### ################## */
/* Check whether a box to be searched existed on the current simulation grid. */
bool dld_sim::exists(int s, int x, int y){
	bool exists = true;
	// When we don't impose period boundary condition, particle could be on the edge
	// with pbc, it's just gonna reemerge on the other side
	if((s%nh+x<0) or (s%nh+x>nh-1)) exists =false;
	if( (int)(s/nh)+y<0 or (int)(s/nh)+y >nv-1 ) exists =false;
	return exists;
}

/** Search one box for 1st and 2nd nearest neighbors. Different from the other min_distance() in that
 * in sweep method, there's a possibility that second nearest neighbor lie just outside the search threshold,
 * thus invisible to the search.
 * the safe approach for now is to only use pacman when both neighbors are found, use empty circle otherwise
 * \param[in] (x,y) the current position of the walker.
 * \param[in] s the block to search.
 * \param[in] stepx number of steps to take in x direction, can be negative
 * \param[in] stepy number of steps to take in y direction, can be negative
 * \param[in,out] dsq on input, the maximum threshold distance; on output, the minimum distance to the cluster, up to the threshold cutoff.
 * \param[in,out] sec_d on input, the large number to be updated; on output, the distance from walker to second nearest neighbor, minus 2 times particle radius.
 * \param[in] verbose if true debug messages are printed.*/
particle * dld_sim::min_distance(double x, double y, int s, int stepx, int stepy, double &dsq, double &sec_d, bool verbose){
	const char * FUNC_NAME = "min_distance()";
	particle *nearest=0;
	int sb=s+stepx+nh*stepy;
	if(verbose) printf("%s: walker at (%f %f) in box #%d, searching box #%d\n", FUNC_NAME, x, y, s, sb);
	for(int i=0; i<co[sb]; i++){
		particle &bb=ba[sb][i];
		if(verbose) printf("%s: Checking particle at (%f, %f)\n", FUNC_NAME, bb.x, bb.y);
		double rsq=dis_sq(x,y,bb.x,bb.y);
		if(rsq<dsq){
			sec_d=dsq;
			dsq=rsq;
			nearest = &ba[sb][i];

		}
		else {
			if(rsq<sec_d) {
				sec_d=rsq;
			}
		}

	}
	if(verbose) printf("%s: dsq %f, sec_d %f)\n", FUNC_NAME, dsq, sec_d);
	// If there isn't a closer neighbor in this box, the pointer returned is NULL
	return nearest;

}

/** Search the concentric square perimeter blocks around the block of the random walker to find a nearest neighbor.
 * Use with worklist generated by calculating bounding boxes.
 * \param[in,out] x current x coordinate of random walker, can be changed due to pbc.
 * \param[in] y current y coordinate of random walker.
 * \param[in,out] dsq the distance from random walker to the nearest neighbor squared.
 * \param[in,out] sec_d the distance from walker to the second nearest neighbor.
 * \param[in,out] rot_angle the amount to rotate compare to abstract pacman setup.
 * \param[in] verbose turns on debug messages.
 * \param[in] pbc turns on periodic boundary condition. */
particle * dld_sim::sweep(double& x, double y, double &dsq, double &sec_d, double &rot_angle, bool verbose){
	// Compute which block the particle is within
	const char * FUNC_NAME = "sweep(double, double, double, double, double, bool)";
	if(verbose) printf("%s\n%s: Random walker at %f, %f\n", DIV_LINE, FUNC_NAME, x, y);
	if (is_flat){
		if(x<-rad) x=fmod(x+rad, twidth)+redge;
		else if(x>=redge) x=fmod(x-redge, twidth)-rad;
		if(verbose) printf("%s: After pbc wrap-around, random walker at %f, %f\n", FUNC_NAME, x, y);
	}
	double ax = (-1)*rad;
	double ay = ax;
	int nx=int((x-ax)*xsp);if(nx<0) nx=0;if(nx>=nh) nx=nh-1;
	int ny=int((y-ay)*xsp);if(ny<0) ny=0;if(ny>=nv) ny=nv-1;
	int s=nx+nh*ny;

	// Compute which cell the particle is at within its own box
	// lower corner of particle home box
	double bx = nx*dx+ax;
	double by = ny*dx+ay;
	int cellx = int((x-bx)*dspcell); if(cellx<0) cellx=0; if(cellx>=nocell) cellx=nocell-1;
	int celly = int((y-by)*dspcell); if(celly<0) celly=0; if(celly>=nocell) celly=nocell-1;
	int c=cellx+nocell*celly;

	particle *min_pp=0;
	particle *nearest=0;
	std::vector<unsigned long int>::iterator it;

	int stepx=0, stepy=0;
	double bound_thr=1e6;
	double curr_thr=0;
	if(verbose) printf("%s: Particle is at row %d, column %d, dsq =%f, sec_d=%f\n", FUNC_NAME, ny, nx, dsq, sec_d);

	//move the iterator to the position after the worklist marker
	//end the iterator at the position before the next worklist marker
	for(it=wl.begin()+wlpos[c]+1;it!=wl.begin()+wlpos[c+1];++it){

		decoder((*it), stepx, stepy, curr_thr);
		if(verbose) printf("%s: checking worklist (%d, %d)\n", FUNC_NAME, stepx, stepy);
		curr_thr = curr_thr *cell_normalizer;

		//if the current threshold of a box if further away than the second nearest neighbor
		//no need to search that box at all
		if(curr_thr>bound_thr) {
			break;
		}

		if (!is_flat){
			if(exists(s,stepx,stepy)) {
				if(verbose) msg(FUNC_NAME, "Using min_distance()");
				min_pp=min_distance(x,y,s,stepx, stepy, dsq, sec_d, verbose);
				if(min_pp!=NULL) nearest=min_pp;
				bound_thr = sec_d;
			}
		}

		// Use periodic boundary condition, we decouple x,y steps
		// y is business as usual, x wraps around
		else {
			if(verbose) msg(FUNC_NAME, "Using PBC");
			if(exists(s,0,stepy)){
				if(verbose) printf("%s: Random walker at blc %d, taking step in x, y (%d, %d)\n",FUNC_NAME, s, stepx, stepy);

				// wrap x around as needed
				double ximage=x;
				if(!exists(s,stepx,0)){
					if(verbose) msg(FUNC_NAME, XSTEP_OUT_OF_BOUND);

					// the x grid position of the box to search after taking a step
					int tmpx=nx+stepx;

					// after wraping the box coordinate around, the walker's x coodinate also need to be mapped
					// by +2rad if box went out of bound on the left
					if(verbose) printf("%s: Stepping to x=%d\n", FUNC_NAME, tmpx);
					if (tmpx<=-1) {
						stepx = (nh+tmpx%nh)%nh-nx;
						if(verbose) printf("%s: Stepping off the grid on the left, wraps around by taking step %d in x\n ", FUNC_NAME, stepx);
						ximage = x+twidth;
					}
					else if (tmpx >= nh)  {
						stepx = tmpx%nh-nx;
						if(verbose) printf("%s: Stepping off the grid on the right, wraps around by taking step %d in x\n ", FUNC_NAME, stepx);
						ximage = x-twidth;
					}

				}
				if(verbose) printf("%s: Before finding min_dist: dsq = %f, sec_d = %f.\n", FUNC_NAME, dsq, sec_d);
				min_pp=min_distance(ximage,y,s,stepx, stepy, dsq, sec_d, verbose);
				if(verbose) printf("%s: After finding min_dist: dsq = %f, sec_d = %f.\n", FUNC_NAME, dsq, sec_d);
				if(min_pp!=NULL) {

					nearest=min_pp;
					rot_angle = pacman_tilt(ximage,y,nearest);
					if(verbose) printf("%s: rotation angle %f\n", FUNC_NAME, rot_angle);
				}

				// This bound_thr update might be redundant
				bound_thr = sec_d;

				// This search might be redundant too
				if (stepx==0  and nx==nh-1) {
					if(verbose) msg(FUNC_NAME, SPECIAL_SEARCH);
					if(verbose) printf("%s: Before finding min_dist: dsq = %f, sec_d = %f.\n",FUNC_NAME, dsq, sec_d);

					// If we search the first column, because their image live on the last column
					// need to transform x coord, otherwise nothing closer would be found
					min_pp = min_distance(x-twidth,y,s,-nh+1, stepy,dsq,sec_d, verbose);
					if(verbose) printf("%s: After finding min_dist: dsq = %f, sec_d = %f.\n",FUNC_NAME, dsq, sec_d);
					if(min_pp!=NULL){
						nearest=min_pp;
						rot_angle = pacman_tilt(ximage,y,nearest);
						if(verbose) printf("%s: Rotation angle %f\n", FUNC_NAME, rot_angle);
					}
					bound_thr = sec_d;
				}
			}
		}
	}
	if(verbose) printf("%s: Done. dsq=%f, sec_d =%f, rotation angle %f\n", FUNC_NAME, dsq, sec_d, rot_angle);

	// If second neighbor wasn't found, return no pacman stats, instead random walker would walk to the edge of the circle
	if (curr_thr<bound_thr) return NULL;
	else{
		double R = sqrt(dsq)+2*sfac;
		sec_d = sqrt(sec_d)-2*sfac;
		if (sec_d >= R) sec_d = R - 0.001;
		if(nearest!=NULL and !is_flat) rot_angle = pacman_tilt(x,y,nearest);
		if(verbose) printf("%s: After adjustment dsq=%f, sec_d =%f, rotation angle = %f\n",FUNC_NAME, dsq, sec_d, rot_angle);
		return nearest;
	}
}

/** Deposit values into wlpos, an array storing the starting position of each worklist,
  in the scenario with subgrids in the walker box. */
void dld_sim::wl_pos(){
	int pos;
	for(int i=1; i<=nocell*nocell; i++){
		pos = find(wl.begin(), wl.end(), (unsigned int) i) - wl.begin();
		wlpos[i-1] = pos;
	}
	wlpos[nocell*nocell] = wl.end()-wl.begin();
}

/** Decodes the worklist, which are 3 numbers encoded in a single unsigned int
 x step, y step, and the threshold at that step.
 * \param[in] encoded the integer with encoded worklist info.
 * \param[in,out] x step to take in x direction as a result of decoding.
 * \param[in,out] y step to take in y direction as a result of decoding.
 * \param[in,out] thr threshold of the current worklist as a result of decoding. */
void dld_sim::decoder(unsigned long int encoded, int &x_, int &y_, double &thr_){
	int x, y, sx, sy;
	unsigned long int mask = 65535;//First 2 bytes
	thr_ = encoded & mask;
	thr_ = double(thr_)*dx*dx;
	mask = 127; //  first byte from the right, minus the leftest bit in that byte
	y =encoded>>16 & mask;
	sy=encoded>>23 & 1;
	x =encoded>>24 & mask;
	sx =encoded>>31 & 1;
	if(sx==1) x*=(-1);
	if(sy==1) y*=(-1);
	x_=x; y_=y;
}


/*  #################### #################### ##################### #################### ################## */


/*  #################### #################### For the KDTree #################### ##################*/

/** Extras the information about simulation grid. */
Grid dld_sim::extract_grid(){
	//!!! MODIFY THIS TO ACCOMODATE NEW RECTANGULAR STRUCTURE
	Grid g = Grid(nh, co, rad, dx, ba);
	return g;
}


/** Performs a search for nearest neighbor using KDTree object. */
particle * dld_sim::treesearch(double x, double y, KDTree qt, double &dsq, double &sec_d, double &rot_angle){
	double original_secd = sec_d;
	particle * nearest = qt.closest(x, y, dsq, sec_d);
	if(nearest != NULL and sec_d<original_secd) {
		double R = sqrt(dsq)+2*sfac;
		sec_d = sqrt(sec_d)-2*sfac;
		if (sec_d >= R) sec_d = R - 0.001;
		if(nearest!=0) rot_angle = pacman_tilt(x, y, nearest);
		return nearest;
	}
	else return NULL;
}
/*  #################### #################### ##################### #################### ################## */
// End of Class DLD
