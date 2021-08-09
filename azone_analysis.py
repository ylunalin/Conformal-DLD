import os
import re
import sys
import time
import getopt
import fnmatch
import argparse
import numpy as np
from scipy import optimize
from os import listdir
from os.path import isfile, isdir, join
from matplotlib.pyplot import cm
from matplotlib import pyplot as plt

# Using argparse module to parse the commandline argument
def main_parse():
	"""
	Create parser and parse the commandline arguments.

	Parameters
	----------
	commandline arguments, simulation directory is required.

	Returns
	-------
	parser : Argument parser

	args : Namespace
		The arguments list resulting from parser parsing commandline arguments.
	"""
	parser = argparse.ArgumentParser()

	parser.add_argument("sim_dir", help="the path to simulation directory")
	parser.add_argument("-a", "--all",action="store_true", default=False, help="process all the active zone files")
	parser.add_argument("-c", "--cluster",action="store_true", default=False, help="the data direcotry is from cluster simulations")
	parser.add_argument("-f", "--fname_prefix", metavar = 'FN', default="azone_N", help="the filename prefix")
	parser.add_argument("-s", "--single", metavar = 'N', default=[], nargs='+', help="process specific active zone file(s)")
	parser.add_argument("-r", "--aspect_ratio", metavar = 'R', default=1, help="aspect ratio")
	parser.add_argument("-v", "--verbosity", action="count", default=0, help="levels of output verbosity")
	parser.add_argument("-x", "--show_plot", action="store_true", default=False, help="show the plots")
	try:
		args = parser.parse_args()
	except SystemExit:
		print("Arguments can't be parsed. Try agian.");
		return parser, None

	# if no single files are specified, process all existing files
	if len(args.single)==0:
		args.all = True

	# if verbose, print a summary
	if args.verbosity>0:
		print("Summary:")
		print("        simulation directory  ==>",args.sim_dir)
		print("        prefix of filename    ==>", args.fname_prefix)
		if args.all:
			print("        processing all frames")
		else:
			for s in args.single:
				print("        processing frame      ==> \"{}/{}{}\"".format(args.sim_dir, args.fname_prefix, s))
	return parser, args;

def stat(args):
	"""
	Get particle, mean, width, and skewness from data files specified by args.

	Parameters
	----------
	args : Namespace
		The arguments list resulting from parser parsing commandline arguments.

	Returns
	-------
	name : ndarray
		The active zone collection points. Number of particles remaining in clusters, number of particles dissolved in flat.

	mean, width, skew : ndarray
		The mean, width, skewness of active zone at the collection points.
	"""
	# parse the arguments
	sim_dir = args.sim_dir
	fname = args.fname_prefix
	fnums = args.single
	proc_all = args.all
	if len(fnums)==0: proc_all = True;
	show = args.show_plot
	# we get all the files to process
	files = []
	if proc_all:
		files = [join(sim_dir,f) for f in listdir(sim_dir) if isfile(join(sim_dir, f)) and fnmatch.fnmatch(f,fname+"*")]
	# or get the specified files to process
	else:
		fullname = join(sim_dir,fname)
		files = [fullname+num for num in fnums if isfile(fullname+num)]

	# get a sample file, since the number of realizations are the same,
	# we can create subsequent data array accordingly
	ex = np.array([])
	ex = np.fromfile(files[0], dtype = np.float64)
	azone_data = np.empty([len(files), len(ex)])
	name = np.empty(len(files))
	for i, f in enumerate(files):
		azone_data[i] = np.fromfile(f, dtype=np.float64)
		fid = re.findall('\d+',f)[-1]
		name[i] = fid

	mean = np.mean(azone_data, axis=1)
	var = np.var(azone_data, axis=1)
	width = np.sqrt(var)
	skew = skewness(azone_data, mean, var)
	if(show):
		plt.loglog(name,mean,'bo',label='mean')
		plt.loglog(name,width,'gv',label='width')
		plt.loglog(name,skew, 'rd',label='skewness')
		plt.xlabel("time (number of particle remaining (cluster) or dissolved (flat))")
		plt.legend(loc='best')
		plt.show()

	name = np.reshape(name, (-1,1))
	mean = np.reshape(mean, (-1,1))
	width = np.reshape(width, (-1,1))
	skew = np.reshape(skew, (-1,1))
	output = np.concatenate((name, mean), axis=1)
	output = np.concatenate((output, width), axis=1)
	output = np.concatenate((output, skew), axis=1)

	fn = "aggregate_"+fname+".dat"
	np.savetxt(join(sim_dir, fn), output, header="N, mean, std, skew")
	return name, mean,width,skew

def skewness(data,mean,var):
	"""
	Calculate skewness of a data set.

	Parameters
	----------
	data, mean, var : ndarray
		The data set, a 2D for all the data at each active zone point. The mean and variance (1D).

	Returns
	-------
	skewness : ndarray
		The skewness for each active zone point. Same dimension as mean and variance.
	"""
	m,n = data.shape
	mmt3 = np.mean((data-mean.reshape(m,1))**3, axis=1)
	return mmt3 / var**(3/2)


# Using the commandline args, plot out active zone data files as histograms
def hist(args):
	# parse the arguments
	sim_dir = args.sim_dir
	fname = args.fname_prefix
	fnums = args.single
	proc_all = args.all
	if len(fnums)==0: proc_all = True;
	show = args.show_plot
	# we get all the files to process
	files = []
	if proc_all:
		files = [join(sim_dir,f) for f in listdir(sim_dir) if isfile(join(sim_dir, f)) and fnmatch.fnmatch(f,fname+"*")]
	# or get the specified files to process
	else:
		fullname = join(sim_dir,fname)
		files = [fullname+num for num in fnums if isfile(fullname+num)]

	# get a sample file, since the number of realizations are the same,
	# we can create subsequent data array accordingly
	ex = np.array([])
	ex = np.fromfile(files[0], dtype = np.float64)
	azone_data = np.empty([len(files), len(ex)])
	name = np.empty(len(files))
	for i, f in enumerate(files):
		azone_data[i] = np.fromfile(f, dtype=np.float64)
		fid = re.findall('\d+',f)[2]
		plt.hist(azone_data[i], bins=20, normed=True, histtype='stepfilled', alpha=0.5, label="N="+fid)
	plt.legend(loc='upper right', fontsize='xx-small', ncol=5)
	if(show):plt.show()

def load(args):
	"""
	Load the dataset into memory, directory, filename prefix, specified in args.

	Parameters
	----------
	args : Namespace
		The arguments list resulting from parser parsing commandline arguments.

	Returns
	-------
	name : ndarray
		The index of the particle where azone data is recorded.
	azone_data : ndarray
		The active zone dataset.

	"""
	# parse the arguments
	sim_dir = args.sim_dir
	fname = args.fname_prefix
	fnums = args.single
	proc_all = args.all
	if len(fnums)==0: proc_all = True;

	# we get all the files to process
	files = []
	if proc_all:
		files = [join(sim_dir,f) for f in listdir(sim_dir) if isfile(join(sim_dir, f)) and fnmatch.fnmatch(f,fname+"*")]
	# or get the specified files to process
	else:
		fullname = join(sim_dir,fname)
		files = [fullname+num for num in fnums if isfile(fullname+num)]

	# get a sample file, since the number of realizations are the same,
	# we can create subsequent data array accordingly
	ex = np.array([])
	ex = np.fromfile(files[0], dtype = np.float64)
	azone_data = np.empty([len(files), len(ex)])
	name = np.empty(len(files))
	for i, f in enumerate(files):
		azone_data[i] = np.fromfile(f, dtype=np.float64)
		fid = re.findall('\d+',f)[-1]
		name[i] = fid
	return name, azone_data;

def grid_alignment(args):
	"""
	Load the dataset into memory, directory, filename prefix, specified in args.
	Then compute the grid alignment parameter <cos(4theta)>

	Parameters
	----------
	args : Namespace
		The arguments list resulting from parser parsing commandline arguments.

	Returns
	-------
	name : ndarray
		The index of the particle where azone data is recorded.
	rcos4t : ndarray
		The alignment parameter

	"""
	args.fname_prefix = "azone_X"
	x, X = load(args)

	args.fname_prefix = "azone_Y"
	x, Y = load(args)
	# we compute <cos(4theta)>
	data_shape = X.shape
	rcos4t = np.zeros(data_shape[0])
	theta = np.empty(data_shape[1])
	for i in range(data_shape[0]):
		theta = np.arctan2(Y[i], X[i])
		rcos4t[i] = np.mean(np.cos(4*theta))
	return x, rcos4t

def ellipse_azone_width(args):
	aspect_ratio = np.float(args.aspect_ratio)
	args.fname_prefix = "azone_O"
	x, R = load(args)
	sorted_ind = np.argsort(x)
	R = R[sorted_ind]
	meanR = np.mean(R, axis=1)
	majorAxis = np.array([solve_ellipse_axes(r, aspect_ratio) for r in meanR])
	minorAxis = majorAxis/aspect_ratio

	args.fname_prefix = "azone_X"
	x, X = load(args)
	sorted_ind = np.argsort(x)
	X = X[sorted_ind]

	args.fname_prefix = "azone_Y"
	x, Y = load(args)
	sorted_ind = np.argsort(x)
	Y = Y[sorted_ind]
	x = x[sorted_ind]

	rbar = np.zeros_like(R)
	for i in range(len(x)):
		rbar[i] = ellipse_radius(X[i], Y[i], majorAxis[i], minorAxis[i])

	var = np.mean((R-rbar)**2, axis=1)

	x = x.reshape(-1,1)
	std = np.sqrt(var).reshape(-1,1)
	output = np.concatenate( (x, std), axis=1)
	output = np.concatenate( (output, majorAxis.reshape(-1,1)), axis=1)
	output = np.concatenate( (output, minorAxis.reshape(-1,1)), axis=1)
	fn = "ellipse_width.dat"
	np.savetxt(join(args.sim_dir, fn), output, header="N, std, mean_major, mean_minor")
	return x, std, majorAxis, minorAxis

def ellipse_radius(x, y, a, b):
	r2 = x*x + y*y
	costheta2 = x*x / r2
	sintheta2 = y*y / r2
	radius = a*b / np.sqrt(a*a*sintheta2 + b*b*costheta2)
    # TODO this is tricky, when (x,y)=(0,0), we can't calculate angle, so no ellipse radius
	radius[r2==0] = (a+b)*0.5
	return radius


# A few functions to process ellipse clusters
def solve_ellipse_axes(rbar, aspect_ratio):
	x0 = rbar
	args = (aspect_ratio, rbar, )
	guess = optimize.newton(ellipse_circum_newton, x0, args=args)
	return guess

def ellipse_circum_newton(a,aspect_ratio, rbar):
	# Note we have divided the circumference formula by 2PI
	b = a/aspect_ratio
	three_h = (a-b)**2 / (a+b)**2 * 3
	return rbar - 0.5* (a+b) * (1 + three_h / (10+np.sqrt(4-three_h)))

def ellipse_circumference(a,aspect_ratio):
	# Note we have divided the circumference formula by 2PI
	b = a/aspect_ratio
	three_h = (a-b)**2 / (a+b)**2 * 3
	return  np.pi * (a+b) * (1 + three_h / (10+np.sqrt(4-three_h)))

# give log x and log data, perform least square fit to get a coefficient for the line
def lst_sqr_fit(x,data):
	"""
	Given log x and log of data (y), perform linear least square fit to get a coefficient.

	Parameters
	----------
	x : ndarray
		Log of independent variable.
	data : ndarray
		Log of dependent variable.

	Returns
	-------
	x : ndarray
		Coefficents, i.e. constant term and the slope of the linear fit.
	"""
	a11 = np.sum(x**2)
	a12 = np.sum(x)
	a22 = len(x)
	b1	= np.sum(x*data)
	b2	= np.sum(data)
	A	= np.array([[a11,a12],[a12,a22]])
	b	= np.array([b1,b2])
	x	= np.linalg.solve(A,b)
	return x

def proc_flat(args):
	"""
	Get particle, mean, width, and skewness from data subdirectory in dir specified by args.

	Parameters
	----------
	args : Namespace
		The arguments list resulting from parser parsing commandline arguments.

	Returns
	-------
	true_angles : ndarray
		The angles of flat substrates.

	avg_ws, avg_ss : ndarray
		The average saturation width and skewness for each angle.

	xs : ndarray
		The active zone collection points. Number of particles remaining in clusters, number of particles dissolved in flat.

	means, widths, skews : ndarray
		The mean, width, skewness of active zone at the collection points.
	"""
	n			= 0
	local_time	= time.time()
	sim_dir		= args.sim_dir
	sub_dirs	= listdir(sim_dir)
	for s in sub_dirs:
		stmp = join(sim_dir, s)
		# check that it's a diretory, there's the size_angle file too
		if isdir(stmp):
			n +=1

	xs		= np.empty(n, dtype=np.ndarray)
	means	= np.copy(xs)
	widths	= np.copy(xs)
	skews	= np.copy(xs)
	angles	= np.empty(n)
	sub_args= []
	n		= 0

	for s in sub_dirs:
		stmp = join(sim_dir, s)
		# check that it's a diretory, there's the size_angle file too
		if not isdir(stmp):
			continue
		# get the args that contains a subdirectory
		tmp_args = parser.parse_args(["-a", "-f", args.fname_prefix, stmp])
		angles[n] =re.findall("\d+", stmp)[-1]
		xs[n], means[n], widths[n], skews[n] = stat(tmp_args)

		sub_args.append(tmp_args)
		n += 1

	avg_ws = np.empty(n)
	avg_ss = np.empty(n)
	# get the averages
	for i in range(n):
		avg_ws[i] = np.mean(widths[i])
		avg_ss[i] = np.mean(skews[i])

	size_angle	= np.loadtxt(join(sim_dir, "size_angle"), skiprows=1, unpack=True)
	true_angles = np.copy(angles)
	true_widths = np.copy(angles)
	for i in range(n):
		true_widths[i] = size_angle[2][size_angle[1]==angles[i]][0]
		true_angles[i] = size_angle[-1][size_angle[1]==angles[i]][0]

	if args.show_plot:
		plt.plot(true_angles, avg_ws/np.log(true_widths), 'ro')
		plt.xlabel("angles (degs)")
		plt.ylabel("std. dev. in width")
		plt.title("Normalized saturation width of flat subtrate via active zone")
		plt.show()
		plt.plot(true_angles, avg_ss/np.log(true_widths), 'ro')
		plt.xlabel("angles (degs)")
		plt.ylabel("skewness")
		plt.title("Normalized saturation skewness of flat subtrate via active zone")
		plt.show()

	# reorder the array in acsending order
	angle_order	= np.argsort(true_angles)
	true_angles	= true_angles[angle_order]
	avg_ws		= avg_ws[angle_order]
	avg_ss		= avg_ss[angle_order]
	xs			= xs[angle_order]
	means		= means[angle_order]
	widths		= widths[angle_order]
	skews		= skews[angle_order]
	saveout		= np.concatenate((true_angles, avg_ws, avg_ss), axis=0)
	saveout.tofile(join(sim_dir, "avg_stat.bin{}".format(saveout.shape[0])))
	# finish timing, report
	local_time	= time.time() - local_time;
	print("Time elapsed {:10.2f} seconds".format(local_time));

	return true_angles, avg_ws, avg_ss, xs, means, widths, skews

if __name__ == "__main__":
	parser, args = main_parse()
	# if cluster simulation (circular, adld shapes)
	# or if a subdirectory of flat is specified
	if args is not None and args.cluster :
		dazone_data = load(args)
		x, mean, width, skew = stat(args)


	# if not, then it's the entire flat geometry simulation
	# in that case there are subdirectories for each angle
	elif args is not None:
		true_angles, avg_ws, avg_ss, xs, means, widths, skews = proc_flat(args)
