import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

def ndft(tn, pn):
	# sort first
	l		= tn.shape[0]
	order	= np.argsort(tn)
	T		= tn[order]
	P		= pn[order]
	A		= np.zeros((l,l))

	# put together the matrix, i.e. the row vectors
	dk		= np.pi/T[-1]
	omegas = np.zeros_like(T,dtype=np.complex)
	for m in range(l):
		exponent = -1j*m*dk*T
		row_vec = np.exp(exponent)
		omegas[m] = np.dot(row_vec.T, P)

	return omegas


if __name__=="__main__":
	tn, wn, sn = np.loadtxt("nflat_w100_nr150k/avg_stat", unpack=True)
	omegas100 = ndft(tn,wn)
	tn, wn, sn = np.loadtxt("nflat_w200_nr150k/avg_stat", unpack=True)
	omegas200 = ndft(tn,wn)

	color = iter(cm.Accent(np.linspace(0,1,6)))
	c = next(color)
	plt.title("Non-uniform DFT of flat width data, w=100, 200")
	plt.plot(range(len(omegas100)), np.real(omegas100), 'o-', c=c,label="w=100, real component")
	c = next(color)
	plt.plot(range(len(omegas100)), np.imag(omegas100), 'o-', c=c, label="w=100, imag component")
	c = next(color)
	plt.plot(range(len(omegas100)), np.abs(omegas100), 'o-', c=c, label="w=100, magnitude")
	c = next(color)
	plt.plot(range(len(omegas200)), np.real(omegas200), 'o-', c=c, label="w=200, real component")
	c = next(color)
	plt.plot(range(len(omegas200)), np.imag(omegas200), 'o-', c=c, label="w=200, imag component")
	c = next(color)
	plt.plot(range(len(omegas200)), np.abs(omegas200), 'o-', c=c, label="w=200, magnitude")
	plt.legend(loc='best')
	plt.show()
