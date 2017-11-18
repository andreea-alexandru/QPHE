#!/usr/bin/env python3
import numpy

def RandOrthMat(n):
	tol=10**(-6)
	M = numpy.zeros((n,n)) # prealloc
	# gram-schmidt on random column vectors
	vi = numpy.random.normal(0,1,n)
	# the n-dimensional normal distribution has spherical symmetry, which implies
	# that after normalization the drawn vectors would be uniformly distributed on the
	# n-dimensional unit sphere.

	M[:][0] = vi/numpy.linalg.norm(vi);
	for i in range(1,n):
		nrm = 0
		while nrm<tol:
			vi = numpy.random.normal(0,1,n)
			x = numpy.dot(M[:][0:i],vi)
			vi = vi - numpy.dot(x,M[:][0:i])
			nrm = numpy.linalg.norm(vi);
		M[:][i] = vi / nrm;
	return M


def write_matrices(n,m):
	rangeMINeig = 1;
	rangeMAXeig = 6;
	rangeMIN = -5;
	rangeMAX = 5;
	fileQ = 'Q'+str(n)+'_'+str(m)+'.txt'
	fileA = 'A'+str(n)+'_'+str(m)+'.txt'
	fileb = 'b'+str(n)+'_'+str(m)+'.txt'
	filec = 'c'+str(n)+'_'+str(m)+'.txt'
	fileparam = 'param'+str(n)+'_'+str(m)+'.txt'

	r = numpy.random.rand(n)*(rangeMAXeig - rangeMINeig) + rangeMINeig;
	D = numpy.zeros((n,n))
	numpy.fill_diagonal(D, r)
	V = RandOrthMat(n);
	Q = numpy.dot(numpy.dot(V,D),numpy.transpose(V))
	numpy.savetxt(fileQ, Q, fmt='%7.4f', delimiter=",")
	r = numpy.random.rand(n*m)*(rangeMAX - rangeMIN) + rangeMIN;
	A = r.reshape(m,n)
	numpy.savetxt(fileA, A, fmt='%7.4f', delimiter=",")
	b = numpy.random.rand(m)*(rangeMAX - rangeMIN) + rangeMIN;
	numpy.savetxt(fileb, b, fmt='%7.4f')
	c = numpy.random.rand(n)*(rangeMAX - rangeMIN) + rangeMIN;
	numpy.savetxt(filec, c, fmt='%7.4f')

	Qinv = numpy.linalg.inv(Q);
	eps = 0.001;
	AQA = numpy.dot(numpy.dot(A,Qinv),numpy.transpose(A))
	spec,w = numpy.linalg.eig(AQA);
	maxEig = numpy.real(numpy.amax(spec))
	minEig = numpy.real(numpy.amin(spec))
	if m<=n:
	    k_number = maxEig/minEig;
	    K = numpy.floor(numpy.log(1/eps)*k_number/2);
	else:
		K = numpy.floor(numpy.log(1/eps)*maxEig);
	eta = 1/maxEig;
	print(K)
	numpy.savetxt(fileparam, (K, eta), fmt='%-7.6f')


def main():
	# for n in range(75,81,10):
	# 	for m in range(45,46,8):
	# 		print(n,m)
	# 		write_matrices(n,m)
	n = 5
	m = 5
	write_matrices(n,m)

main()
