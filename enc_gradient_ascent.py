#!/usr/bin/env python3

from gmpy2 import mpz
import numpy
import time
import random
import paillier

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False

DEFAULT_MSGSIZE = 32
DEFAULT_PRECISION = int(DEFAULT_MSGSIZE/2) # of fractional bits
DEFAULT_KEYSIZE = 512

def encrypt_vector(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.encrypt(y) for y in x]
	else: return [pubkey.encrypt(y,coins.pop()) for y in x]

def decrypt_vector(privkey, x):
    return numpy.array([privkey.decrypt(i) for i in x])

def fixed_point(scalar,prec=DEFAULT_PRECISION):
	return mpz(scalar*(2**prec))

def fixed_point_vector(vec,prec=DEFAULT_PRECISION):
	if numpy.size(vec)>1:
		return [fixed_point(x) for x in vec]
	else:
		return fixed_point(vec)

def fixed_point_matrix(mat,prec=DEFAULT_PRECISION):
	return [fixed_point_vector(x) for x in mat]

def retrieve_fixed_point(scalar,prec=DEFAULT_PRECISION):
	# return mpz(scalar/(2**prec))
	return gmpy2.f_div_2exp(mpz(scalar),prec)

def retrieve_fixed_point_vector(vec,prec=DEFAULT_PRECISION):
	return [retrieve_fixed_point(x) for x in vec]

def sum_encrypted_vectors(x, y):
	return [x[i] + y[i] for i in range(numpy.size(x))]

def diff_encrypted_vectors(x, y):
	return [x[i] - y[i] for i in range(len(x))] 

def mul_sc_encrypted_vectors(x, y): # x is encrypted, y is plaintext
    return [y[i]*x[i] for i in range(len(x))]    

def dot_sc_encrypted_vectors(x, y): # x is encrypted, y is plaintext
    return sum(mul_sc_encrypted_vectors(x,y))

def dot_m_encrypted_vectors(x, A):
    return [dot_sc_encrypted_vectors(x,vec) for vec in A]

def grad_sol_strong(A,b,c,m,Q,K,eta):
	At = A.transpose()
	minvQ = -numpy.linalg.inv(Q)			### - Q^{-1}
	AinvQ = numpy.dot(A,minvQ)				### - AQ^{-1}
	coeff_mu = fixed_point_matrix(numpy.identity(m) + eta*numpy.dot(AinvQ,At)) ### I-\eta AQ^{-1}A'
	coeff_c = fixed_point_matrix(eta*AinvQ) 					### -\etaAQ^{-1}
	coeff_x = fixed_point_matrix(numpy.dot(minvQ,At))
	minvQ = fixed_point_matrix(minvQ)
	meta = fixed_point(-eta)

	keypair = paillier.generate_paillier_keypair(n_length=DEFAULT_KEYSIZE)
	pubkey, privkey = keypair
	enc_b = encrypt_vector(pubkey,b)
	enc_c = encrypt_vector(pubkey,c)

	# mu = numpy.random.randint(-1,1, size=m)
	mu = numpy.zeros(m).astype(int)
	for k in range (0,K):
		emu = encrypt_vector(pubkey,mu)
		# mu_bar = numpy.dot(coeff_mu,mu)+numpy.dot(coeff_c,c)+[x*meta for x in b]	### mu_bar * 2**(2f)
		emu_bar = sum_encrypted_vectors(numpy.dot(coeff_mu,emu),numpy.dot(coeff_c,enc_c))
		emu_bar = sum_encrypted_vectors(emu_bar,[x*meta for x in enc_b])
		mu_bar = decrypt_vector(privkey,emu_bar)
		mu_bar = retrieve_fixed_point_vector([int(x) for x in mu_bar])		### mu_bar * 2**f
		# # mu_bar = [mpz(x/(2**DEFAULT_PRECISION))for x in mu_bar]
		# print(mu_bar)
		mu = numpy.maximum([0 for x in range(0,m)],mu_bar)
		# print(mu)
	x = numpy.dot(coeff_x,mu)+numpy.dot(minvQ,c)
	x = x/2**(2*DEFAULT_PRECISION)
	return x


def main():
	n = 20;
	m = 13;
	fileA = "Data/A"+str(n)+"_"+str(m)+".txt"
	fileQ = "Data/Q"+str(n)+"_"+str(m)+".txt"
	fileb = "Data/b"+str(n)+"_"+str(m)+".txt"
	filec = "Data/c"+str(n)+"_"+str(m)+".txt"
	fileparam = "Data/param"+str(n)+"_"+str(m)+".txt"
	A = numpy.loadtxt(fileA, delimiter=',')
	Q = numpy.loadtxt(fileQ, delimiter=',')
	param = numpy.loadtxt(fileparam, delimiter='\n')
	# K = int(param[0])
	K = 10
	eta = param[1]
	b_A = numpy.loadtxt(fileb, delimiter='\n')
	c_A = numpy.loadtxt(filec, delimiter='\n')
	b_A = fixed_point_vector(b_A)
	c_A = fixed_point_vector(c_A)

	x = grad_sol_strong(A,b_A,c_A,m,Q,K,eta)
	# print("%.4f" % x)
	print(["%.8f"% i for i in x])

main()