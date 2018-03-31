#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
import numpy as np
import time
import testDGK
import random


DEFAULT_KEYSIZE = 512
DEFAULT_MSGSIZE = 32 # The message size of DGK has to be greater than 2*log2(DEFAULT_MSGSIZE), check u in DGK_pubkey
DEFAULT_PRECISION = int(DEFAULT_MSGSIZE/2) # of fractional bits
DEFAULT_SECURITYSIZE = 80
DEFAULT_STATISTICAL = 100
NETWORK_DELAY = 0 #ms

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False

seed = 42

def encrypt_vector(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.encrypt(y) for y in x]
	else: return [pubkey.encrypt(y,coins.pop()) for y in x]

def encrypt_matrix(pubkey, x, coins=None):
	if (coins==None):
		return [[pubkey.encrypt(y) for y in z] for z in x]
	else: return [[pubkey.encrypt(y,coins.pop()) for y in z] for z in x]

def decrypt_vector(privkey, x):
    return [privkey.decrypt(i) for i in x]

def sum_encrypted_vectors(x, y):
	return [x[i] + y[i] for i in range(np.size(x))]

def diff_encrypted_vectors(x, y):
	return [x[i] - y[i] for i in range(len(x))] 

def mul_sc_encrypted_vectors(x, y): # x is encrypted, y is plaintext
    return [y[i]*x[i] for i in range(len(x))]    

def dot_sc_encrypted_vectors(x, y): # x is encrypted, y is plaintext
    return sum(mul_sc_encrypted_vectors(x,y))

def dot_m_encrypted_vectors(x, A):
    return [dot_sc_encrypted_vectors(x,vec) for vec in A]

def encrypt_vector_DGK(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.raw_encrypt(y) for y in x]
	else: return [pubkey.raw_encrypt(y,coins.pop()) for y in x]

def decrypt_vector_DGK(privkey, x):
    return np.array([privkey.raw_decrypt0(i) for i in x])

####### We take the convention that a number x < N/3 is positive, and that a number x > 2N/3 is negative. 
####### The range N/3 < x < 2N/3 allows for overflow detection.

def fp(scalar,prec=DEFAULT_PRECISION):
	if prec < 0:
		return gmpy2.t_div_2exp(mpz(scalar),-prec)
	else: return mpz(gmpy2.mul(scalar,2**prec))

def fp_vector(vec,prec=DEFAULT_PRECISION):
	if np.size(vec)>1:
		return [fp(x,prec) for x in vec]
	else:
		return fp(vec,prec)

def fp_matrix(mat,prec=DEFAULT_PRECISION):
	return [fp_vector(x,prec) for x in mat]

def retrieve_fp(scalar,prec=DEFAULT_PRECISION):
	return scalar/(2**prec)
	# return gmpy2.div(scalar,2**prec)

def retrieve_fp_vector(vec,prec=DEFAULT_PRECISION):
	return [retrieve_fp(x,prec) for x in vec]

def retrieve_fp_matrix(mat,prec=DEFAULT_PRECISION):
	return [retrieve_fp_vector(x,prec) for x in mat]


class Agents:
    def __init__(self, pubkey,fileb,filec):
    	self.pubkey = pubkey
    	b_A = np.loadtxt(fileb, delimiter='\n')
    	c_A = np.loadtxt(filec, delimiter='\n')
    	self.enc_b_A = encrypt_vector(pubkey,fp_vector(b_A))
    	self.enc_c_A = encrypt_vector(pubkey,fp_vector(c_A))
    	self.m = np.size(b_A)
    	self.n = np.size(c_A)

    def send_data(self):
    	return self.enc_b_A, self.enc_c_A, self.m, self.n

class Cloud:
	def __init__(self, pubkey, DGK_pubkey, fileA, fileQ, fileparam):
		self.pubkey = pubkey
		self.DGK_pubkey = DGK_pubkey
		self.N = pubkey.n
		self.N2 = pubkey.nsquare
		self.N_len = (self.N).bit_length()
		self.l = DEFAULT_MSGSIZE
		self.sigma = DEFAULT_SECURITYSIZE		
		A = np.loadtxt(fileA, delimiter=',')
		Q = np.loadtxt(fileQ, delimiter=',')
		self.A = A
		m = np.size(A,0)
		self.m = m
		At = A.transpose()
		n = np.size(Q,0)
		self.Q = Q
		invQ = np.linalg.inv(Q)			### Q^{-1}
		AinvQ = np.dot(A,invQ)			###  AQ^{-1}
		AinvQA = np.dot(AinvQ,At)		### AQ^{-1}A'
		eigs = np.linalg.eigvals(AinvQA)	
		eta = 1/np.real(max(eigs))
		self.delta_A = [0]*m
		# param = np.loadtxt(fileparam, delimiter='\n')
		# self.K = param[0]
		# self.K = int(self.K)
		self.K = 30

		coeff_mu = fp_matrix(np.identity(m) - eta*AinvQA) ### I-\eta AQ^{-1}A'
		self.coeff_mu = coeff_mu
		coeff_c = fp_matrix(-eta*AinvQ) 					### -\etaAQ^{-1}
		self.coeff_c = coeff_c
		coeff_muK = fp_matrix(np.dot(-invQ,At))			### -Q^{-1}A'
		self.coeff_muK = coeff_muK
		coeff_cK = fp_matrix(-invQ)
		self.coeff_cK = coeff_cK
		etabar = fp(-eta)	
		self.etabar = etabar	
		self.gen_rands()

	def gen_rands(self):
		lf = DEFAULT_PRECISION
		m = self.m
		l = self.l
		sigma = self.sigma
		lambd = DEFAULT_STATISTICAL
		K = self.K
		random_state = gmpy2.random_state(seed)
		mu = np.zeros(m).astype(int)
		# mu = fp_vector([gmpy2.mpz_urandomb(random_state,self.l-DEFAULT_PRECISION-1) for i in range(0,m)])
		self.mu = encrypt_vector(self.pubkey, mu)
		# Noise for blinding mu in the update step
		rn = [[[gmpy2.mpz_urandomb(random_state,l+lambd),gmpy2.mpz_urandomb(random_state,l + lambd)] for i in range(0,m)] for k in range(0,K)]
		self.obfuscations = rn
		# Noise for comparison
		rn = [[gmpy2.mpz_urandomb(random_state,l+lambd) for i in range(0,m)] for k in range(0,K)]
		self.rn = rn
		random_state = gmpy2.random_state(seed)
		coinsP = [gmpy2.mpz_urandomb(random_state,self.N_len-1) for i in range(0,4*m*K)]
		coinsP = [gmpy2.powmod(x, self.N, self.N2) for x in coinsP]
		coinsDGK = [gmpy2.mpz_urandomb(random_state,2*sigma) for i in range(0,2*(l+1)*m*K)]
		coinsDGK = [gmpy2.powmod(self.DGK_pubkey.h, x, self.DGK_pubkey.n) for x in coinsDGK]
		self.coinsDGK = coinsDGK
		# Noise for truncation
		rn = [gmpy2.mpz_urandomb(random_state,lambd+l+lf) for i in range(0,m*K)]
		self.fixedNoise = encrypt_vector(self.pubkey, rn, coinsP[-m*K:])
		er = [-fp(x,-lf) for x in rn]
		er = encrypt_vector(self.pubkey,er,coinsP[-2*m*K:-m*K])
		self.er = er
		coinsP = coinsP[:-2*m*K]
		self.coinsP = coinsP


	def compute_grad(self,b_A,c_A):
		mu_bar = sum_encrypted_vectors(np.dot(self.coeff_mu,self.mu),np.dot(self.coeff_c,c_A))
		mu_bar = sum_encrypted_vectors(mu_bar,[x*self.etabar for x in b_A])
		self.mu_bar = mu_bar ### \mu_bar*2^{2*lf}

	def temporary_prec_mu(self):
		m = self.m
		pubkey = self.pubkey
		r = [self.fixedNoise.pop() for i in range(0,m)]
		temp_mu = sum_encrypted_vectors(self.mu_bar,r)		### mu_bar*2**(2*lf)+r
		return temp_mu
	
	def compute_primal_optimum(self,c_A):
		x = np.dot(self.coeff_muK,self.mu)
		x = sum_encrypted_vectors(x,np.dot(self.coeff_cK,c_A))
		return x

	def randomize(self):
		m = self.m
		a = [0]*m
		b = [0]*m
		for i in range(0,m):
			a[i],b[i] = np.random.permutation([self.pubkey.encrypt(0),self.mu_bar[i]])
		self.a = a
		self.b = b
		### HAVE TO BE NUMBERS OF l BITS
		return a,b

	def init_comparison_cloud(self):
		m = self.m
		l = self.l
		pubkey = self.pubkey
		r = self.r
		a,b = self.randomize()
		z = diff_encrypted_vectors(b,a)
		z = sum_encrypted_vectors(z,encrypt_vector(pubkey,r,self.coinsP[-m:]))
		z = sum_encrypted_vectors(z,encrypt_vector(pubkey,[2**l]*m,self.coinsP[-2*m:-m]))
		self.coinsP = self.coinsP[:-2*m]
		alpha = [gmpy2.t_mod_2exp(x,l) for x in r]
		alpha = [x.digits(2) for x in alpha]
		for i in range(0,m):
			if (len(alpha[i]) < l):
				alpha[i] = "".join(['0'*(l-len(alpha[i])),alpha[i]])
		self.alpha = alpha
		return z

	def obfuscate(self):
		m = self.m
		self.a2 = [0]*m
		self.b2 = [0]*m
		for i in range(0,m):
			r = self.obfuscation[i]
			self.a2[i] = self.a[i]+self.pubkey.encrypt(r[0])
			self.b2[i] = self.b[i]+self.pubkey.encrypt(r[1])
		return self.a2, self.b2

	def update(self,v):
		for i in range(0,self.m):
			r = self.obfuscation[i]
			self.mu[i] = v[i] + (self.t[i]-1)*r[0] + self.t[i]*(-r[1]) ### mu = mpz(mu*2**lf)

	def DGK_cloud(self,b):
		l = self.l
		m = self.m
		self.delta_A = [0]*m
		c_all = [[0]*l]*m;
		for k in range(0,m):
			beta = b[k]
			alpha = self.alpha[k]
			DGK_pubkey = self.DGK_pubkey
			delta_A = np.random.randint(0,2)
			self.delta_A[k] = delta_A
			prod = [0]*l
			c = [DGK_pubkey.raw_encrypt(0)]*l
			for i in range(0,l):
				if (int(alpha[i]) == 0):
					prod[i] = beta[i]
				else: prod[i] = testDGK.diff_encrypted(DGK_pubkey.raw_encrypt(1,self.coinsDGK.pop()),beta[i],DGK_pubkey)
				if (int(delta_A)==int(alpha[i])):
					if i==0: c[i] = DGK_pubkey.raw_encrypt(0,self.coinsDGK.pop())
					else: 
						for iter in range(0,i):
							c[i] = testDGK.add_encrypted(c[i],prod[iter],DGK_pubkey)
					if (int(delta_A) == 0):
						diff = testDGK.diff_encrypted(DGK_pubkey.raw_encrypt(1,self.coinsDGK.pop()),beta[i],DGK_pubkey)
						c[i] = testDGK.add_encrypted(c[i],diff,DGK_pubkey)
					else: c[i] = testDGK.add_encrypted(c[i],beta[i],DGK_pubkey)
			for i in range(0,l):
				if (int(delta_A)==int(alpha[i])):
					r = gmpy2.mpz_urandomb(gmpy2.random_state(),self.sigma+self.sigma)
					c[i] = testDGK.mul_sc_encrypted(c[i],r,DGK_pubkey)
				else: 
					c[i] = DGK_pubkey.raw_encrypt(gmpy2.mpz_urandomb(gmpy2.random_state(),self.sigma+self.sigma),self.coinsDGK.pop()) 
			c_all[k] = np.random.permutation(c)
		return c_all

	def compute_t(self,delta_B,zdivl):
		t = [0]*self.m
		for i in range(0,self.m):
			if (self.delta_A[i] == 1):
				t[i] = delta_B[i]
			else: t[i] = self.pubkey.encrypt(1) - delta_B[i]
			t[i] = zdivl[i] - self.pubkey.encrypt(mpz(gmpy2.f_div_2exp(self.r[i],self.l))) - t[i]
		self.t = t
		return t

def key(serialised):
	received_dict = json.loads(serialised)
	pk = received_dict['public_key']
	n = int(pk['n'])
	public_key = paillier.PaillierPublicKey(n=n)
	pk = received_dict['public_key_DGK']
	n = mpz(pk['n']); g = mpz(pk['g']); h = mpz(pk['h']); u = mpz(pk['u']);
	DGK_pubkey = testDGK.DGKpubkey(n=n,g=g,h=h,u=u)
	return public_key, DGK_pubkey

def send_encr_data(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	enc_with_one_pub_key = {}
	enc_with_one_pub_key = [str(x.ciphertext()) for x in encrypted_number_list]
	return json.dumps(enc_with_one_pub_key)

def send_plain_data(data):
	time.sleep(NETWORK_DELAY)
	return json.dumps([str(x) for x in data])

def recv_size(the_socket):
	#data length is packed into 4 bytes
	total_len=0;total_data=[];size=sys.maxsize
	size_data=sock_data=bytes([]);recv_size=4096
	while total_len<size:
		sock_data=the_socket.recv(recv_size)
		if not total_data:
			if len(sock_data)>4:
				size=struct.unpack('>i', sock_data[:4])[0]
				recv_size=size
				if recv_size>262144:recv_size=262144
				total_data.append(sock_data[4:])
			else:
				size_data+=sock_data

		else:
			total_data.append(sock_data)
		total_len=sum([len(i) for i in total_data ])
	return b''.join(total_data)

def get_enc_data(received_dict,pubkey):
	return [paillier.EncryptedNumber(pubkey, int(x)) for x in received_dict]

def send_DGK_data(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	encrypted = {}
	encrypted = [str(x) for x in encrypted_number_list]
	return json.dumps(encrypted)

def send_DGK_matrix(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	encrypted = {}
	encrypted = [[str(y) for y in x] for x in encrypted_number_list]
	return json.dumps(encrypted)

def get_DGK_data(received_dict):
	return [mpz(x) for x in received_dict]

def get_DGK_matrix(received_dict):
	return [[mpz(y) for y in x] for x in received_dict]

def main():
	# In order to run more instances consecutively, change n_final and m_final
	n_initial = 10
	m_initial = 5
	n_final = n_initial
	m_final = m_initial
	# Create a TCP/IP socket
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	print('Cloud: Socket successfully created')
	port = 10000
	# Bind the socket to the port
	localhost = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2] if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]]) if l][0][0]
	server_address = (localhost, port)
	print('Cloud: Starting up on {} port {}'.format(*server_address))
	sock.bind(server_address)

	# Listen for incoming connections
	sock.listen(1)      
	print('Cloud: Socket is listening')
	# Wait for a connection
	print('Cloud: Waiting for a connection')
	connection, client_address = sock.accept()	
	try:
		for n in range(n_initial,n_final+1,10):
			for m in range(m_initial,m_final+1,8):
				time.sleep(1)
				print('Cloud: Connection from', client_address)
				data = recv_size(connection)
				if data:
					pubkey,DGK_pubkey = key(data)
					fileA = "Data/A"+str(n)+"_"+str(m)+".txt"
					fileQ = "Data/Q"+str(n)+"_"+str(m)+".txt"
					fileb = "Data/b"+str(n)+"_"+str(m)+".txt"
					filec = "Data/c"+str(n)+"_"+str(m)+".txt"
					fileparam = "Data/param"+str(n)+"_"+str(m)+".txt"
					v = []; t = []
					agents = Agents(pubkey,fileb,filec)
					cloud = Cloud(pubkey,DGK_pubkey,fileA,fileQ,fileparam)
					b_A, c_A, m, n = agents.send_data()
					# Send m and K
					K = cloud.K
					data = send_plain_data([m,K])
					connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					# Iterations of the projected gradient descent
					print(n,m)
					for k in range(0,K):
						print(k)
						cloud.obfuscation = cloud.obfuscations[k]
						cloud.r = cloud.rn[k]		
						cloud.compute_grad(b_A,c_A)
						temp_mu = cloud.temporary_prec_mu()
						# Send temp_mu to the target
						data = send_encr_data(temp_mu)
						connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
						# Receive trunc_mu_r
						data = json.loads(recv_size(connection))
						trunc_mu_r = get_enc_data(data,pubkey)
						cloud.mu_bar = sum_encrypted_vectors(trunc_mu_r,[cloud.er.pop() for i in range(0,m)]) ### mu_bar = int(mu_bar*2**16)

						# Begin comparison procedure
						# Send z
						z = cloud.init_comparison_cloud()
						data = send_encr_data(z)
						connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
						# Receive b = bits of beta
						data = json.loads(recv_size(connection))
						b = get_DGK_matrix(data)
						c = cloud.DGK_cloud(b)
						# Send c
						serialized_data = send_DGK_matrix(c)
						connection.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
						# Receive delta_B, zvdil
						data = json.loads(recv_size(connection))
						merged = get_enc_data(data,pubkey)
						delta_B = merged[:m];zdivl = merged[m:]
						t = cloud.compute_t(delta_B,zdivl)
						# Send t,a2,b2
						a2,b2 = cloud.obfuscate()
						data = send_encr_data(t+a2+b2)
						connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
						# Receive v
						data = json.loads(recv_size(connection))
						v = get_enc_data(data,pubkey)
						cloud.update(v)
					x = cloud.compute_primal_optimum(c_A)
					# Send x
					data = send_encr_data(x)
					connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					# Wait for the target to finish its tasks -- this is for when consecutive problems are run
					data = json.loads(recv_size(connection))
					# Send 1 if the target should keep the connection open and 0 otherwise
					if(n==n_final and m==m_final):
						data = json.dumps(0)
						connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					else:
						data = json.dumps(1)
						connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))	
				else:
					print('Cloud: No data from', client_address)
					break

	finally:
	# Clean up the connection
		print('Cloud: Closing connection')
		connection.close()

if __name__ == '__main__':
	main()