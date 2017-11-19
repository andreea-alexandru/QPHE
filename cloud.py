#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
import numpy
import time
import testDGK
import random


DEFAULT_KEYSIZE = 512
DEFAULT_MSGSIZE = 32 # The message size of DGK has to be greater than 2*log2(DEFAULT_MSGSIZE), check u in DGK_pubkey
DEFAULT_SECURITYSIZE = 80
DEFAULT_PRECISION = int(DEFAULT_MSGSIZE/2) # of fractional bits
NETWORK_DELAY = 0 #0.01 # 10 ms

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
    return numpy.array([privkey.decrypt(i) for i in x])

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

def encrypt_vector_DGK(pubkey, x, coins=None):
	if (coins==None):
		return [pubkey.raw_encrypt(y) for y in x]
	else: return [pubkey.raw_encrypt(y,coins.pop()) for y in x]

def decrypt_vector_DGK(privkey, x):
    return numpy.array([privkey.raw_decrypt0(i) for i in x])

####### We take the convention that a number x < N/3 is positive, and that a number x > 2N/3 is negative. 
####### The range N/3 < x < 2N/3 allows for overflow detection.

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

class Agents:
    def __init__(self, pubkey,fileb,filec):
    	self.pubkey = pubkey
    	b_A = numpy.loadtxt(fileb, delimiter='\n')
    	c_A = numpy.loadtxt(filec, delimiter='\n')
    	self.enc_b_A = encrypt_vector(pubkey,fixed_point_vector(b_A))
    	self.enc_c_A = encrypt_vector(pubkey,fixed_point_vector(c_A))
    	self.m = numpy.size(b_A)
    	self.n = numpy.size(c_A)

    def send_data(self):
    	return self.enc_b_A, self.enc_c_A, self.m, self.n

class Cloud:
	def __init__(self, pubkey, DGK_pubkey, fileA, fileQ, fileparam):
		self.pubkey = pubkey
		self.DGK_pubkey = DGK_pubkey
		A = numpy.loadtxt(fileA, delimiter=',')
		Q = numpy.loadtxt(fileQ, delimiter=',')
		self.A = A
		m = numpy.size(A,0)
		At = A.transpose()
		self.At = At
		n = numpy.size(Q,0)
		self.Q = Q
		minvQ = -numpy.linalg.inv(Q)			### - Q^{-1}
		self.minvQ = fixed_point_matrix(minvQ)
		AinvQ = numpy.dot(A,minvQ)				### - AQ^{-1}
		self.AinvQ = AinvQ
		self.N = pubkey.n
		self.N2 = pubkey.nsquare
		self.m = m
		self.N_len = (self.N).bit_length()
		self.l = DEFAULT_MSGSIZE
		self.sigma = DEFAULT_SECURITYSIZE
		self.delta_A = [0]*m
		mu = numpy.zeros(m).astype(int)
		# mu = [x*(2**DEFAULT_PRECISION) for x in numpy.random.randint(-1,1, size=m).astype(int)]
		self.mu = encrypt_vector(pubkey, mu)
		param = numpy.loadtxt(fileparam, delimiter='\n')
		# self.K = param[0]
		# self.K = int(self.K)
		self.K = 30
		self.eta = param[1] 
		coeff_mu = numpy.identity(m) + self.eta*numpy.dot(AinvQ,At) ### I-\eta AQ^{-1}A'
		self.coeff_mu = fixed_point_matrix(coeff_mu)
		coeff_c = self.eta*AinvQ ### -\etaAQ^{-1}
		self.coeff_c = fixed_point_matrix(coeff_c)
		coeff_x = numpy.dot(minvQ,At)
		self.coeff_x = fixed_point_matrix(coeff_x)
		self.gen_rands()

	def gen_rands(self):
		m = self.m
		l = self.l
		sigma = self.sigma
		K = self.K
		random_state = gmpy2.random_state(seed)
		rn = [[[gmpy2.mpz_urandomb(random_state,l + sigma),gmpy2.mpz_urandomb(random_state,l + sigma)] for i in range(0,m)] for k in range(0,K)]
		self.obfuscations = rn
		rn = [[gmpy2.mpz_urandomb(random_state,l + sigma) for i in range(0,m)] for k in range(0,K)]
		self.rn = rn
		random_state = gmpy2.random_state(seed)
		coinsP = [gmpy2.mpz_urandomb(random_state,self.N_len-1) for i in range(0,4*m*K)]
		coinsP = [gmpy2.powmod(x, self.N, self.N2) for x in coinsP]
		self.coinsP = coinsP
		coinsDGK = [gmpy2.mpz_urandomb(random_state,2*sigma) for i in range(0,2*(l+1)*m*K)]
		coinsDGK = [gmpy2.powmod(self.DGK_pubkey.h, x, self.DGK_pubkey.n) for x in coinsDGK]
		self.coinsDGK = coinsDGK
		random_state = gmpy2.random_state(seed)
		rn = [2**DEFAULT_PRECISION*gmpy2.mpz_urandomb(random_state,sigma+l) for i in range(0,m*K)]
		# rn = [2**(l) for i in range(0,m*K)]
		self.fixedNoise = encrypt_vector(self.pubkey, rn, coinsP)
		er = encrypt_vector(self.pubkey,retrieve_fixed_point_vector([-x for x in rn]),coinsP)
		self.er = er

	def compute_grad(self,b_A,c_A):
		mu_bar = sum_encrypted_vectors(numpy.dot(self.coeff_mu,self.mu),numpy.dot(self.coeff_c,c_A))
		mu_bar = sum_encrypted_vectors(mu_bar,[x*fixed_point(-self.eta) for x in b_A])
		self.mu_bar = mu_bar ### \mu_bar*2^{32}

	def temporary_prec_mu(self):
		m = self.m
		pubkey = self.pubkey
		# for i in range(0,m):
			# r[i] = gmpy2.mpz_urandomb(gmpy2.random_state(),self.l + 1) ### FLOAT ONLY HAS 53 BITS OF PRECISION
		r = [self.fixedNoise.pop() for i in range(0,m)]
		temp_mu = sum_encrypted_vectors(self.mu_bar,r)		### mu_bar*2**(2f)+r
		# er = encrypt_vector(pubkey,fixed_point_vector([-x for x in r]), self.coinsP)
		return temp_mu
	
	def compute_primal_optimum(self,c_A):
		x = numpy.dot(self.coeff_x,self.mu)
		x = sum_encrypted_vectors(x,numpy.dot(self.minvQ,c_A))
		return x

	def randomize(self):
		m = self.m
		a = [0]*m
		b = [0]*m
		for i in range(0,m):
			a[i],b[i] = numpy.random.permutation([self.pubkey.encrypt(0),self.mu_bar[i]])
		self.a = a
		self.b = b
		#### ATTENTION AT NUMBERS OVER 2**l
		return a,b

	def init_comparison_cloud(self):
		m = self.m
		l = self.l
		pubkey = self.pubkey
		r = self.r
		a,b = self.randomize()
		z = diff_encrypted_vectors(b,a)
		z = sum_encrypted_vectors(z,encrypt_vector(pubkey,r,self.coinsP))
		z = sum_encrypted_vectors(z,encrypt_vector(pubkey,[2**l]*m,self.coinsP))
		alpha = [gmpy2.f_mod_2exp(x,l) for x in r]
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

	def update(self,v,t):
		for i in range(0,self.m):
			r = self.obfuscation[i]
			self.mu[i] = v[i] + (t[i]-1)*r[0] + t[i]*(-r[1]) ### mu = mpz(mu*2**16)

	def correct_d(self,d):
		for i in range(0,self.m):
			# if (self.r[i] < mpz(self.DGK_pubkey.n-1)/2):
			if (self.r[i] < mpz(self.N-1)/2):
				d[i] = self.DGK_pubkey.raw_encrypt(0,self.coinsDGK.pop())
		return d

	def DGK_cloud(self,b,d): ### OPTIMIZE THIS TO DO IT FOR THE WHOLE VECTOR
		l = self.l
		m = self.m
		c_all = [[0]*l]*m
		for i in range(0,m):
			beta = b[i]
			di = d[i]
			alpha = self.alpha[i]
			DGK_pubkey = self.DGK_pubkey
			delta_A = numpy.random.randint(0,1)
			self.delta_A[i] = delta_A
			s = 1 - 2*delta_A
			w = [0]*l
			c = [0]*l
			t_alpha = gmpy2.f_mod_2exp(self.r[i]-self.N,l) #### FIGURE HOW TO REPRESENT NEGATIVE BINARY NUMBERS
			if (t_alpha < 0): 
				neg = 1
			else: neg = 0
			t_alpha = abs(t_alpha).digits(2)
			if (len(t_alpha) < l and neg == 0):
				t_alpha = "".join(['0'*(l-len(t_alpha)),t_alpha])
			if (len(t_alpha) < l and neg == 1):
				t_alpha = "".join(['1','0'*(l-len(t_alpha)-1),t_alpha])
			# l-1 is the MSB
			for j in range(0,l):
				l1j = l-1-j
				if (int(alpha[l1j]) == 0):
					prod = beta[l1j]
				else: prod = testDGK.diff_encrypted(DGK_pubkey.raw_encrypt(1,self.coinsDGK.pop()),int(beta[l1j]),DGK_pubkey)
				if (int(alpha[l1j]) == int(t_alpha[l1j])):
					w[l1j] = prod
				else:
					w[l1j] = testDGK.diff_encrypted(prod,di,DGK_pubkey)
					w[l1j] = gmpy2.powmod(w[l1j],l,DGK_pubkey.n)
				# w[l1j] = testDGK.mul_sc_encrypted(w[l1j],2**j,DGK_pubkey)
			w3 = [0]*l
			for j in range(0,l):
				l1j = l-1-j
				if (j==l-1):
					w3[l1j] = DGK_pubkey.raw_encrypt(0,self.coinsDGK.pop())
				else: w3[l1j] = w[l1j-1]
				for k in range(j+2,l):
					w3[l1j] = testDGK.add_encrypted(w3[l1j],w[l-1-k],DGK_pubkey)
				w3[l1j] = testDGK.mul_sc_encrypted(w3[l1j],3,DGK_pubkey)
				c[l1j] = testDGK.add_encrypted(DGK_pubkey.raw_encrypt(s+int(alpha[l1j]),self.coinsDGK.pop()),
					testDGK.mul_sc_encrypted(di,int(alpha[l1j])-int(t_alpha[l1j]),DGK_pubkey),DGK_pubkey)
				c[l1j] = testDGK.diff_encrypted(c[l1j],beta[l1j],DGK_pubkey)
				c[l1j] = testDGK.add_encrypted(c[l1j],w3[l1j],DGK_pubkey)
				c2 = c
				r = gmpy2.mpz_urandomb(gmpy2.random_state(),self.sigma+self.sigma)
				c[l1j] = testDGK.mul_sc_encrypted(c[l1j],r,DGK_pubkey)
			c_all[i] = numpy.random.permutation(c)
		return c_all

	def compute_t(self,delta_B,zdivl):
		t = [0]*self.m
		for i in range(0,self.m):
			if (self.delta_A[i] == 1):
				t[i] = delta_B[i]
			else: t[i] = self.pubkey.encrypt(1) - delta_B[i]
			t[i] = zdivl[i] - self.pubkey.encrypt(mpz(gmpy2.f_div_2exp(self.r[i],self.l))) - t[i]
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
		# print('Received {!r}'.format(data))
		# n = 3
		# m = 2
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
						# Receive msgf
						data = json.loads(recv_size(connection))
						msgf = get_enc_data(data,pubkey)
						cloud.mu_bar = sum_encrypted_vectors(msgf,[cloud.er.pop() for i in range(0,m)]) ### mu_bar = int(mu_bar*2**16)
						# ### DEBUG: send cloud.mu_bar
						# data = send_encr_data(cloud.mu_bar)
						# connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
						# Send z
						z = cloud.init_comparison_cloud()
						data = send_encr_data(z)
						connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
						# Receive di
						data = json.loads(recv_size(connection))
						di = get_DGK_data(data)
						di = cloud.correct_d(di)
						# Receive b
						data = json.loads(recv_size(connection))
						b = get_DGK_matrix(data)
						c = cloud.DGK_cloud(b,di)
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
						cloud.update(v,t)
						# ### DEBUG: send cloud.mu
						# data = send_encr_data(cloud.mu)
						# connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
						# # cloud.mu = retrieve_fixed_point_vector(cloud.mu)		### No need, mu = int(mu*2**16)
					# 	# print(decrypt_vector(target.privkey,cloud.mu))
					x = cloud.compute_primal_optimum(c_A)
					# Send x
					data = send_encr_data(x)
					connection.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
					# Wait for the target to finish its tasks
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
# main()
if __name__ == '__main__':
	main()