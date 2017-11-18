#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
# from phe import util
import numpy
import time
import testDGK
import random
import os


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

class Target:
	def __init__(self, l=DEFAULT_MSGSIZE,t_DGK=2*DEFAULT_SECURITYSIZE):
		keypair = paillier.generate_paillier_keypair(n_length=DEFAULT_KEYSIZE)
		self.pubkey, self.privkey = keypair
		self.l = l
		self.t_DGK = t_DGK
		self.generate_DGK()


	def params(self,m,K):
		self.m = m
		self.K = K
		t2 = 2*self.t_DGK
		N_len = self.pubkey.n.bit_length()
		random_state = gmpy2.random_state(seed)
		coinsP = [gmpy2.mpz_urandomb(random_state,N_len-1) for i in range(0,5*m*K)]
		coinsP = [gmpy2.powmod(x, self.pubkey.n, self.pubkey.nsquare) for x in coinsP]
		self.coinsP = coinsP
		coinsDGK = [gmpy2.mpz_urandomb(random_state,t2) for i in range(0,(self.l+1)*m*K)]
		coinsDGK = [gmpy2.powmod(self.DGK_pubkey.h, x, self.DGK_pubkey.n) for x in coinsDGK]
		self.coinsDGK = coinsDGK
		self.delta_B = [0]*self.m

	def init_comparison_target(self,msg):
		l = self.l
		z = decrypt_vector(self.privkey,msg)
		z = [mpz(x) for x in z]
		self.z = z
		self.d = [(x >= mpz(self.pubkey.n-1)/2)*1 for x in z]
		beta = [gmpy2.f_mod_2exp(x,l) for x in z]
		beta = [x.digits(2) for x in beta]
		for i in range(0,self.m):
			if (len(beta[i]) < l):
				beta[i] = "".join(['0'*(l-len(beta[i])),beta[i]])
		self.beta = beta


	def generate_DGK(self):
		file = 'DGK_keys.txt'
		p,q,u,vp,vq,fp,fq,g,h = testDGK.loadkey(file)
		n = p*q
		self.DGK_pubkey = testDGK.DGKpubkey(n,g,h,u)
		self.DGK_privkey = testDGK.DGKprivkey(p,q,vp,self.DGK_pubkey)


	def DGK_target(self,c_all):
		l = self.l
		m = self.m
		for i in range(0,m):
			c = c_all[i]
			self.delta_B[i] = 0
			for j in range(0,l):
				if (int(self.DGK_privkey.raw_decrypt0(c[j])) == 0):
					self.delta_B[i] = 1
					break
		return encrypt_vector(self.pubkey,self.delta_B,self.coinsP), encrypt_vector(self.pubkey,[mpz(gmpy2.f_div_2exp(self.z[i],l)) for i in range(0,m)],self.coinsP)


	def choose(self,a,b):
		v = [0]*self.m
		for i in range(0,self.m):
			if self.t[i]==0: 
				v[i] = a[i] + self.pubkey.encrypt(0,self.coinsP.pop())
			else: v[i] = b[i] + self.pubkey.encrypt(0,self.coinsP.pop())
		return v, encrypt_vector(self.pubkey,[int(self.t[i]) for i in range(0,self.m)],self.coinsP)

def keys(pubkey,DGK_pubkey):
	pubkeys = {}
	pubkeys['public_key'] = {'n': pubkey.n}
	pubkeys['public_key_DGK'] = {'n': int(DGK_pubkey.n), 'g':int(DGK_pubkey.g),'h':int(DGK_pubkey.h), 'u':int(DGK_pubkey.u)}
	serialized_pubkeys = json.dumps(pubkeys)
	return serialized_pubkeys

def get_enc_data(received_dict,pubkey):
	return [paillier.EncryptedNumber(pubkey, int(x)) for x in received_dict]

def get_plain_data(data):
	return [int(x) for x in data]

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

def send_encr_data(encrypted_number_list):
	time.sleep(NETWORK_DELAY)
	encrypted = {}
	encrypted = [str(x.ciphertext()) for x in encrypted_number_list]
	return json.dumps(encrypted)

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

		start = time.time()
	# Create a TCP/IP socket
		sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		port = 10000

		# Connect the socket to the port where the server is listening
		localhost = [l for l in ([ip for ip in socket.gethostbyname_ex(socket.gethostname())[2] if not ip.startswith("127.")][:1], [[(s.connect(('8.8.8.8', 53)), s.getsockname()[0], s.close()) for s in [socket.socket(socket.AF_INET, socket.SOCK_DGRAM)]][0][1]]) if l][0][0]
		server_address = (localhost, port)
		print('Target: Connecting to {} port {}'.format(*server_address))
		sock.connect(server_address)			
		target = Target()
		pubkey = target.pubkey
		privkey = target.privkey
		DGK_pubkey = target.DGK_pubkey
		serialized_pubkeys = keys(pubkey,DGK_pubkey)
		cont = True
		off = time.time() - start
		try:		
			while cont:
				time.sleep(1)
				start = time.time()
				# Send public key
				# print('sending {!r}'.format(serialized_pubkeys))
				sock.sendall(struct.pack('>i', len(serialized_pubkeys))+serialized_pubkeys.encode('utf-8'))					
				# Receive m and K
				data = json.loads(recv_size(sock))
				m,K = get_plain_data(data)
				offline = off + time.time()-start
				start = time.time()
				l = target.l
				target.params(m,K)
				target.delta_B = [0]*m
				for k in range(0,K):
					# Receive temp_mu
					data = json.loads(recv_size(sock))
					temp_mu = get_enc_data(data,pubkey)
					temp_mu = decrypt_vector(privkey,temp_mu)
					#msgf = [mpz(x/(2**DEFAULT_PRECISION))for x in temp_mu] ### int((mu_bar*2**(2f) + r)/2**f)
					msgf = retrieve_fixed_point_vector([int(x) for x in temp_mu])
					# print(msgf) ### if r=0, supposed to be mu_bar*2**f
					msgf = encrypt_vector(target.pubkey,msgf,target.coinsP)
					# Send msgf		
					serialized_data = send_encr_data(msgf)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# ### Debug: receive cloud.mu_bar
					# data = json.loads(recv_size(sock))
					# mu_bar = get_enc_data(data,pubkey)
					# print(decrypt_vector(privkey,mu_bar))
					# Receive z
					data = json.loads(recv_size(sock))
					z = get_enc_data(data,pubkey)
					target.init_comparison_target(z)
					di = encrypt_vector_DGK(DGK_pubkey,target.d,target.coinsDGK)
					# Send di
					serialized_data = send_DGK_data(di)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					b = [[0]*l]*m
					b = [encrypt_vector_DGK(DGK_pubkey,[int(target.beta[i][j]) for j in range(0,l)],target.coinsDGK) for i in range(0,m)]
					# Send b
					serialized_data = send_DGK_matrix(b)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# Receive c
					data = json.loads(recv_size(sock))
					c = get_DGK_matrix(data)
					delta_B, zdivl = target.DGK_target(c)
					# Send delta_B, zdivl
					serialized_data = send_encr_data(delta_B+zdivl)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# Receive t,a2,bs
					data = json.loads(recv_size(sock))
					merged = get_enc_data(data,pubkey)
					t = merged[:m]; a2 = merged[m:2*m]; b2 = merged[2*m:]
					target.t = decrypt_vector(target.privkey,t)
					v,t = target.choose(a2,b2)
					# Send v
					serialized_data = send_encr_data(v)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
					# ### DEBUG: receive cloud.mu
					# data = json.loads(recv_size(sock))
					# mu = get_enc_data(data,pubkey)
					# print(decrypt_vector(privkey,mu))
				# Receive x
				data = json.loads(recv_size(sock))
				x = get_enc_data(data,pubkey)
				x = decrypt_vector(privkey,x)/2**(2*DEFAULT_PRECISION)
				print(x)
				end = time.time()
				sec = end-start
				print("%.2f" % sec)
				n = len(x)
				sys.stdout.flush()
				with open(os.path.abspath('Data/'+str(DEFAULT_KEYSIZE)+'_'+str(DEFAULT_MSGSIZE)+'_results_'+str(K)+'.txt'),'a+') as f: f.write("%d, %d, %.2f, %.2f\n" % (n,target.m,sec,offline))
				# Let the cloud know that it is ready
				data = json.dumps(1)
				sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
				# Receive 1 if to continue and 0 if to end
				data = json.loads(recv_size(sock))
				cont = bool(int(data))

		finally:
			print('Target: Closing socket')
			sock.close()				


# main()
if __name__ == '__main__':
	main()