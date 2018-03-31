#!/usr/bin/env python3

import socket
import sys,struct
import json
from gmpy2 import mpz
import paillier
# from phe import util
import numpy as np
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
    return [privkey.decrypt(i) for i in x]

def sum_encrypted_vectors(x, y):
	return [x[i] + y[i] for i in range(np.size(x))]

def diff_encrypted_vectors(x, y):
	return [x[i] - y[i] for i in range(len(x))] 

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

def retrieve_fp_vector(vec,prec=DEFAULT_PRECISION):
	return [retrieve_fp(x,prec) for x in vec]

def retrieve_fp_matrix(mat,prec=DEFAULT_PRECISION):
	return [retrieve_fp_vector(x,prec) for x in mat]


class Target:
	def __init__(self, l=DEFAULT_MSGSIZE,t_DGK=2*DEFAULT_SECURITYSIZE):
		keypair = paillier.generate_paillier_keypair(n_length=DEFAULT_KEYSIZE)
		self.pubkey, self.privkey = keypair
		# filepub = "Keys/pubkey"+str(DEFAULT_KEYSIZE)+".txt"
		# with open(filepub, 'r') as fin:
		# 	data=[line.split() for line in fin]
		# Np = int(data[0][0])
		# pubkey = paillier.PaillierPublicKey(n=Np)

		# filepriv = "Keys/privkey"+str(DEFAULT_KEYSIZE)+".txt"
		# with open(filepriv, 'r') as fin:
		# 	data=[line.split() for line in fin]
		# p = mpz(data[0][0])
		# q = mpz(data[1][0])
		# privkey = paillier.PaillierPrivateKey(pubkey, p, q)		
		# self.pubkey = pubkey; self.privkey = privkey
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
		beta = [gmpy2.t_mod_2exp(x,l) for x in z]
		beta = [x.digits(2) for x in beta]
		for i in range(0,self.m):
			if (len(beta[i]) < l):
				beta[i] = "".join(['0'*(l-len(beta[i])),beta[i]])
		self.beta = beta


	def generate_DGK(self):
		file = 'Keys/DGK_keys512_16.txt'
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
		db = encrypt_vector(self.pubkey,self.delta_B,self.coinsP[-m:]); z = encrypt_vector(self.pubkey,[mpz(gmpy2.f_div_2exp(self.z[i],l)) for i in range(0,m)],self.coinsP[-2*m:-m])
		self.coinsP = self.coinsP[:-2*m]
		return db,z


	def choose(self,a,b):
		m = self.m
		v = [0]*m
		for i in range(0,m):
			if self.t[i]==0: 
				v[i] = a[i] + self.pubkey.encrypt(0,self.coinsP.pop())
			else: v[i] = b[i] + self.pubkey.encrypt(0,self.coinsP.pop())
		return v

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
		lf = DEFAULT_PRECISION
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
				sock.sendall(struct.pack('>i', len(serialized_pubkeys))+serialized_pubkeys.encode('utf-8'))					
				# Receive m and K
				data = json.loads(recv_size(sock))
				m,K = get_plain_data(data)
				offline = off + time.time()-start
				start = time.time()
				l = target.l
				target.params(m,K)
				for k in range(0,K):
					# Receive temp_mu
					data = json.loads(recv_size(sock))
					temp_mu = get_enc_data(data,pubkey)
					temp_mu = decrypt_vector(privkey,temp_mu)
					msgf = fp_vector(temp_mu,-lf)
					msgf = encrypt_vector(target.pubkey,msgf,target.coinsP[-m:])
					target.coinsP = target.coinsP[:-m]
					# Send msgf		
					serialized_data = send_encr_data(msgf)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))

					# Begin comparison procedure
					# Receive z
					data = json.loads(recv_size(sock))
					z = get_enc_data(data,pubkey)
					target.init_comparison_target(z)
					b = [[0]*l]*m
					b = [encrypt_vector_DGK(DGK_pubkey,[int(target.beta[i][j]) for j in range(0,l)],target.coinsDGK[-(i+1)*l:-i*l] or target.coinsDGK[-l:]) for i in range(0,m)]
					target.coinsDGK = target.coinsDGK[:-l*m]
					# Send b = bits of beta
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
					v = target.choose(a2,b2)
					# Send v
					serialized_data = send_encr_data(v)
					sock.sendall(struct.pack('>i', len(serialized_data))+serialized_data.encode('utf-8'))
				# Receive x
				data = json.loads(recv_size(sock))
				x = get_enc_data(data,pubkey)
				x = retrieve_fp_vector(decrypt_vector(privkey,x),2*lf)
				print(["%.8f"% i for i in x])
				end = time.time()
				sec = end-start
				print("%.2f" % sec)
				n = len(x)
				sys.stdout.flush()
				with open(os.path.abspath('Results/'+str(DEFAULT_KEYSIZE)+'_'+str(DEFAULT_MSGSIZE)+'_results_'+str(K)+'.txt'),'a+') as f: f.write("%d, %d, %.2f, %.2f\n" % (n,target.m,sec,offline))
				# Let the cloud know that it is ready
				data = json.dumps(1)
				sock.sendall(struct.pack('>i', len(data))+data.encode('utf-8'))
				# Receive 1 if to continue and 0 if to end
				data = json.loads(recv_size(sock))
				cont = bool(int(data))

		finally:
			print('Target: Closing socket')
			sock.close()				


if __name__ == '__main__':
	main()