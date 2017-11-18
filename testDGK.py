#!/usr/bin/env python3
import random
from phe.util import invert, powmod, getprimeover, isqrt
from gmpy2 import mpz


DEFAULT_KEYSIZE = 512
DEFAULT_SECURITYSIZE = 160

try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False

# pubkey = (n,g,h,u)
# privkey = (p,q,v)

class DGKpubkey:
	def __init__(self,n,g,h,u,t=DEFAULT_SECURITYSIZE):
		self.n = n
		self.g = g;
		self.h = h;
		self.u = u;
		self.t = t;

	def raw_encrypt(self, plaintext, r_value=None):
		"""DGK encryption of a positive integer plaintext .
		You probably should be using :meth:`encrypt` instead, because it handles positive and negative ints and floats.

		Args:
			plaintext (int): a positive integer < :attr:`n` to be DGK
			encrypted. Typically this is an encoding of the actual 
			number you want to encrypt.
			r_value (int): obfuscator for the ciphertext; by default (i.e.
			r_value is None), a random value of 2t bits is used.
		Returns:
			int: DGK encryption of plaintext.

		Raises:
			TypeError: if plaintext is not an int.
		"""
		if not isinstance(plaintext, int):
			raise TypeError('Expected int type plaintext but got: %s' %
                            type(plaintext))

		nude_ciphertext = powmod(self.g, plaintext, self.n)
		# nude_ciphertext = self.g**plaintext

<<<<<<< HEAD
		# r = r_value or self.get_random_lt_2t()
		# obfuscator = powmod(self.h, r, self.n)
		r = r_value or powmod(self.h, self.get_random_lt_2t(), self.n) # Pass the precomputed obfuscator
		obfuscator = r		
=======
		r = r_value or self.get_random_lt_2t()
		obfuscator = powmod(self.h, r, self.n)
>>>>>>> ee98cce2f2beaf0080cf452f430c6f220d17f086

		return (nude_ciphertext * obfuscator) % self.n

	def get_random_lt_2t(self):
		"""Return a cryptographically random number less than :attr:`n`"""
		# return random.SystemRandom().randrange(1, 2**(2*self.t))
		t2 = 2*DEFAULT_SECURITYSIZE
		return random.SystemRandom().randrange(1, 2**t2)

class DGKprivkey:
	def __init__(self,p,q,v,pubkey):
		self.p = p
		self.q = q
		self.v = v
		self.pubkey = pubkey

	def raw_decrypt0(self, ciphertext):
		"""Decrypt raw ciphertext and return raw plaintext.

		Args:
			ciphertext (int): (usually from :meth:`EncryptedNumber.ciphertext()`)
				that is to be DGK decrypted.

		Returns:
			int: DGK decryption of ciphertext. This is a positive
			integer < :attr:`public_key.n`.

		Raises:
			TypeError: if ciphertext is not an int.

		Only see if the plaintext is 0 by checking cyphertext^v mod n == 1
		"""
		# if not isinstance(ciphertext, int):
		# 	raise TypeError('Expected ciphertext to be an int, not: %s' %
		# 		type(ciphertext))

		# c = powmod(ciphertext, self.v, self.pubkey.n)
		c = powmod(ciphertext, self.v, self.p)
		if c==1:
			return 0
		else: return 1



def loadkey(file):
	with open(file, 'r') as fin:
		data=[line.split() for line in fin]
	# data = numpy.loadtxt(file, delimiter=' ')
	p = data[0][0]
	q = data[1][0]
	u = data[2][0]
	vp = data[3][0]
	vq = data[4][0]
	fp = data[5][0]
	fq = data[6][0]	
	g = data[7][0]
	h = data[8][0]		
	p = mpz(p)
	q = mpz(q)
	u = mpz(u)
	vp = mpz(vp)
	vq = mpz(vq)
	fp = mpz(fp)
	fq = mpz(fq)	
	g = mpz(g)
	h = mpz(h)					
	return p,q,u,vp,vq,fp,fq,g,h

def add_encrypted(a,b,pubkey):
	return gmpy2.t_mod(gmpy2.mul(a,b),pubkey.n)

def diff_encrypted(a,b,pubkey):
	return add_encrypted(a, gmpy2.invert(b,pubkey.n), pubkey)

def mul_sc_encrypted(a,b,pubkey):
	return gmpy2.powmod(a,b,pubkey.n)


# def main():
# 	file = 'DGK_keys.txt'
# 	p,q,u,vp,vq,fp,fq,g,h = loadkey(file)
# 	n = p*q
# 	v = vp*vq

# 	pubkey = DGKpubkey(n,g,h,u)
# 	privkey = DGKprivkey(p,q,vp,pubkey)
# 	plaintext = 0
# 	ciphertext = pubkey.raw_encrypt(plaintext)
# 	print(privkey.raw_decrypt0(ciphertext))
# 	plaintext2 = 1
# 	ciphertext2 = pubkey.raw_encrypt(plaintext2)
# 	print(privkey.raw_decrypt0(ciphertext2))
# 	sume = diff_encrypted(ciphertext,ciphertext2,pubkey)
# 	print(privkey.raw_decrypt0(sume))

# main()