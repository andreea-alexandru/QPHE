#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys
 
"""Run the dual gradient ascent for a strongly convex optimization
	problem on encrypted data in a Two-Server architecture. 
	The two servers are emulated by different threads.
	See https://arxiv.org/abs/1809.02267 for more details."""

if __name__ == '__main__':
	proc1 = Popen(['python','cloud.py'],stdout=PIPE)
	proc2 = Popen(['python','target.py'])
	# stdout_value = proc2.communicate()[0]
	stdout_value = proc1.communicate()[0]

