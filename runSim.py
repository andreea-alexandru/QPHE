#!/usr/bin/env python3
from subprocess import Popen, PIPE
import sys
 
if __name__ == '__main__':
	proc1 = Popen(['python','cloud.py'],stdout=PIPE)
	proc2 = Popen(['python','target.py'])
	# stdout_value = proc2.communicate()[0]
	stdout_value = proc1.communicate()[0]


