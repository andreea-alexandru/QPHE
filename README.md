# QPHE

Quick and dirty code for implementing dual gradient ascent for a strongly convex optimization problem on encrypted data in Python. More details can be found in the paper "Cloud-based Quadratic Optimization with Partially Homomorphic Encryption" https://arxiv.org/abs/1809.02267.

This code is intended to give an idea on the execution speed and accuracy, and not for application that require highly optimized code.

The file runSim.py and runSS.py is the file that runs a thread for each of the two computing parties.

The matrices and vectors for the optimization problem as in the paper have to be written in files in the folder Data, with the name having name-of-matrix_number-of-variables_number-of-constrains.txt. The file write_matrices.py can be used to generate random instances for given dimensions.

If the Paillier and DGK keys are written in files in the folder Keys, then the code will read them. Otherwise, if the key files do not exist, then the code will generate the keys.
