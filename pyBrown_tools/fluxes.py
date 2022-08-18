import numpy as np

def compute_flux_single(r0, r1, plane_normal_vector, plane_point):

	f0 = np.dot( plane_normal_vector, (r0 - plane_point) ) > 0.0
	f1 = np.dot( plane_normal_vector, (r1 - plane_point) ) > 0.0

	if f0 == f1: return 0
	elif f0: return -1
	else: return 1