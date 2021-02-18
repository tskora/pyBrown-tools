# pyBrown is a bound of tools useful for Brownian and Stokesian dynamics simulations
# Copyright (C) 2018  Tomasz Skora (tskora@ichf.edu.pl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses.

import math
import numpy as np
from scipy.special import erfc

class Hydrodynamics():

	def __init__(self):

		return None

class Hydrodynamics_RPY(Hydrodynamics):

	def __init__(self):

		return self.super()

def Mii_rpy(a):

	return np.identity(3) / ( 6 * np.pi * a )

def Mij_rpy(ai, aj, pointer):

	Rh_larger = max( ai, aj )
	Rh_smaller = min( ai, aj )

	dist2 = pointer[0]**2 + pointer[1]**2 + pointer[2]**2
	dist = math.sqrt( dist2 )
	outer = np.outer(pointer, pointer)/dist2

	aij2 = ai**2 + aj**2

	if dist > ( ai + aj ):
	
		coef_1 = 1.0 / ( 8 * np.pi * dist )
		coef_2 = 1.0 + aij2 / ( 3 * dist2 )
		coef_3 = 1.0 - aij2 / dist2
	
		answer = coef_2 * np.identity(3)
		answer += coef_3 * outer
		answer *= coef_1
	
		return answer
        
	elif dist <= ( Rh_larger - Rh_smaller ):

		return np.identity(3) / ( 6 * np.pi * Rh_larger )

	else:

		dist3 = dist * dist2

		coef_1 = 1.0 / ( 6 * np.pi * ai * aj )
		coef_2 = 16 * dist3 * ( ai + aj )
		coef_3 = (ai - aj)**2 + 3 * dist2
		coef3 *= coef3
		coef_4 = ( coef_2 - coef_3 ) / ( 32 * dist3 )
		coef_5 = 3 * ( (ai - aj)**2 - dist2 )**2
		coef_6 = coef_5 / ( 32 * dist3 )

		answer = coef_4 * np.identity(3)
		answer += coef_6 * outer
		answer *= coef_1
        
		return answer

def M_rpy(beads, pointers):

	M = [ [ None for j in range( len(beads) ) ] for i in range( len(beads) ) ]

	for i, bi in enumerate(beads):

		M[i][i] = Mii_rpy(bi.a)

		for j in range(i):

			bj = beads[j]

			M[i][j] = Mij_rpy(bi.a, bj.a, pointers[i][j])
			M[j][i] = np.transpose( M[i][j] )

	return np.block(M)

# def M_rpy_smith(beads, points, box_length, alpha, m, n):

# 	M = [ [ None for j in range( len(beads) ) ] for i in range( len(beads) ) ]

# 	for i, pi in enumerate(beads):

# 		M[i][i] = Mii_rpy_smith(pi, box_length, alpha, m, n)

# 		for j in range(i):

# 			pj = ps[j]

# 			M[i][j] = Mij_rpy_smith(pi, pj, box_length, alpha, m, n)
# 			M[j][i] = np.transpose( M[i][j] )

# 	return np.block(M)

# def O(r):

# 	dist = math.sqrt( r[0]**2 + r[1]**2 + r[2]**2 )

# 	return ( np.identity(3) + np.outer( r / dist, r / dist ) ) / dist

# def Q(r):

# 	dist = math.sqrt( r[0]**2 + r[1]**2 + r[2]**2 )

# 	return ( np.identity(3) - 3.0 * np.outer( r / dist, r / dist ) ) / dist**3

# def Oii_pbc_smith(p, box_length, alpha, m, n ):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) if not (mi==0 and mj==0 and mk==0 ) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		mlength = math.sqrt( mvec[0]**2 + mvec[1]**2 + mvec[2]**2 )

# 		answer += erfc( alpha * mlength ) * O(mvec)

# 		answer += 2.0 * alpha / math.sqrt(np.pi) * math.exp( - alpha**2 * mlength**2 ) * np.outer(mvec/mlength, mvec/mlength)

# 	for nvec in ns:

# 		nlength = math.sqrt( nvec[0]**2 + nvec[1]**2 + nvec[2]**2 )

# 		mult = np.identity(3) - ( 1 + np.pi**2 * nlength**2 / alpha**2 ) * np.outer(nvec/nlength, nvec/nlength)

# 		answer += 2.0 / ( np.pi * nlength**2 ) * math.exp( - np.pi**2 * nlength**2 / alpha**2 ) * mult

# 	return answer - 3.0 * alpha * p.a / ( 2.0 * math.sqrt(np.pi) * box_length ) * np.identity(3)

# def Oij_pbc_smith(sigma, alpha, m, n):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		msvec = mvec + sigma

# 		mslength = math.sqrt( msvec[0]**2 + msvec[1]**2 + msvec[2]**2 )

# 		answer += erfc( alpha * mslength ) * O( mvec + sigma )

# 		answer += 2.0 * alpha / math.sqrt(np.pi) * math.exp( - alpha**2 * mslength**2 ) * np.outer(mvec+sigma, mvec+sigma)/mslength**2

# 	for nvec in ns:

# 		nlength = math.sqrt( nvec[0]**2 + nvec[1]**2 + nvec[2]**2 )

# 		mult = 2.0 / ( np.pi * nlength**2 ) * math.exp( - np.pi**2 * nlength**2 / alpha**2 ) * np.exp( 2 * np.pi * 1j * np.dot(nvec, sigma)  )

# 		mult_real = mult.real

# 		answer += mult_real * ( np.identity(3) - (1 + np.pi**2 * nlength**2 / alpha**2 ) * np.outer(nvec / nlength, nvec / nlength) )

# 	return answer

# def Qii_pbc_smith(p, box_length, alpha, m, n):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) if not (mi==0 and mj==0 and mk==0 ) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		mlength = math.sqrt( mvec[0]**2 + mvec[1]**2 + mvec[2]**2 )

# 		mult = erfc( alpha * mlength ) + 2.0 * alpha / math.sqrt(np.pi) * mlength * math.exp( -alpha**2 * mlength**2 )

# 		answer += mult * Q(mvec)

# 		answer -= 4.0 * alpha**3 / math.sqrt(np.pi) * math.exp(-alpha**2 * mlength**2) * np.outer(mvec/mlength, mvec/mlength)

# 	for nvec in ns:

# 		nlength = math.sqrt( nvec[0]**2 + nvec[1]**2 + nvec[2]**2 )

# 		answer += 4.0 * np.pi * math.exp( -np.pi**2 * nlength**2 / alpha**2) * np.outer(nvec/nlength, nvec/nlength)

# 	return answer - 1.0 / ( 3.0 * math.sqrt(np.pi) ) * (alpha * p.a / box_length)**3 * np.identity(3)

# def Qij_pbc_smith( sigma, alpha, m, n ):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		msvec = mvec + sigma

# 		mslength = math.sqrt( msvec[0]**2 + msvec[1]**2 + msvec[2]**2 )

# 		mult = erfc( alpha * mslength ) + 2.0 * alpha / np.sqrt( np.pi ) * mslength * math.exp(- alpha**2 * mslength**2)

# 		answer += mult * Q(mvec + sigma)

# 		answer -= 4.0 * alpha**3 / np.sqrt(np.pi) * math.exp( -alpha**2 * mslength**2 ) * np.outer(mvec+sigma,mvec+sigma)/mslength**2

# 	for nvec in ns:

# 		nlength = math.sqrt( nvec[0]**2 + nvec[1]**2 + nvec[2]**2 )

# 		addi = 4.0 * np.pi * math.exp( -np.pi**2 * nlength**2 / alpha**2 ) * np.exp(2*np.pi*1j*np.dot(nvec,sigma)) * np.outer(nvec/nlength, nvec/nlength)

# 		addi_real = addi.real

# 		answer += addi_real

# 	return answer

# def Mii_rpy_smith(p, box_length, alpha, m, n):

# 	coef1 = 1.0 / ( 6 * np.pi * p.a )

# 	coef2 = 3.0 * p.a / ( 4.0 * box_length )

# 	coef3 = ( p.a / box_length )**3 / 2.0

# 	comp1 = np.identity(3)

# 	comp2 = Oii_pbc_smith( p, box_length, alpha, m, n )

# 	comp3 = Qii_pbc_smith( p, box_length, alpha, m, n )

# 	return coef1 * ( comp1 + coef2 * comp2 + coef3 * comp3 )

# def Mij_rpy_smith(pi, pj, box_length, alpha, m, n):

# 	sigma = np.array( (pj.r - pi.r) / box_length, float )

# 	coef1 = 1.0 / ( 6 * np.pi * pi.a )

# 	coef2 = 3.0 * pi.a / ( 4.0 * box_length )

# 	coef3 = ( pi.a / box_length )**3 / 2.0

# 	comp1 = Oij_pbc_smith( sigma, alpha, m, n )

# 	comp2 = Qij_pbc_smith( sigma, alpha, m, n )

# 	return coef1 * ( coef2 * comp1 + coef3 * comp2 )

# def Bij_rpy_smith(pi, pj):

# 	rij = pj.r - pi.r

# 	rijlength = math.sqrt( rij[0]**2 + rij[1]**2 + rij[2]**2 )

# 	coef1 = 1.0 / ( 8 * np.pi * rijlength )

# 	coef2 = 1.0 + ( pi.a**2 + pj.a**2 ) / ( 3.0 * rijlength**2 )

# 	coef3 = 1.0 - ( pi.a**2 + pj.a**2 ) / rijlength**2

# 	return coef1 * ( coef2 * np.identity(3) + coef3 * np.outer(rij, rij) / rijlength**2 )

# def M_rpy_smith(ps, box_length, alpha, m, n):

# 	M = [ [ None for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		M[i][i] = Mii_rpy_smith(pi, box_length, alpha, m, n)

# 		for j in range(i):

# 			pj = ps[j]

# 			M[i][j] = Mij_rpy_smith(pi, pj, box_length, alpha, m, n)
# 			M[j][i] = np.transpose( M[i][j] )

# 	return np.block(M)

# def B_rpy_smith(ps):

# 	B = [ [ None for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		B[i][i] = np.zeros((3,3))

# 		for j in range(i):

# 			pj = ps[j]

# 			B[i][j] = Bij_rpy_smith(pi, pj)
# 			B[j][i] = np.transpose( B[i][j] )

# 	return np.block(B)

# def X_f_poly(l, rank):

# 	if rank == 0: return 1
# 	if rank == 1: return 3 * l
# 	if rank == 2: return 9 * l
# 	if rank == 3: return -4 * l + 27 * l**2 - 4 * l**3
# 	if rank == 4: return -24 * l + 81 * l**2 + 36 * l**3
# 	if rank == 5: return 72 * l**2 + 243 * l**3 + 72 * l**4
# 	if rank == 6: return 16 * l + 108 * l**2 + 281 * l**3 + 648 * l**4 + 144 * l**5
# 	if rank == 7: return 288 * l**2 + 1620 * l**3 + 1515 * l**4 + 1620 * l**5 + 288 * l**6
# 	if rank == 8: return 576 * l**2 + 4848 * l**3 + 5409 * l**4 + 4524 * l**5 + 3888 * l**6 + 576 * l**7
# 	if rank == 9: return 1152 * l**2 + 9072 * l**3 + 14752 * l**4 + 26163 * l**5 + 14752 * l**6 + 9072 * l**7 + 1152 * l**8
# 	if rank == 10: return 2304 * l**2 + 20736 * l**3 + 42804 * l**4 + 115849 * l**5 + 76176 * l**6 + 39264 * l**7 + 20736 * l**8 + 2304 * l**9
# 	if rank == 11: return 4608 * l**2 + 46656 * l**3 + 108912 * l**4 + 269100 * l**5 + 319899 * l**6 + 269100 * l**7 + 108912 * l**8 + 46656 * l**9 + 4608 * l**10
# 	else: return None

# def Y_f_poly(l, rank):

# 	if rank == 0: return 1
# 	if rank == 1: return 3 / 2 * l
# 	if rank == 2: return 9 / 4 * l
# 	if rank == 3: return 2 * l + 27 / 8 * l**2 + 2 * l**3
# 	if rank == 4: return 6 * l + 81 / 16 * l**2 + 18 * l**3
# 	if rank == 5: return 63 / 2 * l**2 + 243 / 32 * l**3 + 63 / 2 * l**4
# 	if rank == 6: return 4 * l + 54 * l**2 + 1241 / 64 * l**3 + 81 * l**4 + 72 * l**5
# 	if rank == 7: return 144 * l**2 + 1053 / 8 * l**3 + 19083 / 128 * l**4 + 1053 / 8 * l**5 + 144 * l**6
# 	if rank == 8: return 279 * l**2 + 4261 / 8 * l**3 + 126369 / 256 * l**4 - 117 / 8 * l**5 + 648 * l**6 + 288 * l**7
# 	if rank == 9: return 576 * l**2 + 1134 * l**3 + 60443 / 32 * l**4 + 766179 / 512 * l**5 + 60443 / 32 * l**6 + 1134 * l**7 + 576 * l**8
# 	if rank == 10: return 1152 * l**2 + 7857 / 4 * l**3 + 98487 / 16 * l**4 + 10548393 / 1024 * l**5 + 67617 / 8 * l**6 - 351 / 2 * l**7 + 3888 * l**8 + 1152 * l**9
# 	if rank == 11: return 2304 * l**2 + 7128 * l**3 + 22071 / 2 * l**4 + 2744505 / 128 * l**5 + 95203835 / 2048 * l**6 + 2744505 / 128 * l**7 + 22071 / 2 * l**8 + 7128 * l**9 + 2304 * l**10
# 	else: return None

# #-------------------------------------------------------------------------------

# def X_g_poly(l, rank):

# 	if rank == 1: return 2 * l**2 * ( 1 + l )**(-3)
# 	if rank == 2: return 1 / 5 * l * ( 1 + 7 * l + l**2 ) * ( 1 + l )**(-3)
# 	if rank == 3: return 1 / 42 * ( 1 + 18 * l - 29 * l**2 + 18 * l**3 + l**4 ) * ( 1 + l )**(-3)
# 	else: return None

# #-------------------------------------------------------------------------------

# def Y_g_poly(l, rank):

# 	if rank == 2: return 4 / 15 * l * ( 2 + l + 2 * l**2 ) * ( 1 + l )**(-3)
# 	if rank == 3: return 2 / 375 * ( 16 - 45 * l + 58 * l**2 - 45 * l**3 + 16 * l**4 ) * ( 1 + l )**(-3)
# 	else: return None

# #-------------------------------------------------------------------------------

# def XA11(s, l):

# 	answer = 0.0

# 	answer += X_g_poly(l, 1) * ( 1 - 4 * s**(-2) )**(-1)

# 	answer -= X_g_poly(l, 2) * np.log( 1 - 4 * s**(-2) )

# 	answer -= X_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( 1 - 4 * s**(-2) )

# 	answer += X_f_poly(l, 0) - X_g_poly(l, 1)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 0 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * X_f_poly(l, m) - X_g_poly(l, 1) )

# 		answer += mult * ( 4 * m**(-1) * m1**(-1) * X_g_poly(l, 3) - 2 * m**(-1) * X_g_poly(l, 2) )

# 	return answer

# #-------------------------------------------------------------------------------

# def YA11(s, l):

# 	answer = 0.0

# 	answer -= Y_g_poly(l, 2) * np.log( 1 - 4 * s**(-2) )

# 	answer -= Y_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( 1 - 4 * s**(-2) )

# 	answer += Y_f_poly(l, 0)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 0 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * Y_f_poly(l, m) - 2 * m**(-1) * Y_g_poly(l, 2) )

# 		answer += mult * 4 * m**(-1) * m1**(-1) * Y_g_poly(l, 3)

# 	return answer

# #-------------------------------------------------------------------------------

# def XA12(s, l):

# 	answer = 0.0

# 	answer += 2 * s**(-1) * X_g_poly(l, 1) * ( 1 - 4 * s**(-2) )**(-1)

# 	answer += X_g_poly(l, 2) * np.log( ( s + 2 ) / ( s - 2 ) )

# 	answer += X_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( ( s + 2 ) / ( s - 2 ) ) + 4 * X_g_poly(l, 3) * s**(-1)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 1 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * X_f_poly(l, m) - X_g_poly(l, 1) )

# 		answer += mult * ( 4 * m**(-1) * m1**(-1) * X_g_poly(l, 3) - 2 * m**(-1) * X_g_poly(l, 2) )

# 	divisor = -1 / 2 * ( 1 + l )

# 	return answer / divisor

# #-------------------------------------------------------------------------------

# def YA12(s, l):

# 	answer = 0.0

# 	answer += Y_g_poly(l, 2) * np.log( ( s + 2 ) / ( s - 2 ) )

# 	answer += Y_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( ( s + 2 ) / ( s - 2 ) )

# 	answer += 4 * Y_g_poly(l, 3) * s**(-1)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 1 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * Y_f_poly(l, m) - 2 * m**(-1) * Y_g_poly(l, 2) )

# 		answer += mult * ( 4 * m**(-1) * m1**(-1) * Y_g_poly(l, 3) )

# 	divisor = -1 / 2 * ( 1 + l )

# 	return answer / divisor

# #-------------------------------------------------------------------------------

# def R_jeffrey(p1, p2):

# 	r = p2.r - p1.r

# 	dist = math.sqrt( r[0]**2 + r[1]**2 + r[2]**2 )

# 	s = 2 * dist / ( p1.a + p2.a )

# 	l = p2.a / p1.a

# 	R = [ [ None , None ], [ None, None ] ]

# 	R[0][0] = XA11(s, l) * np.outer(r/dist, r/dist) + YA11(s, l) * ( np.identity(3) - np.outer(r/dist, r/dist) )

# 	R[1][1] = XA11(s, 1/l) * np.outer(r/dist, r/dist) + YA11(s, 1/l) * ( np.identity(3) - np.outer(r/dist, r/dist) )

# 	R[0][1] = XA12(s, l) * np.outer(r/dist, r/dist) + YA12(s, l) * ( np.identity(3) - np.outer(r/dist, r/dist) )

# 	R[1][0] = XA12(s, 1/l) * np.outer(r/dist, r/dist) + YA12(s, 1/l) * ( np.identity(3) - np.outer(r/dist, r/dist) )

# 	return 3 * np.pi * ( p1.a + p2.a ) * np.block(R)

# #-------------------------------------------------------------------------------

# def R_lub_corr(ps):

# 	corr = [ [ np.zeros((3,3)) for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		for j in range(i):

# 			pj = ps[j]

# 			nf2b = R_jeffrey( pi, pj )

# 			ff2b = np.linalg.inv( M_rpy( [pi, pj] ) )

# 			lub_corr = nf2b - ff2b

# 			corrii = lub_corr[0:3,0:3]
# 			corrjj = lub_corr[3:6,3:6]
# 			corrij = lub_corr[3:6,0:3]

# 			corr[i][i] += corrii
# 			corr[j][j] += corrjj
# 			corr[i][j] += corrij
# 			corr[j][i] += np.transpose( corrij )

# 	return np.block(corr)

#####

# #-------------------------------------------------------------------------------

# class Particle:

# 	def __init__(self, r, Rh):

# 		self.r = np.array(r, float)

# 		self.Rh = Rh

# 	#---------------------------------------------------------------------------

# 	def __eq__(self, p):

# 		if isinstance( p, Particle ):
# 			return ( np.all( self.r == p.r ) ) and ( self.Rh == p.Rh )
# 		return False

# 	#---------------------------------------------------------------------------

# 	def distance(self, p):

# 		squared_differences = (self.r - p.r)**2
        
# 		return np.sqrt( np.sum( squared_differences ) )

# 	#---------------------------------------------------------------------------

# 	def outer_norm(self, p):

# 		difference = self.r - p.r
        
# 		return np.outer(difference, difference) / self.distance(p)**2

# 	#---------------------------------------------------------------------------

# #-------------------------------------------------------------------------------

# def Mii_rpy(p):

# 	return np.identity(3) / ( 6 * np.pi * p.Rh * VISCOSITY )

# #-------------------------------------------------------------------------------

# def Mij_rpy(pi, pj):

# 	Rh_larger = max( pi.Rh, pj.Rh )
# 	Rh_smaller = min( pi.Rh, pj.Rh )
	
# 	if pi.distance(pj) > ( pi.Rh + pj.Rh ):
	
# 		coef_1 = 1.0 / ( 8 * np.pi * pi.distance(pj) * VISCOSITY )
# 		coef_2 = 1.0 + ( pi.Rh**2 + pj.Rh**2) / ( 3 * pi.distance(pj)**2 )
# 		coef_3 = 1.0 - ( pi.Rh**2 + pj.Rh**2) / pi.distance(pj)**2
	
# 		answer = coef_2 * np.identity(3)
# 		answer += coef_3 * pi.outer_norm(pj)
# 		answer *= coef_1
	
# 		return answer
        
# 	elif pi.distance(pj) <= ( Rh_larger - Rh_smaller ):

# 		return np.identity(3) / ( 6 * np.pi * Rh_larger * VISCOSITY )

# 	else:

# 		coef_1 = 1.0 / ( 6 * np.pi * pi.Rh * pj.Rh * VISCOSITY )
# 		coef_2 = 16 * pi.distance(pj)**3 * ( pi.Rh + pj.Rh )
# 		coef_3 = ( (pi.Rh - pj.Rh )**2 + 3 * pi.distance(pj)**2 )**2
# 		coef_4 = ( coef_2 - coef_3 ) / ( 32 * pi.distance(pj)**3 )
# 		coef_5 = 3 * ( (pi.Rh - pj.Rh)**2 - pi.distance(pj)**2 )**2
# 		coef_6 = coef_5 / ( 32 * pi.distance(pj)**3 )

# 		answer = coef_4 * np.identity(3)
# 		answer += coef_6 * pi.outer_norm(pj)
# 		answer *= coef_1
        
# 		return answer

# #-------------------------------------------------------------------------------

# def M_rpy(ps):

# 	M = [ [ None for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		M[i][i] = Mii_rpy(pi)

# 		for j in range(i):

# 			pj = ps[j]

# 			M[i][j] = Mij_rpy(pi, pj)
# 			M[j][i] = np.transpose( M[i][j] )

# 	return np.block(M)

# #-------------------------------------------------------------------------------

# def length(r):

# 	return np.sqrt( np.sum( np.array(r, float)**2 ) )

# #-------------------------------------------------------------------------------

# def O(r):

# 	return ( np.identity(3) + np.outer( r / length(r), r / length(r) ) ) / length(r)

# #-------------------------------------------------------------------------------

# def Q(r):

# 	return ( np.identity(3) - 3.0 * np.outer( r / length(r), r / length(r) ) ) / length(r)**3

# #-------------------------------------------------------------------------------

# def Oii_pbc_smith( Rh, L, alpha, m, n ):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) if not (mi==0 and mj==0 and mk==0 ) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		answer += erfc( alpha * length(mvec) ) * O(mvec)

# 		answer += 2.0 * alpha / np.sqrt(np.pi) * np.exp( - alpha**2 * length(mvec)**2 ) * np.outer(mvec/length(mvec), mvec/length(mvec))

# 	for nvec in ns:

# 		mult = np.identity(3) - ( 1 + np.pi**2 * length(nvec)**2 / alpha**2 ) * np.outer(nvec/length(nvec), nvec/length(nvec))

# 		answer += 2.0 / ( np.pi * length(nvec)**2 ) * np.exp( - np.pi**2 * length(nvec)**2 / alpha**2 ) * mult

# 	return answer - 3.0 * alpha * Rh / ( 2.0 * np.sqrt(np.pi) * L ) * np.identity(3)

# #-------------------------------------------------------------------------------

# def Oij_pbc_smith( sigma, alpha, m, n ):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		answer += erfc( alpha * length( mvec + sigma ) ) * O( mvec + sigma )

# 		answer += 2.0 * alpha / np.sqrt(np.pi) * np.exp( - alpha**2 * length( mvec + sigma )**2 ) * np.outer(mvec+sigma, mvec+sigma)/length(mvec+sigma)**2

# 	for nvec in ns:

# 		mult = 2.0 / ( np.pi * length(nvec)**2 ) * np.exp( - np.pi**2 * length(nvec)**2 / alpha**2 ) * np.exp( 2 * np.pi * 1j * np.dot(nvec, sigma)  )

# 		mult_real = mult.real

# 		answer += mult_real * ( np.identity(3) - (1 + np.pi**2 * length(nvec)**2 / alpha**2 ) * np.outer(nvec / length(nvec), nvec / length(nvec)) )

# 	return answer

# #-------------------------------------------------------------------------------

# def Qii_pbc_smith( Rh, L, alpha, m, n ):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) if not (mi==0 and mj==0 and mk==0 ) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		mult = erfc( alpha * length(mvec) ) + 2.0 * alpha / np.sqrt(np.pi) * length(mvec) * np.exp( -alpha**2 * length(mvec)**2 )

# 		answer += mult * Q(mvec)

# 		answer -= 4.0 * alpha**3 / np.sqrt(np.pi) * np.exp(-alpha**2 * length(mvec)**2) * np.outer(mvec/length(mvec), mvec/length(mvec))

# 	for nvec in ns:

# 		answer += 4.0 * np.pi * np.exp( -np.pi**2 * length(nvec)**2 / alpha**2) * np.outer(nvec/length(nvec), nvec/length(nvec))

# 	return answer - 1.0 / ( 3.0 * np.sqrt(np.pi) ) * (alpha * Rh / L)**3 * np.identity(3)

# #-------------------------------------------------------------------------------

# def Qij_pbc_smith( sigma, alpha, m, n ):

# 	ms = [ np.array([ mi,mj,mk ], float) for mi in range(m+1) for mj in range(m+1) for mk in range(m+1) if (mi+mj+mk<=m) ]

# 	ns = [ np.array([ ni,nj,nk ], float) for ni in range(n+1) for nj in range(n+1) for nk in range(n+1) if (ni+nj+nk<=n) if not (ni==0 and nj==0 and nk==0 ) ]

# 	answer = 0.0

# 	for mvec in ms:

# 		mult = erfc( alpha * length( mvec + sigma ) ) + 2.0 * alpha / np.sqrt( np.pi ) * length( mvec + sigma ) * np.exp(- alpha**2 * length(mvec + sigma)**2)

# 		answer += mult * Q(mvec + sigma)

# 		answer -= 4.0 * alpha**3 / np.sqrt(np.pi) * np.exp( -alpha**2 * length(mvec+sigma)**2 ) * np.outer(mvec+sigma,mvec+sigma)/length(mvec+sigma)**2

# 	for nvec in ns:

# 		addi = 4.0 * np.pi * np.exp( -np.pi**2 * length(nvec)**2 / alpha**2 ) * np.exp(2*np.pi*1j*np.dot(nvec,sigma)) * np.outer(nvec/length(nvec), nvec/length(nvec))

# 		addi_real = addi.real

# 		answer += addi_real

# 	return answer

# #-------------------------------------------------------------------------------

# def Mii_rpy_smith(p, L, alpha, m, n):

# 	coef1 = 1.0 / ( 6 * np.pi * p.Rh * VISCOSITY )

# 	coef2 = 3.0 * p.Rh / ( 4.0 * L )

# 	coef3 = ( p.Rh / L )**3 / 2.0

# 	comp1 = np.identity(3)

# 	comp2 = Oii_pbc_smith( p.Rh, L, alpha, m, n )

# 	comp3 = Qii_pbc_smith( p.Rh, L, alpha, m, n )

# 	return coef1 * ( comp1 + coef2 * comp2 + coef3 * comp3 )

# #-------------------------------------------------------------------------------

# def Mij_rpy_smith(pi, pj, L, alpha, m, n):

# 	sigma = np.array( (pj.r - pi.r) / L, float )

# 	coef1 = 1.0 / ( 6 * np.pi * pi.Rh * VISCOSITY )

# 	coef2 = 3.0 * pi.Rh / ( 4.0 * L )

# 	coef3 = ( pi.Rh / L )**3 / 2.0

# 	comp1 = Oij_pbc_smith( sigma, alpha, m, n )

# 	comp2 = Qij_pbc_smith( sigma, alpha, m, n )

# 	return coef1 * ( coef2 * comp1 + coef3 * comp2 )

# #-------------------------------------------------------------------------------

# def Bij_rpy_smith(pi, pj):

# 	rij = pj.r - pi.r

# 	coef1 = 1.0 / ( 8 * np.pi * length(rij) * VISCOSITY )

# 	coef2 = 1.0 + ( pi.Rh**2 + pj.Rh**2 ) / ( 3.0 * length(rij)**2 )

# 	coef3 = 1.0 - ( pi.Rh**2 + pj.Rh**2 ) / length(rij)**2

# 	return coef1 * ( coef2 * np.identity(3) + coef3 * np.outer(rij, rij) / length(rij)**2 )

# #-------------------------------------------------------------------------------

# def M_rpy_smith(ps, L, alpha, m, n):

# 	M = [ [ None for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		M[i][i] = Mii_rpy_smith(pi, L, alpha, m, n)

# 		for j in range(i):

# 			pj = ps[j]

# 			M[i][j] = Mij_rpy_smith(pi, pj, L, alpha, m, n)
# 			M[j][i] = np.transpose( M[i][j] )

# 	return np.block(M)

# #-------------------------------------------------------------------------------

# def B_rpy_smith(ps):

# 	B = [ [ None for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		B[i][i] = np.zeros((3,3))

# 		for j in range(i):

# 			pj = ps[j]

# 			B[i][j] = Bij_rpy_smith(pi, pj)
# 			B[j][i] = np.transpose( B[i][j] )

# 	return np.block(B)

# #-------------------------------------------------------------------------------
# # LUBRICATION
# #-------------------------------------------------------------------------------

# def X_f_poly(l, rank):

# 	if rank == 0: return 1
# 	if rank == 1: return 3 * l
# 	if rank == 2: return 9 * l
# 	if rank == 3: return -4 * l + 27 * l**2 - 4 * l**3
# 	if rank == 4: return -24 * l + 81 * l**2 + 36 * l**3
# 	if rank == 5: return 72 * l**2 + 243 * l**3 + 72 * l**4
# 	if rank == 6: return 16 * l + 108 * l**2 + 281 * l**3 + 648 * l**4 + 144 * l**5
# 	if rank == 7: return 288 * l**2 + 1620 * l**3 + 1515 * l**4 + 1620 * l**5 + 288 * l**6
# 	if rank == 8: return 576 * l**2 + 4848 * l**3 + 5409 * l**4 + 4524 * l**5 + 3888 * l**6 + 576 * l**7
# 	if rank == 9: return 1152 * l**2 + 9072 * l**3 + 14752 * l**4 + 26163 * l**5 + 14752 * l**6 + 9072 * l**7 + 1152 * l**8
# 	if rank == 10: return 2304 * l**2 + 20736 * l**3 + 42804 * l**4 + 115849 * l**5 + 76176 * l**6 + 39264 * l**7 + 20736 * l**8 + 2304 * l**9
# 	if rank == 11: return 4608 * l**2 + 46656 * l**3 + 108912 * l**4 + 269100 * l**5 + 319899 * l**6 + 269100 * l**7 + 108912 * l**8 + 46656 * l**9 + 4608 * l**10
# 	else: return None

# def Y_f_poly(l, rank):

# 	if rank == 0: return 1
# 	if rank == 1: return 3 / 2 * l
# 	if rank == 2: return 9 / 4 * l
# 	if rank == 3: return 2 * l + 27 / 8 * l**2 + 2 * l**3
# 	if rank == 4: return 6 * l + 81 / 16 * l**2 + 18 * l**3
# 	if rank == 5: return 63 / 2 * l**2 + 243 / 32 * l**3 + 63 / 2 * l**4
# 	if rank == 6: return 4 * l + 54 * l**2 + 1241 / 64 * l**3 + 81 * l**4 + 72 * l**5
# 	if rank == 7: return 144 * l**2 + 1053 / 8 * l**3 + 19083 / 128 * l**4 + 1053 / 8 * l**5 + 144 * l**6
# 	if rank == 8: return 279 * l**2 + 4261 / 8 * l**3 + 126369 / 256 * l**4 - 117 / 8 * l**5 + 648 * l**6 + 288 * l**7
# 	if rank == 9: return 576 * l**2 + 1134 * l**3 + 60443 / 32 * l**4 + 766179 / 512 * l**5 + 60443 / 32 * l**6 + 1134 * l**7 + 576 * l**8
# 	if rank == 10: return 1152 * l**2 + 7857 / 4 * l**3 + 98487 / 16 * l**4 + 10548393 / 1024 * l**5 + 67617 / 8 * l**6 - 351 / 2 * l**7 + 3888 * l**8 + 1152 * l**9
# 	if rank == 11: return 2304 * l**2 + 7128 * l**3 + 22071 / 2 * l**4 + 2744505 / 128 * l**5 + 95203835 / 2048 * l**6 + 2744505 / 128 * l**7 + 22071 / 2 * l**8 + 7128 * l**9 + 2304 * l**10
# 	else: return None

# #-------------------------------------------------------------------------------

# def X_g_poly(l, rank):

# 	if rank == 1: return 2 * l**2 * ( 1 + l )**(-3)
# 	if rank == 2: return 1 / 5 * l * ( 1 + 7 * l + l**2 ) * ( 1 + l )**(-3)
# 	if rank == 3: return 1 / 42 * ( 1 + 18 * l - 29 * l**2 + 18 * l**3 + l**4 ) * ( 1 + l )**(-3)
# 	else: return None

# #-------------------------------------------------------------------------------

# def Y_g_poly(l, rank):

# 	if rank == 2: return 4 / 15 * l * ( 2 + l + 2 * l**2 ) * ( 1 + l )**(-3)
# 	if rank == 3: return 2 / 375 * ( 16 - 45 * l + 58 * l**2 - 45 * l**3 + 16 * l**4 ) * ( 1 + l )**(-3)
# 	else: return None

# #-------------------------------------------------------------------------------

# def XA11(s, l):

# 	answer = 0.0

# 	answer += X_g_poly(l, 1) * ( 1 - 4 * s**(-2) )**(-1)

# 	answer -= X_g_poly(l, 2) * np.log( 1 - 4 * s**(-2) )

# 	answer -= X_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( 1 - 4 * s**(-2) )

# 	answer += X_f_poly(l, 0) - X_g_poly(l, 1)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 0 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * X_f_poly(l, m) - X_g_poly(l, 1) )

# 		answer += mult * ( 4 * m**(-1) * m1**(-1) * X_g_poly(l, 3) - 2 * m**(-1) * X_g_poly(l, 2) )

# 	return answer

# #-------------------------------------------------------------------------------

# def YA11(s, l):

# 	answer = 0.0

# 	answer -= Y_g_poly(l, 2) * np.log( 1 - 4 * s**(-2) )

# 	answer -= Y_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( 1 - 4 * s**(-2) )

# 	answer += Y_f_poly(l, 0)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 0 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * Y_f_poly(l, m) - 2 * m**(-1) * Y_g_poly(l, 2) )

# 		answer += mult * 4 * m**(-1) * m1**(-1) * Y_g_poly(l, 3)

# 	return answer

# #-------------------------------------------------------------------------------

# def XA12(s, l):

# 	answer = 0.0

# 	answer += 2 * s**(-1) * X_g_poly(l, 1) * ( 1 - 4 * s**(-2) )**(-1)

# 	answer += X_g_poly(l, 2) * np.log( ( s + 2 ) / ( s - 2 ) )

# 	answer += X_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( ( s + 2 ) / ( s - 2 ) ) + 4 * X_g_poly(l, 3) * s**(-1)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 1 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * X_f_poly(l, m) - X_g_poly(l, 1) )

# 		answer += mult * ( 4 * m**(-1) * m1**(-1) * X_g_poly(l, 3) - 2 * m**(-1) * X_g_poly(l, 2) )

# 	divisor = -1 / 2 * ( 1 + l )

# 	return answer / divisor

# #-------------------------------------------------------------------------------

# def YA12(s, l):

# 	answer = 0.0

# 	answer += Y_g_poly(l, 2) * np.log( ( s + 2 ) / ( s - 2 ) )

# 	answer += Y_g_poly(l, 3) * ( 1 - 4 * s**(-2) ) * np.log( ( s + 2 ) / ( s - 2 ) )

# 	answer += 4 * Y_g_poly(l, 3) * s**(-1)

# 	for m in [ mi for mi in range(1, 12) if ( mi%2 == 1 ) ]:

# 		if m == 2: m1 = -2
# 		else: m1 = m - 2

# 		mult = ( 2 / s )**m

# 		answer += mult * ( 2**(-m) * ( 1 + l )**(-m) * Y_f_poly(l, m) - 2 * m**(-1) * Y_g_poly(l, 2) )

# 		answer += mult * ( 4 * m**(-1) * m1**(-1) * Y_g_poly(l, 3) )

# 	divisor = -1 / 2 * ( 1 + l )

# 	return answer / divisor

# #-------------------------------------------------------------------------------

# def R_jeffrey(p1, p2):

# 	s = 2 * p1.distance(p2) / ( p1.Rh + p2.Rh )

# 	l = p2.Rh / p1.Rh

# 	r = p2.r - p1.r

# 	R = [ [ None , None ], [ None, None ] ]

# 	R[0][0] = XA11(s, l) * np.outer(r/length(r), r/length(r)) + YA11(s, l) * ( np.identity(3) - np.outer(r/length(r), r/length(r)) )

# 	R[1][1] = XA11(s, 1/l) * np.outer(r/length(r), r/length(r)) + YA11(s, 1/l) * ( np.identity(3) - np.outer(r/length(r), r/length(r)) )

# 	R[0][1] = XA12(s, l) * np.outer(r/length(r), r/length(r)) + YA12(s, l) * ( np.identity(3) - np.outer(r/length(r), r/length(r)) )

# 	R[1][0] = XA12(s, 1/l) * np.outer(r/length(r), r/length(r)) + YA12(s, 1/l) * ( np.identity(3) - np.outer(r/length(r), r/length(r)) )

# 	return 3 * np.pi * ( p1.Rh + p2.Rh ) * VISCOSITY * np.block(R)

# #-------------------------------------------------------------------------------

# def R_lub_corr(ps):

# 	corr = [ [ np.zeros((3,3)) for j in range( len(ps) ) ] for i in range( len(ps) ) ]

# 	for i, pi in enumerate(ps):

# 		for j in range(i):

# 			pj = ps[j]

# 			nf2b = R_jeffrey( pi, pj )

# 			ff2b = np.linalg.inv( M_rpy( [pi, pj] ) )

# 			lub_corr = nf2b - ff2b

# 			corrii = lub_corr[0:3,0:3]
# 			corrjj = lub_corr[3:6,3:6]
# 			corrij = lub_corr[3:6,0:3]

# 			corr[i][i] += corrii
# 			corr[j][j] += corrjj
# 			corr[i][j] += corrij
# 			corr[j][i] += np.transpose( corrij )

# 	return np.block(corr)

# #-------------------------------------------------------------------------------
# # END LUBRICATION
# #-------------------------------------------------------------------------------

# def corr(M):

# 	diagonal = np.diag( M )
# 	diagonal_only = np.diag(1 / np.sqrt(diagonal))

# 	result = np.dot(diagonal_only, M)

# 	return np.dot(result, np.transpose(diagonal_only))

# #-------------------------------------------------------------------------------

# def main():

# 	p1 = Particle([0.0,0.0,0.0], 1.0)
# 	p2 = Particle([3.0,0.0,0.0], 1.0)
# 	p3 = Particle([6.0,0.0,0.0], 1.0)
# 	# p4 = Particle([3.0,0.0,0.0], 1.0)
# 	# p5 = Particle([0.0,0.0,0.0], 1.0)
# 	# p6 = Particle([3.0,0.0,0.0], 1.0)
# 	distance = 3.0
# 	ps = [ Particle([i*3.0,0.0,0.0], 1.0) for i in range(100) ]
	
# 	rpy = np.linalg.inv( M_rpy( ps ) )
# 	# smith = M_rpy_smith( [p1, p2], L = 1000.0, alpha = np.sqrt( np.pi ), m = 1, n = 1 )
# 	# jeffrey = R_jeffrey( *[p1, p2] )
# 	lub_corr = R_lub_corr( ps ) + rpy

# 	print( 'R RPY:\n{}\n'.format( rpy ) )
# 	# print( 'R Smith:\n{}\n'.format( np.linalg.inv(smith) ) )
# 	# print( 'R Jeffrey:\n{}\n'.format( jeffrey ) )
# 	# print( 'R Jeffrey / R RPYL\n{}\n'.format( jeffrey / np.linalg.inv(rpy) ) )
# 	print( 'R lubrication corrected:\n{}\n'.format( lub_corr ) )
# 	# print( 'Relative change:\n{}\n'.format( lub_corr / rpy ) )
# 	# Dlubcorr = KT * np.linalg.inv( np.identity(len(smith)) + smith @ R_lub_corr( [p1, p2] ) ) @ smith
# 	# Dsmith = KT * smith
# 	# print( 'D lub corr:\n{}\n'.format( Dlubcorr ) )
# 	# print( 'D Smith:\n{}\n'.format( Dsmith ) )
# 	# print( 'D lub corr / D Smith:\n{}\n'.format( Dlubcorr/Dsmith ) )
	
# 	# import matplotlib.pyplot as plt
# 	# import matplotlib.animation as animation
# 	# from matplotlib.patches import Circle
	
# 	# steps = 1
# 	# dt = 0.0001
# 	# force = np.array( [0.0, 0.0, 0.0, -0.0, 0.0, 0.0] )
# 	# pointsize = 100
	
# 	# fig,axs = plt.subplots(1,3)
# 	# ax1,ax2,ax3 = axs
# 	# ax1.set_aspect('equal')
# 	# ax1.set_xlim( (0.0, 3.0) )
# 	# ax1.set_ylim( (-1.5, 1.5) )
# 	# ax2.set_aspect('equal')
# 	# ax2.set_xlim( (0.0, 3.0) )
# 	# ax2.set_ylim( (-1.5, 1.5) )
# 	# ax3.set_aspect('equal')
# 	# ax3.set_xlim( (0.0, 3.0) )
# 	# ax3.set_ylim( (-1.5, 1.5) )
# 	# point1, = ax1.plot([], [], 'o', ms=pointsize, color = 'blue')
# 	# point2, = ax1.plot([], [], 'o', ms=pointsize, color = 'red')
# 	# point3, = ax2.plot([], [], 'o', ms=pointsize, color = 'blue')
# 	# point4, = ax2.plot([], [], 'o', ms=pointsize, color = 'red')
# 	# point5, = ax3.plot([], [], 'o', ms=pointsize, color = 'blue')
# 	# point6, = ax3.plot([], [], 'o', ms=pointsize, color = 'red')
# 	# circle1 = Circle((0,0), 1, color = 'green')
# 	# ax3.add_patch(circle1)
	
# 	# def init():
# 	# 	"""initialize animation"""
# 	# 	point1.set_data([], [])
# 	# 	point2.set_data([], [])
# 	# 	point3.set_data([], [])
# 	# 	point4.set_data([], [])
# 	# 	point5.set_data([], [])
# 	# 	point6.set_data([], [])
# 	# 	circle1 = Circle( (0.0,0.0), 1.0, color = 'green' )
# 	# 	return point1,point2,point3,point4,point5,point6,circle1
	
# 	# def animate(i):
# 	# 	"""perform animation step"""
# 	# 	global dt
	
# 	# 	rpy = M_rpy( [p1, p2] )
# 	# 	coords1 = np.array( [ x for x in p1.r ] + [ x for x in p2.r ] )
# 	# 	coords1 += dt * rpy @ force
# 	# 	p1.r = np.array( coords1[0:3] )
# 	# 	p2.r = np.array( coords1[3:6] )
	
# 	# 	smith = M_rpy_smith( [p3, p4], L = 1000.0, alpha = np.sqrt( np.pi ), m = 2, n = 2 )
# 	# 	coords2 = np.array( [ x for x in p3.r ] + [ x for x in p4.r ] )
# 	# 	coords2 += dt * smith @ force
# 	# 	p3.r = np.array( coords2[0:3] )
# 	# 	p4.r = np.array( coords2[3:6] )
	
# 	# 	smith = M_rpy_smith( [p5, p6], L = 1000.0, alpha = np.sqrt( np.pi ), m = 2, n = 2 )
# 	# 	lub = np.linalg.inv( np.identity(len(smith)) + smith @ R_lub_corr( [p5, p6] ) ) @ smith
# 	# 	coords3 = np.array( [ x for x in p5.r ] + [ x for x in p6.r ] )
# 	# 	coords3 += dt * lub @ force
# 	# 	p5.r = np.array( coords3[0:3] )
# 	# 	p6.r = np.array( coords3[3:6] )
	    
# 	# 	point1.set_data([p1.r[0]], [p1.r[1]])
# 	# 	point2.set_data([p2.r[0]], [p2.r[1]])
# 	# 	point3.set_data([p3.r[0]], [p3.r[1]])
# 	# 	point4.set_data([p4.r[0]], [p4.r[1]])
# 	# 	point5.set_data([p5.r[0]], [p5.r[1]])
# 	# 	point6.set_data([p6.r[0]], [p6.r[1]])
	
# 	# 	circle1 = Circle( (0.0,0.0), 1.0, color = 'green' )
	
# 	# 	return point1,point2,point3,point4,point5,point6,circle1
	
# 	# # choose the interval based on dt and the time to animate one step
# 	# from time import time
# 	# t0 = time()
# 	# animate(0)
# 	# t1 = time()
# 	# interval = 100 * dt - (t1 - t0)
	
# 	# ani = animation.FuncAnimation(fig, animate, frames=steps,
# 	#                               interval=interval, blit=True, init_func=init)
	
# 	# plt.show()
	
# 	###
	
# 	# for step in range(steps):
	
# 	# 	smith = M_rpy( [p1, p2] )
# 	# 	coords = np.array( [ x for x in p1.r ] + [ x for x in p2.r ] )
# 	# 	print(coords)
# 	# 	coords += dt * smith @ force
# 	# 	print(coords)
# 	# 	p1.r = np.array( coords[0:3] )
# 	# 	p2.r = np.array( coords[3:6] )
	
	
	
	
	
	
	
# 	# print( KT * ( smith - b ) )
# 	# print( corr( rpy ) )
# 	# print( corr( smith ) )
	
# 	#-------------------------------------------------------------------------------
	
# 	# from sympy import *
# 	# from sympy import Matrix
	
# 	# b, a, a1, a2, r = symbols('b a a1 a2 r')
# 	# A = Matrix( [[b/a1,0,0,3/4*b/r*(2-2*a/3),0,0],
# 	# 			 [0,b/a1,0,0,3/4*b/r*(1+a/3),0],
# 	# 			 [0,0,b/a1,0,0,3/4*b/r*(1+a/3)],
# 	# 			 [3/4*b/r*(2-2*a/3),0,0,b/a2,0,0],
# 	# 			 [0,3/4*b/r*(1+a/3),0,0,b/a2,0],
# 	# 			 [0,0,3/4*b/r*(1+a/3),0,0,b/a2]] ) # Creates a matrix.
# 	# # pprint(A)
# 	# A_inverse = A.inv() #Doesn't work
# 	# # pprint( A_inverse[0:3,0:3] )
# 	# b_val = 1.0 / ( 6.0 * np.pi * VISCOSITY )
# 	# a_val = ( p1.Rh**2 + p2.Rh**2 ) / 3.0**2
# 	# a1_val = p1.Rh
# 	# a2_val = p2.Rh
# 	# result = A_inverse.subs([(b,b_val), (a,a_val), (a1,a1_val), (a2,a2_val), (r,3.0)])
# 	# pprint(result)

# #-------------------------------------------------------------------------------

# if __name__ == '__main__':

# 	main()