#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 13:38:59 2021

@author: elliott
"""

import igrf
import numpy as np

Re = 6357.e3

def norm_theta(theta, thetamax, theta_limit=np.pi):
	
	"""
	PURPOSE:
		This function returns the adjusted values of the 
		angle theta, such that the full range of theta values run
		from 0 to Theta_Limit, where Theta_Limit is either pi or pi/2

	$Log: norm_theta.pro,v $
	Revision 1.1  1997/08/05 14:50:53  baker
	Initial revision
	"""
	
	alpha = theta_limit/thetamax
	theta_prime = alpha * theta
	
	return theta_prime, alpha

def dbang(Mmax):
	
	"""
	FUNCTION EVAL_LEGENDRE

	Purpose:  evalaute all the Associated Legendre Polynomials
				from L=0 to L=Lmax, at a set of points

	Calling Sequence:
		PLM = eval_legendre(Lmax, pos)

		where Lmax give the maximum order of the Legendre
		polynomials, and pos gives the positions where the
		the polynomials are to be evaluated.

		pos = pos(2,N) where the first index indicates
				the theta,phi position and the second index
				lists the points where we have data.


	
	$Log: eval_legendre.pro,v $
	Revision 1.1  1997/08/05 14:42:55  baker
	Initial revision
	"""
	
	result = np.ones(Mmax+1, dtype = np.double)
	for i in range(1, Mmax+1):
		for j in range(i, Mmax+1):
			result[j] = result[j]*(2*i - 1)
			
	return result

def index_legendre(l, m, keyword):
	
	"""
	PURPOSE
	This routine converts a Legendre polynoiomial index pair (l,m)
	into a single index (k).

	FUNCTION:  INDEX_LEGENDRE

	Calling Sequence:
		k = index_legendre(l,m, keyword)

	The keywords SH, SNGL, and DBL have the following
	meanings:
		SH:  We are doing Spherical harmonics where m runs from -l to + l

		SINGLE:  We are doing Associated Legendre Polynomials with m=0,l

		DOUBLE:  We are doing Associated Legendre Polynomials
					but for each polynomial we have two coefficients 
					one for cos(phi) and one for sin(phi), as before, m
					runs from 0 to l.  Basically, /DOUBLE means we
					are doing spherical harmonics for a real valued 
					function using sin(phi) and cos(phi) rather than
					exp(i*phi).

	$Log: index_legendre.pro,v $
	Revision 1.1  1997/08/05 14:47:30  baker
	Initial revision
	"""
	
	if m > l: 
		return -1
	
	if keyword == "SH":
		return int(m + l*(l+1))
	
	elif keyword == "single":
		return int(m + l*(l+1)/2)
	
	elif keyword == "double":
		if l == 0: 
			return 0
		elif m == 0:
			return l*l
		else:
			return l*l + 2*m - 1
	else:
		print("keyword must be SH, single, or double.")
	
	return

def eval_potential(a, plm, phi):

	"""
	PURPOSE:  evaluate the electric potential on a set of
				points, given the coefficients of the spherical
				harmonic expansion.

	Calling Sequence:

	pot = eval_potential,a,plm,phi

	where 'a' is the set of coefficients given as a vector
				indexed by k = index_legendre(l,m,/dbl)

	plm is an array (N,Lmax,Lmax) where N = number of points
		where the potential is to be evaluated, and the
		other two indices give the associated Legendre polynomial
		to be evaluated at each position.

	phi is the azimuthal coordinate at each evaluation point (in radians).
	"""
	
	#dictionary of idl type codes
	type_codes = {"byte":1, "int":2, "long":3, "float":4, "double":5, 
			   "complex":6, "string":7, "struct":8}
	
	#get code for our type
	lmax = type_codes[str(type(plm[0]))[8:-2]]-1
	v = np.zeros(len(phi))
	for m in range(0, lmax+1):
		for l in range(m, lmax+1):
			k = index_legendre(l, m, "double")
			if m == 0:
				v += a[k]*plm[:,l,0]
			else:
				v += a[k]*np.cos(m*phi)*plm[:,l,m] + \
					a[k+1]*np.sin(m*phi)*plm[:,l,m]
					
	return v

def eval_legendre(lmax, x_in, keyword="n/a"):
	
	x = np.atleast_2d(np.double(x_in.reshape(len(x_in)))).T
	N = len(x)
	xx = (1-(x**2))
	extension = np.atleast_2d(np.ones(lmax+1))
	xx_L = np.matmul(xx, extension, dtype=np.double)
	mover2 = np.atleast_2d(np.arange(lmax+1)/2.0)
	extension2 = np.atleast_2d(np.ones(N, dtype = np.double)).T
	mover2 = np.matmul(extension2, mover2)
	xx_Mover2 = xx_L**mover2
	
	# xx_Mover2 is the matrix of all the (1-x^2) values raised to all the
	# powers of m/2 with m running from 0 to lmax
	
	two_m_dbang = np.atleast_2d(dbang(lmax))
	two_m_dbang = np.matmul(extension2, two_m_dbang)
	extension3 = np.atleast_2d(pow(-1, (np.arange(lmax+1))))
	pwrm = np.multiply(extension2, extension3)
	pmm = xx_Mover2*pwrm*two_m_dbang
	
	# we have now computed the value of Pmm at each point for each value
	# of m from 0 to lmax
	
	# we have now computed P(m+1, m) at each point and for each value of m
	# from 0 to lmax
	
	# p(m+1, m) = x(2m+1)P(m, m)
	
	extension4 = np.atleast_2d(np.ones(lmax+1, dtype = np.double))
	
	pmmp1 = np.matmul(x, extension4) * \
		np.matmul(extension2, np.atleast_2d(np.arange(lmax+1.)*2.0+1.0)) * \
			pmm
			
	# Ok, now we have pmm and pmmp1 at every point. Next we have to compute
	# the rest of the plm values from the recursion relation;
	# (l-m)P(l, m) = x(2l-1)p(l-1,m) - (l+m-1)P(l-2,m)
	
	plm = np.zeros([N, lmax+1, lmax+1], dtype = np.double)
	for l in range(0, lmax+1):
		plm[:,l,l] = pmm[:,l]
	for l in range(0, lmax):
		plm[:,l+1,l] = pmmp1[:,l]
	for l in range(0, lmax-1):
		for k in range(l+2, lmax+1):
			plm[:,k,l] = (1.0/(k-l))*((x[:,0]*plm[:,k-1,l]*(2*k -1))-(k+l-1)*plm[:,k-2,l])
			
	if keyword == "pm":
		pm = pmm
	if keyword == "pp":
		pp = pmmp1
		
	return plm

def eval_etheta_coeffs(Pcoeffs, theta, lmax, latmin):
	
	"""
	PURPOSE
	This set of routines is used to evaluate the Electric field
	given a set of coefficients defining the potential.
	"""

	theta_max = np.deg2rad((90-latmin))
	alpha = 1.
	theta_prime, alpha = norm_theta(theta, theta_max)
	
	n = len(theta)
	
	kmax = index_legendre(lmax, lmax, "double")
	ecoeffs = np.zeros([kmax+2, n], float)
	q = np.where(theta_prime != 0.0)
	
	for m in range(0, lmax+1):
		for l in range(m, lmax+1):
			k3 = index_legendre(l, m, "double")
			k4 = index_legendre(l, m, "double")
			if k3 > 0:
				ecoeffs[k4, q] = ecoeffs[k4,q] - (Pcoeffs[k3]*alpha*l*np.cos(theta_prime[q])/np.sin(theta_prime[q]))/Re
			if l < lmax:
				k1 = index_legendre(l+1, m, "double")
			else:
				k1 = -1	
			k2 = index_legendre(l, m, "double")

			if k1 > 0:
				ecoeffs[k2, q] = ecoeffs[k2, q] + (Pcoeffs[k1]*alpha*(l+1+m) / np.sin(theta_prime[q]))/Re
					
			if m > 0:
				if k3 > 0: k3 += 1
				k4 += 1
				if k1 > 0: k1 += 1
				k2 += 1
				if k3 > 0:
					ecoeffs[k4, q] = ecoeffs[k4, q] - (Pcoeffs[k3]*alpha*l*np.cos(theta_prime[q]) / np.sin(theta_prime[q]))/Re
				if k1 > 0:
					ecoeffs[k2, q] = ecoeffs[k2, q] + (Pcoeffs[k1]*alpha*(l+1+m) / np.sin(theta_prime[q]))/Re
						
	return ecoeffs
					
def eval_ephi_coeffs(Pcoeffs, theta, lmax, latmin):
	
	n = len(theta)
	
	kmax = index_legendre(lmax, lmax, "double")
	ecoeffs = np.zeros([kmax+2, n], float)
	q = np.where(theta != 0)
	
	for m in range(1, lmax+1):
		for l in range(m, lmax+1):
			k3 = index_legendre(l, m, "double")
			k4 = index_legendre(l, m, "double")
			if k3 > 0:
				ecoeffs[k4, q] = ecoeffs[k4, q] - Pcoeffs[k3+1]*m/np.sin(theta[q])/Re
				ecoeffs[k4+1, q] = ecoeffs[k4+1, q] + Pcoeffs[k3]*m/np.sin(theta[q])/Re
				
	return ecoeffs

def eval_component(ecoeffs, plm, pos, lmax):
	
	theta = np.deg2rad(90.0-pos[0,])
	phi = np.deg2rad(pos[1,])
	n = len(theta)
	ecomp = np.zeros(n, float)
	
	for m in range(0, lmax+1):
		for l in range(m, lmax+1):
			k = index_legendre(l, m, "double")
			if m == 0:
				ecomp += ecoeffs[k,]*plm[:,l,m]
			else:
				ecomp += ecoeffs[k,]*plm[:,l,m]*np.cos(m*phi) + \
					ecoeffs[k+1,]*plm[:,l,m]*np.sin(m*phi)
					
	return ecomp

def eval_efield(pcoeffs, plm, pos, lmax, latmin):
	
	theta = pos[0,]
	theta = np.deg2rad(90.0-theta)
	
	etc = eval_etheta_coeffs(pcoeffs, theta, lmax, latmin)
	epc = eval_ephi_coeffs(pcoeffs, theta, lmax, latmin)

	etheta = eval_component(etc, plm, pos, lmax)
	ephi = eval_component(epc, plm, pos, lmax)
	
	n = len(theta)
	e_field = np.zeros([2, n], float)
	e_field[0,] = etheta
	e_field[1,] = ephi
	
	return e_field

def cal_efield(pos, solution, latmin, lon_shft, lat_shft, order):
	
	theta = np.deg2rad(90. - pos[0,])
	thetamax = np.deg2rad(90. - latmin)
	theta_prime, alpha = norm_theta(theta, thetamax)
	x = np.cos(theta_prime)
	plm = eval_legendre(order, x)
	lmax = order
	pcoeffs=solution[2,]
	
	etc = eval_etheta_coeffs(pcoeffs, theta, lmax, latmin)
	epc = eval_ephi_coeffs(pcoeffs, theta, lmax, latmin)
	
	etheta = eval_component(etc, plm, pos, lmax)
	ephi = eval_component(epc, plm, pos, lmax)
	
	n = len(theta)
	
	e_field = np.zeros([2, n], float)
	e_field[0,] = etheta
	e_field[1,] = ephi
	
	return e_field
	
def eval_vel(pcoeffs, plm, pos, lmax, latmin, dtime):
	
	e_field = eval_efield(pcoeffs, plm, pos, lmax, latmin)
	
	Altitude = 300*1000
	bpolar = -.62e-4
	phi = 90 - pos[0,]
	#bmag = bpolar*(1.-3.*Altitude/Re)*np.sqrt(3.*np.cos(np.deg2rad(phi))**2 +1.)/2.
	bmag = np.empty(len(e_field[0,]))
	#get magnetic field strength from IGRF model
	for i in range(len(bmag)):
		bfield = igrf.igrf(dtime, pos[0,i], pos[1,i], 150)
		bmag[i] = bfield.total.to_dict()["data"][0]*-1e-9
	vel = np.copy(e_field)
	vel[0,] = e_field[1,]/bmag
	vel[1,] = -e_field[0,]/bmag
	
	return vel

def calc_vels(pos, solution, latmin, order, dtime):
	
	theta = np.deg2rad(90. - pos[0,])
	thetamax = np.deg2rad(90.-latmin)
	theta_prime, alpha = norm_theta(theta, thetamax)
	x = np.cos(theta_prime)
	plm = eval_legendre(order, x)
	vvec = eval_vel(solution[2,], plm, pos, order, latmin, dtime)
	
	return vvec

def find_gradV(pos, solution, latmin, lon_shft, lat_shft, order, dtime):
	
	# First shift coordinates into "model" reference (pole shifted 4 deg nightwards)
	posx = pos # lat, lon
	if lat_shft != 0:
		npnts = len(posx)/2
		kaz = np.zeros(npnts, float)

	# calculate vectors
	vvec = calc_vels(posx, solution, latmin, order, dtime)
	vmag = np.sqrt(vvec[0,]**2 + vvec[1,]**2)
	if len(vmag) == 1:
		vmag = vmag[0]

	q = np.where(vmag != 0)
	qc = len(q[0])
	
	if qc == 0:
		print("all vectors have 0 length")
		return
	
	#print("vmag", vmag)
	if isinstance(vmag, np.ndarray):
		vaz = np.zeros(len(vmag), float)
		vaz[q] = np.rad2deg(np.arctan2(vvec[1,q], -vvec[0,q]))
	else:
		vaz = np.rad2deg(np.arctan2(vvec[1,q], - vvec[0,q]))
		while isinstance(vaz, np.ndarray):
			vaz = vaz[0]
			
	
	# Now shift back into "real world"
	if lat_shft != 0:
		xat_shft = -lat_shft
		npnts = len(vmag)
	
	return vmag, vaz