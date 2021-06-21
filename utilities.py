# Copyright (C) 2020-2021 by Andrew Hoffman <hoffmaao@uw.edu>
#
# This file is part of gravity-inversion code.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The full text of the license can be found in the file LICENSE in the
# gravity-inversion source directory or at <http://www.gnu.org/licenses/>.


def msplit(m,**kwargs):
	r"""
	this function is used to split an earth model

	"""
	if len(kwargs)>1:
		dlevel = varargin{1}
	else:
		dlevel = 0

	m1 = np.ones(m.shape[0],m.shape[1])*1e-10
	m2 = m
	for i in range(m.shape[1]):
		m1[:,i] = m1[:,i]*(m.shape[1]+1-i)
		id1 = m[:,i]<0
		m1[id1,i] = m[id1,i]
		m2[id1,i] = i*1e-10
	m1 = m1+dlevel
	m2 = m2+dlevel
	return m1, m2

def twodpolyA(xobs,M,rho0,rhoback):
	"""
	TWODPOLYA computes gravity anomaly Dg of polygons using Talwani line 
	integral. For details, see Talwani et al. (1959), JGR, 64(1), 49-59
	"""
	N = len(rho0)
	mx = len(xobs)
	drho = lrho0-rhoback
	dg = np.empty((mx,N))
	dg[:]=np.nan
	for i in range(N):
		id = #find function
		xp = M[:,1]
		xp[i+1] = []

		zp = M[:,i+1]
		zp[id+1] = []
		zp[zp==0] = 1e-10
		dg[:,i] = talwani2(xobs,xp,zp,drho[ip])

	return np.sum(dg,2)*1e5



def mmodtwodpoly(m0,**kwargs):
	r"""
	This function modifies an earth-model profily into polygons that are exteded 
	+/- 50 km of either end.

	"""
	a = np.array([m0[0,0]-50000,m0;m0[-1,0]+50000, ])
	m = np.array(a[:,0];np.flipud(a[:,0]);a[0,0]);

	S=m.shape()

	for i in range(1,S[1]-1):
		b = np.array(a[:,i];np.flipud(a[:,i+1]);a[1,i])
		m = np.array(m b);

	return m


