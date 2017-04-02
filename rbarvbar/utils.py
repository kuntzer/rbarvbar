import numpy as np

def unit_vector(vector):
	""" Returns the unit vector of the vector.  """
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
	""" Returns the angle in radians between vectors 'v1' and 'v2'::

			>>> angle_between((1, 0, 0), (0, 1, 0))
			1.5707963267948966
			>>> angle_between((1, 0, 0), (1, 0, 0))
			0.0
			>>> angle_between((1, 0, 0), (-1, 0, 0))
			3.141592653589793
	"""
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	angle = np.arccos(np.dot(v1_u, v2_u))
	if np.isnan(angle):
		if (v1_u == v2_u).all():
			return 0.0
		else:
			return np.pi
	return angle

def vangle_between(v1s, v2s):
	return np.array([angle_between(v1, v2) for v1, v2 in zip(v1s, v2s)])

def signed_angle_between(v1, v2, Vn):
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	angle = angle_between(v1_u, v2_u)
	
	cross = np.cross(v1_u, v2_u)
	if (np.dot(Vn, cross) < 0):
		angle = -angle
	
	return angle

def cart2polar(x, y):
	theta = np.arctan2(y, x)
	r = np.hypot(x, y)
	
	return r, theta 

def polar2cart(r, theta):
	x = np.cos(theta) * r
	y = np.sin(theta) * r
	
	return x, y

def RK4(f):
	return lambda t, y, dt: (
			lambda dy1: (
			lambda dy2: (
			lambda dy3: (
			lambda dy4: (dy1 + 2*dy2 + 2*dy3 + dy4)/6
			)( dt * f( t + dt  , y + dy3   ) )
		)( dt * f(t + dt/2, y + dy2/2 ) )
		)( dt * f(t + dt/2, y + dy1/2 ) )
		)( dt * f(t, y) )
