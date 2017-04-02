import numpy as np

import constants
import utils

class Satellite():
	
	def __init__(self, hp, ha, name=None):
		
		self.hp = hp
		self.ha = ha
		self.Rp = hp * 1e3 + constants.radiusE
		self.Ra = ha * 1e3 + constants.radiusE
		self.a = (self.Rp + self.Ra)/2.#800#(self.hp + self.ha) * 0.5# * 1e3
		self.name = name
		
		self.P = 2. * np.pi * np.sqrt((self.a * self.a * self.a) / constants.muE)
		self.orbital_speed = np.sqrt(constants.muE / self.a)
		self.e = self.Ra / self.a - 1.
		
		self.x = None
		self.y = None
		self.u = None
		self.v = None
		
		self.reset_mem()
		
	def reset_mem(self):
		
		self.trajectory_x = []
		self.trajectory_y = []
		self.speed_x = []
		self.speed_y = []
		
	def save_state(self, x=None, y=None, u=None, v=None):
		
		if x is None and y is None:
			x = self.x 
			y = self.y
		if u is None and v is None:
			u = self.u 
			v = self.v
		
		self.trajectory_x.append(x)
		self.trajectory_y.append(y)
		self.speed_x.append(u)
		self.speed_y.append(v)
		
	def set_init_pos(self, true_anomaly):
		
		alpha = np.deg2rad(true_anomaly)
		
		self.init_true_anomaly = alpha
		
		r = self.a * (1. - self.e * self.e) / (1. + self.e * np.cos(alpha))
		
		self.x = np.cos(alpha) * r
		self.y = np.sin(alpha) * r
		
		self.V = np.sqrt(2. * constants.muE / r - constants.muE / self.a)
		
		self.gamma = np.arctan2(self.e * np.sin(alpha), (1. + self.e * np.cos(alpha)))
		beta = alpha - ( np.pi/2. + self.gamma) 
		
		self.u = self.V * np.cos(beta)
		self.v = self.V * np.sin(beta)
		
		self.save_state()

	def step(self, dt, save_state=True):
		
		# Very very simple integration
		# TODO: Do much better, too much deviation!
		d = np.hypot(self.x, self.y)
		a = -constants.muE / d / d / d
		
		#du = utils.RK4(lambda x, u: a * x)
		#dv = utils.RK4(lambda y, v: a * y)

		#du = utils.RK4(lambda x, u: self.x)
		#dv = utils.RK4(lambda x, u: self.y)
		
		#self.u-= du(a, self.u, dt)
		#self.v -= dv(a, self.v, dt)
		#ax = -F * self.x 
		#ay = -F * self.y
		
		self.u += a * self.x * dt
		self.v += a * self.y * dt
		
		self.x += self.u * dt 
		self.y += self.v * dt
		
		if save_state: 
			self.save_state()
		
	def get_alt(self):
		
		return (np.hypot(self.x, self.y) - constants.radiusE) * 1e-3
	
	def get_speed(self):
		
		return np.hypot(self.u, self.v) * 1e-3
	
	def get_posangle(self):
		
		return np.arctan2(self.y, self.x)
	
	def get_true_anomaly(self):
		
		r = np.hypot(self.x, self.y)
		if self.e > 0:
			costheta = self.Rp * (1. + self.e) / (r * self.e) - 1. / self.e
			costheta = np.clip(costheta, -1., 1.)
			return np.arccos(costheta)
		else:
			return self.get_posangle()
		
		