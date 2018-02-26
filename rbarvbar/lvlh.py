import numpy as np
import constants
import utils

class LVLH():
	
	def __init__(self, target, chaser):
		
		self.target = target 
		self.chaser = chaser 
		
		self.reset_mem()
		
	def reset_mem(self):
		
		self.trajectory_x = []
		self.trajectory_y = []
		
	def save_state(self, x, y):

		self.trajectory_x.append(x)
		self.trajectory_y.append(y)
		
	def inertial2orbital(self):
		"""
		Takes the intertial positions and velocity and turns that into the LVLH distances and relative speeds.
		
		This gives excellent estimations of the horizontal (ie vbar) distance, but fails in the rbar direction.		
		Issue with linearisation of the fact that the trajectory is a curve?
		
		.. note:: This *does work* when the chaser and target are very close
		
		:return: delta r, delta v (both 3d vector)
		
		.. note:: this method is currently unused...
		"""
		
		Rchaser = np.array([0,self.chaser.x, self.chaser.y]).reshape((1,3))
		Rtarget = np.array([0,self.target.x, self.target.y]).reshape((1,3))
		deltaR = (Rchaser - Rtarget)
		
		Vchaser = np.array([0,self.chaser.u, self.chaser.v]).reshape((1,3))
		Vtarget = np.array([0,self.target.u,  self.target.v]).reshape((1,3))
		deltaV = (Vchaser - Vtarget)
		
		Ht = np.cross(Rtarget, Vtarget)
		
		Roi = np.zeros([3,3])
		
		Roi[:,0] = utils.unit_vector(np.cross(Ht,Rtarget))
		Roi[:,1] = -utils.unit_vector(Ht)
		Roi[:,2] = -utils.unit_vector(Rtarget)
		
		omegaI = np.cross(Rtarget, Vtarget) / np.linalg.norm(Rtarget) / np.linalg.norm(Rtarget)
		
		deltaro = np.matmul(deltaR, Roi)

		
		OmegaR = np.cross(omegaI, deltaR)
		
		deltaVo = np.matmul(deltaV, Roi) - np.matmul(OmegaR, Roi)
		
		return deltaro, deltaVo
		
	def compute_distances(self, save_state=True):
		"""
		Computes the LVLH distances, i.e. delta horizontal distances and delta altitude.
		Returns the distances in km.
		"""

		v1 = np.array([self.target.x, self.target.y])
		v2 = np.array([self.chaser.x, self.chaser.y])
		
		delta_angle = utils.angle_between(v1, v2)
		
		dx = 2e-3 * np.pi * (constants.radiusE + self.target.get_alt() * 1e3) * delta_angle / (2. * np.pi)
		dx *= np.sign(self.chaser.get_posangle() - self.target.get_posangle())
		# Fixing the issue going to +0 + epsilon to -2pi - epsilon:
		try:
			if np.abs(dx-self.trajectory_x[-1]) > 20.:
				dx *= -1
		except:
			pass
		
		dy = -self.target.get_alt() + self.chaser.get_alt() # We multiplied by -1
		
		"""
		if np.isnan(dx) or np.isnan(dy):
			r = np.hypot(self.chaser.x, self.chaser.y)
			print self.chaser.Rp * (1. + self.chaser.e) / (r * self.chaser.e) - 1. / self.chaser.e
			print self.target.get_true_anomaly(), 'tgt'
			print np.sign(self.chaser.get_posangle() - self.target.get_posangle())
			print dx, dy
			exit()
		"""
		
		#if len(self.trajectory_x) > 0 and dx - self.trajectory_x[-1] < 0 and self.trajectory_x[-1] > 0:
		#	sign = -1.
		if save_state:
			self.save_state(dx, dy)
		
		return dx, dy
	