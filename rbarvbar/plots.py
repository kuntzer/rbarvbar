import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import os

import constants
import utils as u

class Plot():
	
	def __init__(self, target, chaser, distances, dt, outdir="outputs"):
		
		self.target = target 
		self.chaser = chaser
		self.distances = distances 
		self.dt = dt
		#self.zero_speed = highlight_zero_speed
		#self.new_orbit = highlight_new_orbit
		self.outdir = outdir
		
		# Prep'ing the one-orbit trajectory
		self.idtgt = int(np.ceil(self.target.P / self.dt)*1.01)
		self.idch = int(np.ceil(self.chaser.P / self.dt)*1.01) 
		
	def plot(self, time_ticks=1, save=True, highlight_new_orbit=True):
		"""
		A lot of things here could be computed only once
		"""
		
		t = np.arange(len(self.target.speed_x)) * self.dt / self.chaser.P
		vtgt = np.hypot(self.target.speed_x, self.target.speed_y) * 1e-3
		vcha = np.hypot(self.chaser.speed_x, self.chaser.speed_y) * 1e-3
		
		#################################
		
		chaserspeed = np.array([self.chaser.speed_x, self.chaser.speed_y]).T
			
		vch = np.hypot(chaserspeed[:,0], chaserspeed[:,1])# * 1e-3
		chxy = np.array([self.chaser.trajectory_x, np.array(self.chaser.trajectory_y)])
		
		r = u.unit_vector(chxy).T
		
		aalpha = u.vangle_between(chaserspeed, r)
		
		dxo = 3. * np.pi * (self.target.a - self.chaser.a)
		
		v_bar = np.zeros_like(aalpha)
		for jjk in range(len(chxy[0])):
			try:
				vb = (self.distances.trajectory_x[jjk-1] - self.distances.trajectory_x[jjk]) / self.dt
			except:
				vb = np.nan 
			if jjk == 0:
				vb = np.nan
			v_bar[jjk] = vb
			
				
		vz_bar = vch * np.cos(aalpha) * 1e-3
		
		v_bar *= dxo / (np.nanmean((v_bar[:self.idch])) * 1e3 * self.chaser.P)
		
		##############################
		
		if self.chaser.hp == self.chaser.ha:
			txtch = r"$\mathrm{%s:\,%d\,km\,alt.}$" % (self.chaser.name, self.chaser.hp)
		else:
			txtch = r"$\mathrm{%s}:\,%d\times%d\mathrm{\,km\,alt.}$" % (self.chaser.name, self.chaser.hp, self.chaser.ha)
		txttgt = r"$\mathrm{%s:\,%d\,km\,alt.}$" % (self.target.name, self.target.hp)
		
		nimg = 0
		for ii in range(len(self.target.speed_x)-1):
			
			if not ii % time_ticks == 0: continue
			
			#if ii > len(self.target.speed_x): continue
	
			fig = plt.figure(figsize=(18,10))
			plt.subplots_adjust(wspace=0.0)
			plt.subplots_adjust(bottom=0.04)
			plt.subplots_adjust(top=0.92)
			plt.subplots_adjust(right=0.98)
			plt.subplots_adjust(left=0.)
			
			gs = gridspec.GridSpec(4, 3)
			
			# This where we plot the Earth orbit at the right scale
			ax0 = fig.add_subplot(gs[:2,0])
			#ax0.plot([1,2,3,4,5], [10,5,10,5,10], 'r-')
			ax0.set_aspect("equal")
			
			circle1 = plt.Circle((0, 0), constants.radiusE, color='b', alpha=0.5)
			ax0.add_artist(circle1)
			ax0.plot(self.target.trajectory_x[:self.idtgt], -1 * np.array(self.target.trajectory_y[:self.idtgt]), c='k')
			ax0.plot(self.chaser.trajectory_x[:self.idch], -1 * np.array(self.chaser.trajectory_y[:self.idch]), c='r')
			ax0.scatter(self.target.trajectory_x[ii], -1 * np.array(self.target.trajectory_y[ii]), c='k')
			ax0.scatter(self.chaser.trajectory_x[ii], -1 * np.array(self.chaser.trajectory_y[ii]), c='r')
			ax0.set_xticks([])
			ax0.set_yticks([])
			ax0.axis('off')
			
			ax1 = fig.add_subplot(gs[0,1:])
			#t = np.arange(len(self.target.speed_x)) * self.dt / self.chaser.P
			#vtgt = np.hypot(self.target.speed_x, self.target.speed_y) * 1e-3
			#vcha = np.hypot(self.chaser.speed_x, self.chaser.speed_y) * 1e-3
			ax1.plot(t, vtgt, c='k', label=self.target.name)
			ax1.plot(t, vcha, c='r', label=self.chaser.name)
			xmin1, xmax1 = ax1.get_xlim()
			ymin1, ymax1 = ax1.get_ylim()
			ax1.scatter(t[ii], vtgt[ii], c='k')
			ax1.scatter(t[ii], vcha[ii], c='r')			
			#ax1.plot(t, vcha-vtgt, c='b', label=r"$\Delta v$")
			ax1.legend()
			ax1.xaxis.set_ticklabels([])
			ax1.set_ylabel(r"$v\,\mathrm{[km/s]}$")
			ax1.grid()
			ax1.set_xlim([xmin1, xmax1])
			ax1.set_ylim([ymin1, ymax1])
			
			ax2 = fig.add_subplot(gs[2:,0])
			ax2.set_aspect("equal")
			
			r, theta = u.cart2polar(self.target.trajectory_x, self.target.trajectory_y)
			r -= constants.radiusE
			x, y = u.polar2cart(r, theta)
			ax2.plot(x[:self.idtgt] * 1e-3, -1e-3 * y[:self.idtgt], c='k')
			ax2.scatter(x[ii] * 1e-3, -1e-3 * y[ii], c='k')
			
			r, theta = u.cart2polar(self.chaser.trajectory_x, self.chaser.trajectory_y)
			r -= constants.radiusE
			x, y = u.polar2cart(r, theta)
			
			lim2 = 1.15 * np.amax([self.chaser.ha, self.target.ha])
			ax2.set_xlim([-lim2, lim2])
			ax2.set_ylim([-lim2, lim2])
			
			ax2.plot(x[:self.idtgt] * 1e-3, -1e-3 * y[:self.idtgt], c='r')
			ax2.scatter(x[ii] * 1e-3, -1e-3 * y[ii], c='r')
			ax2.set_xticklabels([])
			ax2.set_yticks([])
			ax2.axis('off')
			
			ax3 = fig.add_subplot(gs[1, 1:])#, sharex=ax1)
			ax3.set_xlabel(r"$\mathrm{Time\,[Orbits\,of\,chaser]}$")
			ax3.plot(t, v_bar*1e3, c='gold', label=r"$\bar{v}$")
			ax3.plot(t, vz_bar*1e3, c='g', label=r"$v_R$")
			ymin3, ymax3 = ax3.get_ylim()
			ax3.scatter(t[ii], v_bar[ii], c='gold')
			ax3.scatter(t[ii], vz_bar[ii], c='g')
			ax3.set_xlim([xmin1, xmax1])
			ax3.set_ylim([ymin3, ymax3])
			ax3.legend()
			ax3.grid()
			ax3.axhline(0, c='k', ls='--')
			ax3.set_ylabel(r"$\bar{v},\,v_R\,\mathrm{[m/s]}$")
			
			ax4 = fig.add_subplot(gs[2:,1:])
			
			ax4.spines['left'].set_position('zero')
			ax4.spines['right'].set_color('none')
			ax4.spines['bottom'].set_position('zero')
			ax4.spines['top'].set_color('none')
			
			# remove the ticks from the top and right edges
			ax4.xaxis.set_ticks_position('bottom')
			ax4.yaxis.set_ticks_position('left')
			ax4.invert_xaxis()
			ax4.invert_yaxis()
			
			ax4.plot(-np.array(self.distances.trajectory_x), -np.array(self.distances.trajectory_y))
			ax4.scatter(-np.array(self.distances.trajectory_x)[ii], -np.array(self.distances.trajectory_y)[ii])
			
			xmin, xmax = ax4.get_xlim()
			ymin, ymax = ax4.get_ylim()
			
			# For correctly plotting the vbar axis
			ax4.plot([0], [0], c='k')
			
			if t[ii] > 1 and highlight_new_orbit:
				ax4.axvline(-self.distances.trajectory_x[0], c='k', ls='--')
				ax4.axvline(-self.distances.trajectory_x[int(np.ceil(self.chaser.P / self.dt))], c='k', ls='--')
				x4t = -(self.distances.trajectory_x[0] + self.distances.trajectory_x[int(np.ceil(self.chaser.P / self.dt))])/2.
				y4t = (ymax - ymin) * 0.05 + ymin
				y4tt = (ymax - ymin) * 0.051 + ymin
				ax4.plot([-self.distances.trajectory_x[0], -self.distances.trajectory_x[int(np.ceil(self.chaser.P / self.dt))]], [y4t, y4t], ls='--', c='k')
				ax4.annotate(r'$1\,\mathrm{orbit\,}\Delta x = 3\pi\Delta a\approx %3.0f\,\mathrm{km}$' % (dxo/1e3), xy=(x4t, y4tt), ha='center')
			
			
			xx = ((1. - xmax / (xmax - xmin)) * 1.005 * (xmax - xmin) + xmin)
			yy = ((1. - ymax / (ymax - ymin)) * 1.015 * (ymax - ymin) + ymin)
			ax4.annotate(r'$\bar{v}\,\mathrm{[km]}$', xy=((xmax - xmin) * 0.01 + xmin, yy))
			ax4.annotate(r'$\bar{R}\,\mathrm{[km]}$', xy=(xx, (ymax - ymin) * 0.02 + ymin))
			
			gs.update(wspace=0.05, hspace=0.3)
			

			plt.suptitle(txtch + r"$\quad$--$\quad$" + txttgt + r"$\quad$--$\quad$" + "Time {:03d} min".format(int(self.dt*ii/60.)))
			
			if not os.path.exists(self.outdir):
				os.mkdir(self.outdir)
				
			if save:
				fig.savefig(os.path.join(self.outdir, "img_{:05d}.png".format(nimg)), dpi=200)
				plt.close()
			else:
				plt.show()
			
			nimg += 1
			
	def make_movie(self, name, system="linux"):
			
		if system == "linux":
			os.system("avconv -y -f image2 -i {}/img_%05d.png -c:v libx264 -r 6 -s hd720 -crf 16 {}/{}.mkv".format(self.outdir, self.outdir, name))
		else:
			raise ValueError("Unknown system")
