from rbarvbar import sat, lvlh, plots
	
###################################################################################################
# Selecting options

import inspect


# time interval 
dt = 0.1
save_every = 3.
n_orbits_to_save = 4.98

# Target 
tgt = sat.Satellite(400, 400, "TGT")
tgt.set_init_pos(true_anomaly=0)
# Chaser
chaser = sat.Satellite(350, 300, "Chaser")
chaser.set_init_pos(true_anomaly=18.)

###################################################################################################
# Running the code

t = 0
nit = int(n_orbits_to_save * chaser.P / dt)
sen = save_every / dt 

distances = lvlh.LVLH(tgt, chaser)


for ii in range(nit):
	
	t = ii * dt
	
	if ii % sen == 0:
		do_save_state = True
	else:
		do_save_state = False

	tgt.step(dt, do_save_state)
	chaser.step(dt, do_save_state)
	if do_save_state:
		distances.compute_distances()

outdir = (inspect.stack()[0][1].split("/")[-1]).split(".py")[0]
pl = plots.Plot(tgt, chaser, distances, dt*sen, outdir=outdir)
pl.plot(time_ticks=10)
pl.make_movie(name=outdir)
