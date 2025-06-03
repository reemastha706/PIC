
import os
import random
import nlopt
from grating_power_meep import main

# lattice parameter (um)
wg_min = 0.2
wg_max = 1

# waveguide width (fraction of lattice parameter)
h_min = 0.2
h_max = 2

# d_min = 0.5
# d_max = 1

# p_min = 0.5
# p_max = 1

opt = nlopt.opt(nlopt.LN_BOBYQA, 2)
opt.set_max_objective(main)
opt.set_lower_bounds([ wg_min, h_min])
opt.set_upper_bounds([ wg_max, h_max])
opt.set_ftol_abs(0.005)
opt.set_xtol_abs(0.02)
opt.set_initial_step(0.04)
opt.max_eval = 50

# random initial parameters
wg_0 = wg_min + (wg_max-wg_min)*random.random()  # rounded to 2nm
h_0 = h_min + (h_max-h_min)*random.random()
# p_0 = p_min + (p_max-p_min)*random.random()  # rounded to 2nm
# d_0 = d_min + (d_max-d_min)*random.random()


x = opt.optimize([wg_0, h_0])
maxf = opt.last_optimum_value()
print("optimum at wg={} um, h={} um ".format(x[0],x[1]))
print("maximum value = {}".format(maxf))
print("result code = {}".format(opt.last_optimize_result()))