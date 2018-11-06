##
# probe-make.py
#
# Jaerin Lee
# Applied Superconductivity Laboratory
# Dept. of Electrical and Computer Engineering
# Seoul National University
##

r_min = 0
r_max = .7
r_step = 1000
z_min = 0
z_max = .7
z_step = 1000

filename = 'probe.txt'

with open(filename, 'w') as f:
    for r_idx in range(r_step):
        for z_idx in range(r_step):
            r = r_min + r_idx * (r_max - r_min) / r_step
            z = z_min + z_idx * (z_max - z_min) / z_step
            buf = "%f\t%f\n" % (r, z)
            f.write(buf)

