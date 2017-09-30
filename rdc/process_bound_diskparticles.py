from matplotlib import pyplot
import numpy
from amuse.lab import *

#from read_protoplanetary_disk import read_protoplanetary_disk
def read_planetary_system(filename): #lim, snapshot_id):
    planets = Particles(0)
    debris = Particles(0)
    source = Particles(0)
    gas = Particles(0)
    shell = Particles(0)
    time = 0 | units.yr
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
      if hasattr(bi, "name"):
        if "gas" in bi.name:
            gas = bi.copy()
        elif "debris" in bi.name:
            debris = bi.copy()
        elif "Jupiter" in bi.name:
            planets = bi.copy()
        elif "shell" in bi.name:
            shell = bi.copy()
        else:
            source = bi.copy()
    print "N=", len(source), len(planets), len(debris), len(gas), len(shell)
    print "Read planetary system at time", time.in_(units.yr)
    return source, planets, debris, gas, shell


def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of particles in a gas)
    X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
    x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception("Error in calculating mu: mass fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)

from amuse.ext.orbital_elements import orbital_elements_from_binary

def calculate_orbital_elements(star, planet):
    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc_out, lan_out, aop_out = orbital_elements_from_binary(p, G=constants.G)
    return M, m, a, e, ta_out, inc_out, lan_out, aop_out

def get_component_binary_parameters(comp1, comp2, conv):
    kep = Kepler(redirection = "none")
    kep.initialize_code()

    print "m=", comp1.mass, comp2.mass
    mass = conv.to_nbody(comp1.mass + comp2.mass)
    pos = conv.to_nbody(comp2.position - comp1.position)
    vel = conv.to_nbody(comp2.velocity - comp1.velocity)
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    r = kep.get_separation()
    E,J = kep.get_integrals()	# per unit reduced mass, note
    kep.stop()

    return mass,a,e,r,E

def select_bound_particles(Sun, particles):
    result = []
    mass = particles.mass
    x_vector = particles.x
    y_vector = particles.y
    z_vector = particles.z

    vx_vector = particles.vx
    vy_vector = particles.vy
    vz_vector = particles.vz

    dx = Sun.x - x_vector
    dy = Sun.y - y_vector
    dz = Sun.z - z_vector
    dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
    dr = (dr_squared).sqrt()
    m_m = Sun.mass * mass

    dvx = Sun.vx - vx_vector
    dvy = Sun.vy - vy_vector
    dvz = Sun.vz - vz_vector
    dv_squared = (dvx * dvx) + (dvy * dvy) + (dvz * dvz)

    Ek = 0.5 * (Sun.mass + mass) * dv_squared
    binding_energy_for_all_pairs = Ek - constants.G*(m_m / dr)
    for bi in binding_energy_for_all_pairs:
        if bi<0|units.erg:
            print "bound pair."
    xxx
    bound = 0 | units.erg
    binding_energies = binding_energy_for_all_pairs[binding_energy_for_all_pairs<bound]
    bound_indices = numpy.arange(0, len(particles))[binding_energy_for_all_pairs<bound]
    bound_particles = particles[bound_indices]

    if len(binding_energies)>0 :
        for b,e,ib in zip(bound_particles, binding_energies,bound_indices) :
            #array of primary, secondary and energy
            result.append([particles,b,e]) 

    return result

def process_data(Sun, disk):
    gas = disk.copy()
#    gas = gas[gas.z<10|units.AU]
#    gas = gas[gas.z>-10|units.AU]
    gas.r = ((gas.x-Sun.x)**2 + (gas.y-Sun.y)**2).sqrt()
    s = gas.sorted_by_attributes("r")
    bound = Particles(0)
    print "Start processing orbital elements for N=", len(disk)
    for si in s:
        M, m, sma, ecc, ta, inc, lan, aop = calculate_orbital_elements(Sun, si)
        si.sma = sma
        si.ecc = ecc
        si.ta = ta
        si.inc = inc
        si.lan = lan
        si.aop = aop
        si.T =  (mu() / constants.kB * si.u)
        if si.ecc<1:
            bound.add_particle(si)
        #print "Oe:", si.a.in_(units.AU), si.e, si.i
    #bound = bound[bound.e<1]
    return bound

def X_main(dirname):
    for filename in os.listdir(dirname):
        process_file(filename)

from multiprocessing import Pool, cpu_count
import os
def main(dirname, filename):
    if filename is None:
        filenames = os.listdir(dirname)
        f = []
        for fi in filenames:
            if "supernova_impact_" in fi:
                f.append(fi)
            elif "ID_Nd5Nr7Ms23Rs06pcff-2_" in fi:
                f.append(fi)
    else:
        f = [filename]        
    print f
    ncpu = min(40, cpu_count())
    print "Ncpu=", ncpu
    pool = Pool(ncpu)
    pool.map(process_file, f) 

def post_process_shell(Sun, shell, zdisk):
    shell.r = ((shell.x-Sun.x)**2 + (shell.y-Sun.y)**2).sqrt()
    s = shell.sorted_by_attributes("r")
    sid = s[s.r<100|units.AU]
    sid = sid[abs(sid.z)<zdisk]
    print "shell particles in disk:", len(sid)
    return sid

def X_post_process_shell(Sun, shell):
    shell.r = ((shell.x-Sun.x)**2 + (shell.y-Sun.y)**2).sqrt()
    bound_particle(Sun, shell)
    s = shell.sorted_by_attributes("Etot")
    sid = s[s.r<0|units.erg]
    sout = shell - sid
    s = soutl.sorted_by_attributes("r")
    sout = s[s.r<10|units.AU]
    sout = out[abs(sid.z)<zdisk]

    sid.add_particles(sout)
    print "shell particles in disk:", len(sid)
    return sid

def bound_particle(Sun, p):
    for pi in p:
        Epot = -constants.G*(Sun.mass*pi.mass)/pi.r
        Ekin = 0.5*Sun.mass*(Sun.velocity-pi.velocity).length()**2
        Etot = Ekin+Epot
        pi.Etot = Etot

def process_file(filename):
    f = filename.split("_")[2]
    outfile = "orbital_elements_"+f
    print "process ", filename, "to", outfile

    source, planets, debris, disk, shell = read_planetary_system(filename)
    Sun = source[source.name=="Sun"]
    bound_gas = process_data(Sun, disk)
    bound_shell = process_data(Sun, shell)
    print "Nbound gas= ", len(bound_gas)
    print "Nbound shell= ", len(bound_shell)
    zdisk = 20|units.AU
#    sid = post_process_shell(Sun, shell)
    sid = post_process_shell(Sun, shell, zdisk)
    write_set_to_file(Sun, outfile, "amuse", append_to_file=False)
    write_set_to_file(bound_gas, outfile, "amuse")
    #write_set_to_file(bound_shell, outfile, "amuse")
    write_set_to_file(sid, outfile, "amuse")

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-d", dest="dirname", default = "./",
                      help="input directory [%default]")
    result.add_option("-f", dest="filename", default = None,
                      help="input filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

