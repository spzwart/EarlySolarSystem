from multiprocessing import Pool, cpu_count
import os
import sys
import numpy
from matplotlib import pyplot
from amuse.lab import *

def calculate_disk_density_profile(bound_gas):
    R = [] | units.AU
    rho_mean = [] | units.g/units.cm**2
    nbin = 500
    for gi in range(len(bound_gas)-nbin):
        R.append(bound_gas[gi: gi+nbin].r.mean())
        S = (bound_gas[gi+nbin].r**2-bound_gas[gi].r**2)
        rho = bound_gas[gi: gi+nbin].mass.sum()/S
        rho_mean.append(rho)
    return R, rho_mean

def inclination_extrema(inc):
    n = numpy.sqrt(len(inc))
    a =numpy.histogram(inc, bins=numpy.arange(0, inc.max(), max(0.2, 0.1*inc.max()/n)))
    imax = 0
    amax = 0
    for i in range(len(a[0])):
        if a[0][i]>amax:
            amax = a[0][i]
            imax = i
    ileft = 0
    for i in range(imax):
        if a[0][i]>0.5*amax:
            ileft = a[1][i]
            break
    iright = 0
    for i in range(len(a[0])-imax):
        if a[0][imax+i]<0.5*amax:
            iright = a[1][imax+i]
            break
    ipeak = a[1][imax]
    return ileft, ipeak, iright

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

def selected_last_filenames(f, index, n):
    indices, filenames = zip(*sorted(zip(index, f)))
    return filenames[-n:], indices[-n:]
                             
def main(dirname, filename, n_last, rho_lim):
    if filename is None:
        filenames = os.listdir(dirname)
        f = []
        index = []
        rho = [] | units.g/units.cm**2
        for fi in filenames:
            if "supernova_impact_i" in fi:
                f.append(fi)
                indx = fi.split("_i")[2].split(".")[0]
                index.append(indx)
                rho.append(rho_lim)
        f, index = selected_last_filenames(f, index, n_last)
        
    else:
        f = [filename]
    print f, index
    ncpu = min(n_last, cpu_count())
    print "Ncpu=", ncpu
    pool = Pool(ncpu)
    params = zip(f, index, rho)
    pool.map(process_file, params) 

def post_process_shell(Sun, shell):
    shell.r = ((shell.x-Sun.x)**2 + (shell.y-Sun.y)**2).sqrt()
    s = shell.sorted_by_attributes("r")
    sid = s[s.r<100|units.AU]
    sid = sid[abs(sid.z)<10|units.AU]
    print "shell particles in disk:", len(sid)
    return sid

def extract_disk_edge(filename, index, rho_lim):
#    gas, shell = read_planetary_system(filename)
    source, planets, debris, gas, shell = read_planetary_system(filename)
    bound_gas = gas[gas.ecc<1.0]
    
    if len(shell)>0:
        #bound_shell = shell[shell.ecc<1.0]
        bound_shell = shell
        m_bound_shell = bound_shell.mass.sum()
        m_shell = m_bound_shell
    else:
        bound_shell = shell
        m_shell = 0 | units.MSun
        m_bound_shell = 0 | units.MSun

    ibmin, ibpeak, ibmax = inclination_extrema(bound_gas.inc)
    m_bound_gas = bound_gas.mass.sum()
    T_bound_gas = bound_gas.T.mean()
    full_disk = bound_gas.copy()
#    full_disk.add_particles(bound_shell)
    R, rho = calculate_disk_density_profile(full_disk)
    r_lim = R[0]
    for i in range(len(R)):
        if rho[i]<rho_lim:
            #print "lowest density at R=", R[i].in_(units.AU)
            r_lim = R[i]
            break

    inner_disk = full_disk[full_disk.r<r_lim]
    ibmin, ibpeak, ibmax = inclination_extrema(inner_disk.inc)
    ifmin, ifpeak, ifmax = inclination_extrema(full_disk.inc)
    dirname = os.getcwd().split('/')[-1]
     
    print "Disk edge: ", dirname, index, len(full_disk), "   ", len(inner_disk), "   ", len(bound_shell), "    ", full_disk.mass.sum().value_in(units.MSun), "   ", inner_disk.mass.sum().value_in(units.MSun), "   ", m_bound_shell.value_in(units.MSun), "   ", r_lim.value_in(units.AU), "    ", inner_disk.ecc.mean(), "    ", ifpeak, "    ", ibpeak
    sys.stdout.flush()

    return #index, len(full_disk), len(inner_disk), len(bound_shell), full_disk.mass.sum(), inner_disk.mass.sum(), m_bound_shell, r_lim, inner_disk.ecc.mean(), ifpeak, ibpeak


def process_file(params):
    filename, index, rho_lim = params
    
    f = filename.split("_")[2]
    outfile = "orbital_elements_"+f
    print "process ", filename, "to", outfile
    source, planets, debris, disk, shell = read_planetary_system(filename)
    Sun = source[source.name=="Sun"]
    bound_gas = process_data(Sun, disk)
    bound_shell = process_data(Sun, shell)
    print "Nbound gas= ", len(bound_gas)
    print "Nbound shell= ", len(bound_shell)
    sid = post_process_shell(Sun, shell)
    write_set_to_file(Sun, outfile, "amuse", append_to_file=False)
    write_set_to_file(bound_gas, outfile, "amuse")
    #write_set_to_file(bound_shell, outfile, "amuse")
    write_set_to_file(sid, outfile, "amuse")

    extract_disk_edge(outfile, index, rho_lim)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_last", type="int", default = 5,
                      help="number of files to process [%default]")
    result.add_option("-d", dest="dirname", default = "./",
                      help="input directory [%default]")
    result.add_option("-f", dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("--rho", unit=units.g/units.cm**2,
                      dest="rho_lim", default =  2.0 | units.g/units.cm**2,
                      help="limiting density [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

