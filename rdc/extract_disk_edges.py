import sys
import numpy
from amuse.lab import *

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

def read_planetary_system(filename): #lim, snapshot_id):
    planets = Particles(0)
    debris = Particles(0)
    source = Particles(0)
    gas = Particles(0)
    shell = Particles(0)
    time = 0 | units.yr
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
      if len(bi)>0:
        if "gas" in bi.name:
            gas.add_particles(bi.copy())
        elif "debris" in bi.name:
            debris.add_particles(bi.copy())
        elif "Sun" in bi.name:
            planets.add_particles(bi.copy())
        elif "shell" in bi.name:
            shell.add_particles(bi.copy())
        else:
            source.add_particles(bi.copy())
#    print gas[1]
#    print shell[1]
            
#    shell = gas.select(lambda n: "shell" in n,["name"])
#    gas -= shell
#    print shell[1]
    print "N=", filename, len(source), len(planets), len(debris), len(gas), len(shell)
#    print "Read planetary system at time", time.in_(units.yr)
    return source, planets, debris, gas, shell

def XX_read_planetary_system(filename): #lim, snapshot_id):
    gas = Particles(0)
    shell = Particles(0)
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
      if hasattr(bi, "name"):
        if "gas" in bi.name:
            gas.add_particles(bi.copy())
        elif "shell" in bi.name:
            shell.add_particles(bi.copy())
    return gas, shell

def process_file(params):
    filename, rho_lim = params
    index = filename.split("_i")[1].split(".")[0]
    print "process filename=", filename, "index=", index
    extract_disk_edge(filename, index, rho_lim)

from multiprocessing import Pool, cpu_count
import os
def main(filename, dirname, ncores, rho_lim):
    if filename is None:
        filenames = os.listdir(dirname)
        f = []
        rho = [] | units.g/units.cm**2
        for fi in filenames:
            if "orbital_elements_i" in fi:
                f.append(fi)
                rho.append(rho_lim)
        print "process:", f

        ncpu = min(ncores, cpu_count())
        print "Ncpu=", ncpu, "Njobs=", len(f)
        pool = Pool(ncpu)
        params = zip(f, rho)
        pool.map(process_file, params) 
    else:
        index = filename.split("_i")[1].split(".")[0]
        extract_disk_edge(filename, index, rho_lim)

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

def extract_disk_edge(filename, index, rho_lim):
#    gas, shell = read_planetary_system(filename)
    source, planets, debris, gas, shell = read_planetary_system(filename)
    if hasattr(gas, "ecc"):
        bound_gas = gas[gas.ecc<1.0]
    else:
        bound_gas = gas
    
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
        
    print "Disk edge: ", index, len(full_disk), "   ", len(inner_disk), "   ", len(bound_shell), "    ", full_disk.mass.sum().value_in(units.MSun), "   ", inner_disk.mass.sum().value_in(units.MSun), "   ", m_bound_shell.value_in(units.MSun), "   ", r_lim.value_in(units.AU), "    ", inner_disk.ecc.mean(), "    ", ifpeak, "    ", ibpeak
    
#    print "Disk edge:", index, "N=", len(full_disk), len(R[:i]), "M=", m_bound_gas.in_(units.MSun), m_bound_shell.in_(units.MSun), "R=", R[i].in_(units.AU), "ecc=", inner_disk.ecc.mean(), "inc=", inner_disk.inc.mean(), "T=", T_bound_gas, "Inclination limits full disk:", ifmin, ifpeak, ifmax, "Inclination limits dense disk:", ibmin, ibpeak, ibmax
    sys.stdout.flush()
    return 
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-c", dest="ncores", default = 40,
                      help="max number of cores used [%default]")
    result.add_option("-d", dest="dirname", default = "./",
                      help="input dirname [%default]")
    result.add_option("-f", dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("--rho", unit=units.g/units.cm**2,
                      dest="rho_lim", default =  2.0 | units.g/units.cm**2,
                      help="limiting density [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
