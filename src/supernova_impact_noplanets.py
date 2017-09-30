import os
import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from optparse import OptionParser
from amuse.couple import bridge
#import make_planets_oligarch
from amuse.ext.orbital_elements import orbital_elements_from_binary
from make_impacting_shell_in_sph  import make_impacting_shell_in_sph
from make_impacting_shell_in_sph  import calculate_supernova_impact_shell_parameters

from hydro_cool import Hydro

def calculate_orbital_elements_on_single_object(star, planet):

    p = Particles()
    p.add_particle(star)
    p.add_particle(planet)
    M, m, a, e, ta_out, inc, lan_out, aop_out = orbital_elements_from_binary(p, G=constants.G)
    return a, e, inc

from amuse.ext.orbital_elements import new_binary_from_orbital_elements
def make_planetesimal_disk(Nplanetesimals, amin, amax):

    emax = 0.1
    imax = 0.1
    a = amin + (amax-amin)*numpy.random.random_sample(Nplanetesimals)
    e = emax*numpy.random.random_sample(Nplanetesimals)
    i = imax*numpy.random.random_sample(Nplanetesimals)
    ta = 0 #numpy.acos(np.random.uniform(0,2*numpy.pi,Nplanetesimals))
    loan = 0
    aof = 0
    mp = 0.1 | units.MEarth
    planetesimals = Particles(Nplanetesimals)
    for i, pi in enumerate(planetesimals):
        b = new_binary_from_orbital_elements(Mstar, mp, a[i], e[i], ta, inc[i],
                                             loan, aof, G=constant.G)
        pi.mass = mp
        pi.position = b[1].position
        pi.velocity = b[1].velocity

    return planetesimals

def initialize_star_and_planetary_system(Mstar, Ndisk, Mdisk, Rmin, Rmax):

    converter=nbody_system.nbody_to_si(Mstar, Rmin)
    disk_massfraction = Mdisk/Mstar
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter,
                              densitypower=1.5, 
                              Rmin=Rmin.value_in(units.AU), 
                              Rmax=Rmax.value_in(units.AU),
                              q_out=1.0, discfraction=disk_massfraction).result
    #disk.h_smooth= Rmin/Ndisk

    star = Particles(1)
    star.mass = Mstar
    star.radius = 1| units.RSun
    star.position = (0,0,0) | units.AU
    star.velocity = (0,0,0) | units.kms
    planets = make_planets_oligarch.new_system(Mstar, star.radius,
                                               Rmin, Rmax, Mdisk)
    star.add_particles(planets[0].planets)
    print star
    
    return star, disk

def read_planetary_system(filename): #lim, snapshot_id):
    star = Particles(0)
    SNsource = Particles(0)
    gasdisk = Particles(0)
    time = 0 | units.yr
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
        if len(bi)>1:
            gasdisk.add_particles(bi.copy())
        elif hasattr(bi, "mass"):
            star.add_particles(bi.copy())
        else:
            SNsource.add_particles(bi.copy())
    print "Read planetary system at time", time.in_(units.yr)
    return star, gasdisk, SNsource

def sphere_volume(r):
    return (4./3.)*numpy.pi*r**3

def create_impact_shell(source, target, n_shell, m_shell, r_shell, thickness, v_shell):
    
    from make_impacting_shell_in_sph import make_impacting_shell_in_sph
    impact_shell = make_impacting_shell_in_sph(source, target, n_shell, m_shell, r_shell, thickness, v_shell)#, inc)

    print "Supernova shell parameters:"
    print "v=", v_shell.in_(units.kms)
    print "r_shell=", r_shell.in_(units.AU)
    print "thickness=", thickness.in_(units.AU)

    yield_Fe60 = 5e-4 | units.MSun
    yield_Al26 = 5e-5 | units.MSun
    Fe60_in_shell = yield_Fe60*m_shell/source.sn_mass
    Al26_in_shell = yield_Al26*m_shell/source.sn_mass
    Fe60_pp = Fe60_in_shell/n_shell
    Al26_pp = Al26_in_shell/n_shell
    impact_shell.Fe60mass = Fe60_pp
    impact_shell.Al26mass = Al26_pp
    impact_shell.name = "shell"
    
    #impact_shell.rho += 1.e-15 | units.g/units.cm**3

    return impact_shell

def calculate_orbital_elements(star, minor_bodies):
    for pi in minor_bodies:
        a, e, inclination = calculate_orbital_elements_on_single_object(star, pi)
        pi.semimajor_axis = a
        pi.eccentricity = e
        pi.inclination = inclination
        #print "Planet", pi.name, "a=", a.in_(units.AU), e, inclination
        #print "Orbital elements of the debris."
    print "Mean orbital parameters for N=", len(minor_bodies), "objects: a=", minor_bodies.semimajor_axis.mean().in_(units.AU), "e=", minor_bodies.eccentricity.mean(), "i=", minor_bodies.inclination.mean() 
        
def main(n_shell, input_filename, t_end):

    ## --- create disk ---
    start_index     = int(input_filename.split("i")[1].split(".")[0])
    index = start_index
    star_and_planets, disk, source = read_planetary_system(input_filename)
    if not hasattr(star_and_planets, "name"):
        star_and_planets.name = "planet"
    if not hasattr(disk, "name"):
        disk.name = "disk"
    if not hasattr(source, "name"):
        source.name = "source"
    disk.Fe60mass = 0 | units.MSun
    disk.Al26mass = 0 | units.MSun

    name = source.name
    tss = source.time
    E_sn = 6.84e+51 | units.erg
    if "10a" in name:
        #for SN PS1-10a
        print "Supernova template: PS1-10a"
        source.sn_mass = 6.0 | units.MSun
    elif "12bku" in name:
        print "Supernova template: PS1-12bku" 
        source.sn_mass = 20.0 | units.MSun
    elif "11aof" in name:
        print "Supernova template: PS1-11aof"
        source.sn_mass = 23.5 | units.MSun
    else:
        print "Supernova template: PS1-10a (default)"
        source.sn_mass = 20.0 | units.MSun
    
    Mdisk = disk.mass.sum()
    Rmin = 1|units.AU
    Rmax = 100|units.AU
    Mstar = star_and_planets[0].mass
    converter=nbody_system.nbody_to_si(Mstar, Rmin)

    Rdisk = 1.5*Rmax
    r_shell = Rdisk
    star = star_and_planets[0]

    distance = (source.position-star.position).length()/source.ff
    print "distance=", distance.in_(units.parsec)

    m_shell = source.sn_mass * (1./4.)*(r_shell/distance)**2
#    thickness = min(distance, 10000. | units.AU)
    thickness = 200 | units.AU
    #v_shell = -10000. | units.kms
    v_shell = -1 * numpy.sqrt(2*E_sn/source[0].sn_mass)
    
    t_impact = -distance/v_shell
    print "Delay time:", t_impact.in_(units.yr), "at v=", v_shell.in_(units.kms)
    print "Time since supernova:", tss.in_(units.yr)
    t_impact -= tss
    if t_impact<0|units.day:
        print "Impact should have occurred earlier."
        print "But we adopt an impact time now."        
        t_impact = 0|units.day
    print "Impact time right now"
    t_impact = 0|units.day
    print "Supernova impact time:", t_impact.in_(units.yr)

#    t_end = t_impact + 6*(Rdisk+2*thickness)/abs(v_shell)
    print "t_end=", t_end.in_(units.yr)

    impact_shell = create_impact_shell(source, star, n_shell, m_shell, r_shell, thickness, v_shell)
    print "Masses:", source.sn_mass.in_(units.MSun), "shell mass:", m_shell.in_(units.MSun), impact_shell.mass.sum().in_(units.MSun), disk[0].mass.in_(units.MSun), disk.mass.sum().in_(units.MSun), impact_shell[0].mass.in_(units.MSun)
    
    
    hydro = Hydro(Fi, disk, impact_shell, converter=converter, dt=0.01|units.hour, Sun=star.as_set())
    planet_attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "semimajor_axis", "eccentricity", "inclination"]

    debris = Particles(0)
    
    index += 1
    filename = "supernova_impact_i{0:04}.amuse".format(index)
    write_set_to_file(source, filename, 'amuse', attribute_names=planet_attributes,
                      append_to_file=False)
    write_set_to_file(debris, filename, 'amuse', attribute_names=planet_attributes)
    write_set_to_file(star_and_planets, filename, 'amuse', attribute_names=planet_attributes)
    write_set_to_file(disk, filename, 'amuse', attribute_names=hydro.attributes)
    write_set_to_file(impact_shell, filename, 'amuse', attribute_names=hydro.attributes)

    Etot_init = hydro.kinetic_energy + hydro.potential_energy
    Etot_prev = Etot_init

    impacting_shell_added = False
    dt_diag = 1 | units.day
    t_diag = dt_diag
    time = 0.0 | t_end.unit
    dt = 1 | units.day
    while time < t_end:
        time += dt
        source.time = time

        hydro.evolve_model(time)

        Etot_prev_se = hydro.kinetic_energy + hydro.potential_energy

#        channel_to_planets.copy()
        #channel_from_hydro_to_framework.copy()

#        print "Orbital elements of the planets at t=", time.in_(units.yr)
#        calculate_orbital_elements(star_and_planets[0], star_and_planets[1:])
#        calculate_orbital_elements(star_and_planets[0], debris)
        
        if time>t_diag:
            t_diag += dt_diag
            index += 1
            filename = "supernova_impact_i{0:04}.amuse".format(index)
            write_set_to_file(source, filename, 'amuse', attribute_names=planet_attributes,
                              append_to_file=False)
            write_set_to_file(star_and_planets, filename, 'amuse', attribute_names=planet_attributes)
            write_set_to_file(debris, filename, 'amuse', attribute_names=planet_attributes)
            write_set_to_file(disk, filename, 'amuse', attribute_names=hydro.attributes)
            write_set_to_file(impact_shell, filename, 'amuse', attribute_names=hydro.attributes)

        Ekin = hydro.kinetic_energy 
        Epot = hydro.potential_energy
        Etot = Ekin + Epot
        print "T=", time.in_(units.day)
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

        if os.path.isfile("STOP"):
            os.remove("STOP")
            print "Terminate run at output:", filename
            return
        
    hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--n_shell", dest="n_shell", type="int", default = 100000,
                      help="number of particles in the shell []")
    result.add_option("--t_end", unit=units.day,
                      dest="t_end", type="float", default = 200|units.day,
                      help="end time of the simulation")
    result.add_option("-f", 
                      dest="input_filename", default = "planetary_system_i0100.amuse",
                      help="input filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


