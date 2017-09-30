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

from hydro_cool_noplanets import Hydro

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
    planets = Particles(0)
    debris = Particles(0)
    source = Particles(0)
    gas = Particles(0)
    time = 0 | units.yr
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
        if "gas" in bi.name:
            gas.add_particles(bi.copy())
        elif "debris" in bi.name:
            debris.add_particles(bi.copy())
        elif "Jupiter" in bi.name:
            planets.add_particles(bi.copy())
        else:
            source.add_particles(bi.copy())

    shell = gas.select(lambda n: "shell" in n,["name"])
    gas -= shell
    print "N=", len(source), len(planets), len(debris), len(gas), len(shell)

    print "Read planetary system at time", time.in_(units.yr)
    return source, planets, debris, gas, shell

def sphere_volume(r):
    return (4./3.)*numpy.pi*r**3

def create_impact_shell(source, target, n_shell, m_shell, r_shell, thickness, v_shell, ff):
    
    from make_impacting_shell_in_sph import make_impacting_shell_in_sph
    impact_shell = make_impacting_shell_in_sph(source, target, n_shell, m_shell, r_shell, thickness, v_shell, ff)#, inc)

    print "Supernova shell parameters:"
    print "v=", v_shell.in_(units.kms)
    print "r_shell=", r_shell.in_(units.AU)
    print "m_shell=", m_shell.in_(units.MSun)
    print "thickness=", thickness.in_(units.AU)
    print "Total mass= ", impact_shell.mass.sum().in_(units.MSun)
    
    M_SN = (20-1.4) | units.MSun
    yield_Fe60 = 5e-4 | units.MSun
    yield_Al26 = 5e-5 | units.MSun
    Fe60_in_shell = yield_Fe60*m_shell/M_SN
    Al26_in_shell = yield_Al26*m_shell/M_SN
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
        
def main(n_shell, input_filename):

    ## --- create disk ---
    #start_index     = int(input_filename.split("i")[1].split(".")[0])
    start_index     = int(input_filename.split("i")[2].split(".")[0])
    index = start_index
    source, star_and_planets, debris, disk, impact_shell = read_planetary_system(input_filename)
    print "N read:", len(source), len(star_and_planets), len(debris), len(disk), len(impact_shell)
    dt_diag = 100 | units.day

    ff = source.ff
    name = source.name
    E_sn = 6.84e+51 | units.erg
    if "10a" in name:
        #for SN PS1-10a
        print "Supernova template: PS1-10a"
        m_shell = 6 | units.MSun
    elif "12bku" in name:
        print "Supernova template: PS1-12bku" 
        m_shell = 20.0 | units.MSun
    elif "11aof" in name:
        print "Supernova template: PS1-11aof"
        m_shell = 23.5 | units.MSun
    else:
        print "Supernova template: PS1-10a (default)"
        self.mass = 20 | units.MSun

    Mdisk = disk.mass.sum()
    Rmin = 1|units.AU
    Rmax = 100|units.AU
    Mstar = star_and_planets[0].mass
    converter=nbody_system.nbody_to_si(Mstar, Rmin)


    Sun = Particles(0)
    if len(star_and_planets)==1:
        Sun.add_particles(star_and_planets)
    """
    if len(star_and_planets)==1:
        from amuse.ext.solarsystem import new_solar_system
        first_contact = 2474649.5 |units.day
        solar_system = new_solar_system(first_contact)
        dx = (star_and_planets.position[0] - solar_system[0].position)
        dv = (star_and_planets.velocity[0] - solar_system[0].velocity)
        solar_system.position += dx
        solar_system.velocity += dv
        star_and_planets.add_particles(solar_system[5:])

    if (debris)<=1:
        print "Make debris disk."
        disk_massfraction = 0.01*Mdisk/Mstar
        Ndebris = 10000
        debris = ProtoPlanetaryDisk(Ndebris, convert_nbody=converter,
                                    Rmin=Rmin.value_in(units.AU), 
                                    Rmax=Rmax.value_in(units.AU),
                                    q_out=25.0, discfraction=disk_massfraction).result
        debris.position += star_and_planets[0].position
        debris.velocity += star_and_planets[0].velocity
        debris.mass = 0 | units.MSun
        debris.name = "debris"
        calculate_orbital_elements(star_and_planets[0], debris)
    """

    Rdisk = 1.5*Rmax
    r_shell = Rdisk
    star = star_and_planets[0]

    distance = (source.position-star.position).length()/ff
    thickness = min(distance, 10000. | units.AU)
    v_shell = -10000. | units.kms
    
    t_impact = -distance/v_shell
    tss = t_impact
    print "Delay time:", t_impact.in_(units.yr)
    print "Time since supernova:", tss.in_(units.yr)
    t_impact -= tss
    if t_impact<0|units.day:
        print "Impact should have occurred earlier."
        print "But we adopt an impact time now."        
        t_impact = 0|units.day
    print "Supernova impact time:", t_impact.in_(units.yr)

    t_end = t_impact + 6*(Rdisk+2*thickness)/abs(v_shell)
    print "t_end=", t_end.in_(units.yr)

    if len(impact_shell)==0:
        impact_shell = create_impact_shell(source, star, n_shell, m_shell, r_shell, thickness, v_shell, ff)
        #disk.add_particles(impact_shell)
        dt_diag = 1 | units.day
        impacting_shell_added = True
    else:
        impacting_shell_added = False

    hydro = Hydro(Fi, disk, Sun, impact_shell, converter, dt=0.01|units.hour)
    if len(sun_and_planets) == 1:
        gravity_hydro = Hydro
    else:
        st_gravity = ph4(converter)
        st_gravity.particles.add_particles(star_and_planets)

        pl_gravity = Huayno(converter, number_of_workers=6)
        KEPLER = 14
        pl_gravity.inttype = KEPLER
        pl_gravity.particles.add_particles(debris)

        planet_attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "semimajor_axis", "eccentricity", "inclination"]
        channel_to_planets = st_gravity.particles.new_channel_to(star_and_planets)
        channel_to_debris = pl_gravity.particles.new_channel_to(debris)

    index += 1
    filename = "supernova_impact_i{0:04}.amuse".format(index)
    write_set_to_file(source, filename, 'amuse', attribute_names=planet_attributes,
                      append_to_file=False)
    write_set_to_file(star_and_planets, filename, 'amuse', attribute_names=planet_attributes)
    write_set_to_file(debris, filename, 'amuse', attribute_names=planet_attributes)
    write_set_to_file(disk, filename, 'amuse', attribute_names=hydro.attributes)
    if len(impact_shell)>0:
        write_set_to_file(impact_shell, filename, 'amuse', attribute_names=hydro.attributes)

    if len(sun_and_planets) > 1:
        gravity_hydro = bridge.Bridge(use_threading=False)
        gravity_hydro.add_system(pl_gravity, (st_gravity, hydro,) )
        gravity_hydro.add_system(st_gravity, (hydro,) )
        gravity_hydro.add_system(hydro, (st_gravity,) )
        gravity_hydro.timestep = 1|units.hour


    Etot_init = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy
    Etot_prev = Etot_init

    t_diag = dt_diag
    time = tss
    dt = 1 | units.day
    while time < t_end:
        time += dt

        """
        if time>=t_impact and not impacting_shell_added:
            impacting_shell_added = True
            print "Add impact shell: N=", len(hydro.gas_particles)
            #disk.add_particles(impact_shell)
            hydro.add_particles(impact_shell)
            print "After adding the impact shell: N=", len(hydro.gas_particles)
            dt_diag = 1 | units.day
        """


        gravity_hydro.evolve_model(time-tss)

        Etot_prev_se = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy

        if len(sun_and_planets) > 1:
            channel_to_planets.copy()
            #channel_from_hydro_to_framework.copy()

            print "Orbital elements of the planets at t=", time.in_(units.yr), tss.in_(units.yr)
            calculate_orbital_elements(star_and_planets[0], star_and_planets[1:])
            calculate_orbital_elements(star_and_planets[0], debris)
        
        if time>t_diag:
            t_diag += dt_diag
            index += 1
            filename = "supernova_impact_i{0:04}.amuse".format(index)
            write_set_to_file(source, filename, 'amuse', attribute_names=planet_attributes,
                              append_to_file=False)
            write_set_to_file(star_and_planets, filename, 'amuse', attribute_names=planet_attributes)
            write_set_to_file(debris, filename, 'amuse', attribute_names=planet_attributes)
            write_set_to_file(disk, filename, 'amuse', attribute_names=hydro.attributes)
            if len(impact_shell)>0:
                write_set_to_file(impact_shell, filename, 'amuse', attribute_names=hydro.attributes)

        Ekin = gravity_hydro.kinetic_energy 
        Epot = gravity_hydro.potential_energy
        Etot = Ekin + Epot
        print "T=", time.in_(units.day)
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot

        if os.path.isfile("STOP"):
            os.remove("STOP")
            print "Terminate run at output:", filename
            return
        
    gravity_hydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--n_shell", dest="n_shell", type="int", default = 100000,
                      help="number of particles in the shell []")
    result.add_option("-f", 
                      dest="input_filename", default = "planetary_system_i0100.amuse",
                      help="input filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


