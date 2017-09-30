import numpy
from matplotlib import pyplot
from amuse.lab import *
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.datamodel.rotation import rotate

"""
set_printing_strategy("custom",\
                      preferred_units = [units.MSun, units.parsec, units.Myr],\
                      precision = 6, prefix = "", separator = "[", suffix = "]")
"""
def travel_time(d, v):
    return d/v

def vstar(Eej, Mej):
    return numpy.sqrt(2*Eej/Mej)

def v_shell_ejecta(Esn, t):
    t = t/(1000.|units.yr)
    E51 = Esn/(1.e+51|units.erg)
    n = 1
    v0  = 1950 |units.kms
    v = v0 * E51**(1./5.) * n**(-1./5.) * t**(-3./5.)
    return  v

def r_shell_ejecta(Esn, t):
    t = t/(1000.|units.yr)
    E51 = Esn/(1.e+51|units.erg)
    n = 1
    r0  = 1.54e+19 | units.cm
    return  r0 * E51**(1./5.) * n**(-1./5.) * t**(2./5.)

def t_shell_ejecta(Esn, r):
    E51 = Esn/(1.e+51|units.erg)
    n = 1
    r0  = 1.54e+19 | units.cm
    t = (r.value_in(units.parsec)/(r0.value_in(units.parsec) * E51**(1./5.) * n**(-1./5.)))**(5./2.)
    return 1000*t | units.yr

def T_shell_ejecta(Esn, t):
    t = t/(1000.|units.yr)
    E51 = Esn/(1.e+51|units.erg)
    n = 1
    r0  = 1.54e+19 | units.cm
    T0 = 5.25e+7 | units.K
    T = T0 * E51**(2./5) * n**(-2./5.) * t**(-6./5.)  
    return T

def tstar(Eej, Mej, Rstar): 
    v = vstar(Eej, Mej)
    return Rstar/v

def rhostar(Mej, Rstar):
    return 3*Mej/(4*numpy.pi*Rstar**3)

def supernova_shell_density(Eej, Mej, Rstar, t_trav, t):
    rho = rhostar(Mej, Rstar)
    tstr = tstar(Eej, Mej, Rstar) 
    rho_e = rho * (tstr/t_trav)**3 * (t_trav/(t + t_trav))**3
    return rho_e

## 1987A: v=4000kms

## Pre-Sedov phase: 1990ApJ...356..549S
# adiabatic or sedov phase (10k t0 20k years)
def supernova_shell_temperature(Rs=0.1|units.parsec, gamma=7./5.):
    T0 = 9.e+7 | units.K
    Ts = T0 * (Rs/(10|units.AU))**(-3*(gamma-1))
    return Ts

def get_rot_matrix_from_v1_to_v2(v1, v2,
                                 eps=1.e-10):
  """
  rotation matrix to rotate vector v1 in the direction of the vector v2
  see e.g. http://math.stackexchange.com/questions/293116/rotating-one-3-vector-to-another
  """
  theta = numpy.arccos(numpy.dot(v1, v2)/ (numpy.linalg.norm(v1)*numpy.linalg.norm(v2)))
  if (numpy.abs(theta)<eps):
    rot_matrix = numpy.eye(3)
    return rot_matrix
  elif (numpy.abs(numpy.pi-theta)<eps):
    i_min = numpy.argmin(numpy.abs(v2))
    v2z = numpy.zeros_like(v2)
    v2z[i_min] = 1.
    v1_cross_v2z = numpy.cross(v1, v2z)
    x = v1_cross_v2z / numpy.linalg.norm(v1_cross_v2z)
  else:
    v1_cross_v2 = numpy.cross(v1, v2)
    x = v1_cross_v2 / numpy.linalg.norm(v1_cross_v2)
  A = numpy.array([[0.,   -x[2],  x[1]], \
                   [x[2],    0., -x[0]], \
                   [-x[1], x[0],    0.]])
  rot_matrix = numpy.eye(3) + numpy.sin(theta)*A + (1.-numpy.cos(theta))*numpy.dot(A,A)
  return rot_matrix

def rotate_disk(pos_list, vel_list, angular_momentum_vec, 
                disk_vec_ini=[0.,0.,1.]):
  angular_momentum_unit = angular_momentum_vec.value_in(units.kg*units.km**2/units.s) / \
                          angular_momentum_vec.length().value_in(units.kg*units.km**2/units.s)
  rot_mat = get_rot_matrix_from_v1_to_v2(disk_vec_ini, angular_momentum_unit)
  pos_unit = pos_list[0].unit
  vel_unit = vel_list[0].unit
  pos_rot_list = []
  vel_rot_list = []
  for pos_i, vel_i in zip(pos_list, vel_list):
    pos_rot = numpy.dot(rot_mat, pos_i.value_in(pos_unit)) | pos_unit
    vel_rot = numpy.dot(rot_mat, vel_i.value_in(vel_unit)) | vel_unit
    pos_rot_list.append(pos_rot)
    vel_rot_list.append(vel_rot)
  return pos_rot_list, vel_rot_list

def create_cylindrical_shell(N, mass, radius, thickness, velocity):
    
    new_slice = Particles(N)
    rho = numpy.sqrt(numpy.random.uniform(0,1, N))*radius
    phi = numpy.random.uniform(0,2.0*numpy.pi, N)

    new_slice.x = rho*numpy.sin(phi)
    new_slice.y = rho*numpy.cos(phi)
    new_slice.z = numpy.random.uniform(thickness.value_in(units.AU), 0, N) |units.AU

    new_slice.vx = 0.0 | units.kms    
    new_slice.vy = 0.0 | units.kms
    new_slice.vz = velocity
    
    new_slice.mass = mass/N
    new_slice.u = 0.5 * velocity**2

    print "N new in slice:", len(new_slice)
    return new_slice

def create_power_law_shell(N, mass, r_sn, radius, thickness, velocity):

    new_slice = Particles(N)
    rho = numpy.sqrt(numpy.random.uniform(0,1, N))*radius
    phi = numpy.random.uniform(0,2.0*numpy.pi, N)

    new_slice.x = rho*numpy.sin(phi)
    new_slice.y = rho*numpy.cos(phi)
#    new_slice.z = random_radial_distribution_of_particles(s=-12, N=N, rmin=thickness, rmax=r_sn) | units.parsec    
    new_slice.z = random_radial_distribution_of_particles(s=-12, N=N, rmin=thickness, rmax=2*thickness) | units.parsec    
    #new_slice.z = numpy.random.uniform(thickness.value_in(units.AU), 0, N) |units.AU

    new_slice.vx = 0.0 | units.kms    
    new_slice.vy = 0.0 | units.kms
    new_slice.vz = velocity

    new_slice.mass = (1.0*mass)/float(N)
    #new_slice.u = 0.5 * velocity**2
    #new_slice.u = 0.0 * velocity**2

    #alterntively use the sound speed (THoimas Wijnen)
    mu = 2.3 * constants.proton_mass #mean molecular weight
    print mu/(1|units.g)
    temp = supernova_shell_temperature(r_sn)
    print "Supernova shell temperature:", temp, "Shell particle mass=", new_slice[0].mass.in_(units.MSun)
    cs = ((temp*constants.kB)/mu).sqrt()
    #cs = 100|units.kms
    print "cs=", cs.in_(units.kms)
    new_slice.u = cs**2

    print "N new in slice:", len(new_slice)
    return new_slice
  
def random_radial_distribution_of_particles(s=-12, N=10, rmin=10|units.AU, rmax=1|units.parsec):
    """Power-law gen for pdf(r)\propto r^{s-1} for rmin<=r<=rmax"""
    r = numpy.random.random(size=N)
    ag, bg = rmin.value_in(units.parsec)**s, rmax.value_in(units.parsec)**s
    return (ag + (bg - ag)*r)**(1./s) - rmin.value_in(units.parsec)

def make_impacting_shell_in_sph(source, target, N, mass, radius, thickness, velocity):#, inclination):  

    r_sn  = source.position.in_(units.parsec).length()/source.ff

    f = (r_sn-1.1*radius)/r_sn
    locus_position_fraction = f*source.position
    print "lpf=", locus_position_fraction
    print locus_position_fraction.in_(units.parsec), locus_position_fraction.in_(units.parsec).length()+1.1*radius
    locus_pos_from_sn = (target.position-locus_position_fraction)

    inclination = numpy.arctan((100*source.x-target.x)/(100*source.z-target.z))
    
    print "inclination:", numpy.degrees(inclination)
    disk_pos = (0.1+numpy.sin(inclination)) * radius * source.position/source.position.length()     
    print "shell_pos=", (disk_pos/(1|units.AU)), "[AU]"

    impact_shell = create_power_law_shell(N, mass, r_sn, radius, thickness, velocity)
    impact_shell.name = "shell"
    impact_shell.xion = 1.0
    
    angular_momentum_vec = locus_pos_from_sn.value_in(units.parsec) | units.kg * units.parsec**2/units.s

    # rotate shell
    phi = numpy.radians(0)
    theta = inclination
    psi = numpy.radians(0)
    rotate(impact_shell, phi, theta, psi)

    impact_shell.position += disk_pos

    return impact_shell

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 1000,
                      help="number of diagnostics time steps [10]")
    result.add_option("--ff", dest="ff", type="float", default = 0.01,
                      help="supernova distance fudge factor [0.01]")
    result.add_option("--Rdisk", unit=units.AU,
                      dest="Rdisk", type="float",default = 100|units.AU,
                      help="radius of the shell part [%default]")
    result.add_option("--Mstar", unit=units.MSun,
                      dest="Mstar", type="float",default = 20|units.MSun,
                      help="mass of exploding star [%default]")
    result.add_option("--N_sph", 
                      dest="N_sph", type="int",default = 10000,
                      help="number of sph particles [%default]")
    result.add_option("--RSN", unit=units.parsec,
                      dest="RSN", type="float",default = 0.1|units.parsec,
                      help="distance to exploding star [%default]")
    result.add_option("--inc", 
                      dest="inc", type="float",default = 30,
                      help="inclination [%default]")
    return result

def calculate_supernova_impact_shell_parameters(distance, Eej,
                                                M_shell, r_shell):

    #M_shell = Mstar_PreSN - (1.4|units.MSun)
#    Eej = 1.e+51 | units.erg
    Rstar = 50 | units.RSun
    """
    print distance.value_in(units.parsec)
    t_shell = t_shell_ejecta(Eej, distance)
    T_shell = T_shell_ejecta(Eej, t_shell)
    r_shell = r_shell_ejecta(Eej, t_shell)
    v_shell = v_shell_ejecta(Eej, t_shell)
    print "shell:", t_shell.in_(units.yr), r_shell.in_(units.parsec), v_shell.in_(units.kms), T_shell
    """
    
    #v_shell = -vstar(Eej, M_shell)
    v_shell = -30000. | units.kms
    t = tstar(Eej, M_shell, Rstar)
    print "v=", v_shell.in_(units.kms)
    print "d=", distance.in_(units.parsec)
    t_trav = travel_time(distance, abs(v_shell))
    rho_shell = supernova_shell_density(Eej, M_shell, Rstar, t_trav, t)
    print "Travel time", t_trav.in_(units.yr)
    print  "Shell density:", rho_shell.in_(units.g/units.cm**3)
    area = numpy.pi*r_shell**2
    mass_faction = (r_shell/distance)**3
    a = area/rho_shell
    print "fm=", (M_shell*mass_faction).in_(units.MSun)
#    thickness = 1.0*((M_shell*mass_faction)/(area*rho_shell))
#    thickness = min(distance, 10000. | units.AU)
    thickness = 200.0 | units.AU
    print "Adopted shell thickness=", thickness.in_(units.AU)
    print "Shell passage timescale=", (thickness/abs(v_shell)).in_(units.day)

    disk_volume = numpy.pi*r_shell**2 * thickness
    m_disk = rho_shell*disk_volume
    print "m_disk=", m_disk.in_(units.MSun)

    return r_shell, thickness, v_shell, m_disk
  
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    r_shell = 2*o.Rdisk

    z_sn = o.RSN * numpy.cos(numpy.radians(o.inc))
    x_sn = numpy.sqrt(o.RSN**2 - z_sn**2)
    source = Particle()
    source.position = (0,0,0) | units.AU
    source.x = x_sn
    source.y = 0.00 |units.parsec 
    source.z = z_sn
    source.velocity = (0, 0, 0) | units.kms
    print source.position

    target = Particle()
    target.position = (0.0, 0.0, 0.0) | units.parsec

    distance = (source.position-target.position).length()/o.ff
    r_shell, thickness, v_shell, m_disk = calculate_supernova_impact_shell_parameters(distance, o.Mstar)
#    thickness *=10

    N = o.N_sph
    print "disk mass=", m_disk.in_(units.MSun), "N=", N
   
    impact_shell = make_impacting_shell_in_sph(source, target, N, m_disk, r_shell, thickness, v_shell, o.ff)

    pyplot.figure()
    pyplot.scatter(source.x.value_in(units.AU), source.y.value_in(units.AU), s=100, marker='o', c='y')
    pyplot.scatter(target.x.value_in(units.AU), target.y.value_in(units.AU), s=100, marker='o', c='r')
#    pyplot.scatter(impact_shell.x.value_in(units.AU), impact_shell.y.value_in(units.AU), s=100, marker='o', alpha=0.1)
    pyplot.scatter(impact_shell.x.value_in(units.AU), impact_shell.z.value_in(units.AU), s=100, marker='o', alpha=0.1)
#    for gi in impact_shell:
#      pyplot.arrow(gi.x.value_in(units.AU), gi.y.value_in(units.AU),
#                   gi.vx.value_in(units.parsec/units.Myr), gi.vy.value_in(units.parsec/units.Myr))
    gi = impact_shell[0]
#    pyplot.arrow(gi.x.value_in(units.AU), gi.y.value_in(units.AU),
#                 gi.vx.value_in(units.parsec/units.Myr), gi.vy.value_in(units.parsec/units.Myr))
    pyplot.arrow(gi.x.value_in(units.AU), gi.z.value_in(units.AU),
                 gi.vx.value_in(units.parsec/units.Myr), gi.vz.value_in(units.parsec/units.Myr))
    pyplot.plot([-100, 100], [0, 0])
    pyplot.xlabel("X")
    pyplot.ylabel("Z")
    pyplot.xlim(-3000, 3000)
    pyplot.ylim(-3000, 3000)
    pyplot.show()

    converter=nbody_system.nbody_to_si(impact_shell.mass.sum(), 100|units.AU)
    hydro=Fi(converter) #, number_of_workers=3)
 
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.integrate_entropy_flag=False
    #hydro.parameters.timestep=0.0625 | units.yr  
    hydro.parameters.timestep=0.1 | units.day

    hydro.gas_particles.add_particles(impact_shell)
    disk_attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "u", "rho", "h_smooth"]
    channel_from_hydro_to_framework = hydro.particles.new_channel_to(impact_shell,
                                                                     attributes=disk_attributes)
    channel_from_hydro_to_framework.copy()
    
    dt = 10 | units.day
    time = 0|units.yr
    t_end = 100|units.yr
    while time < t_end:
        time += dt
        hydro.evolve_model(time)
        channel_from_hydro_to_framework.copy()
        pyplot.scatter(impact_shell.x.value_in(units.AU), impact_shell.z.value_in(units.AU), s=100, marker='o', alpha=0.1)
        pyplot.plot([-100, 100], [0, 0])
        pyplot.xlabel("X")
        pyplot.ylabel("Z")
        pyplot.xlim(-200, 200)
        pyplot.ylim(-200, 200)
        pyplot.show()
