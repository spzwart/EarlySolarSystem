import os
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from matplotlib import pyplot
from supernova_IIp_Lightcurve import Supernova_IIp

class RadHydro:
    def __init__(self, disk, Nray):

        self.Nray = 10**Nray
        self.index = 0
        self.ff = 0.1
        self.time = 0 | units.day
        self.disk = disk
        
        self.star = Particle()
        self.star.name = "Sun"
        self.star.mass = 1|units.MSun
        self.star.position = (0,0,0) | units.AU
        self.star.velocity = (0,0,0) | units.kms
        
        converter=nbody_system.nbody_to_si(1|units.MSun, 1|units.AU)
        
        self.hydro = Fi(converter)#, mode="openmp")
        self.hydro.parameters.use_hydro_flag=True
        self.hydro.parameters.radiation_flag=False
        self.hydro.parameters.self_gravity_flag=True
        self.hydro.parameters.gamma=1.
        self.hydro.parameters.isothermal_flag=True
        self.hydro.parameters.integrate_entropy_flag=False
        self.hydro.parameters.timestep=0.5 | units.hour #0.125 | units.day
        self.hydro.parameters.epsilon_squared=0.1 | units.AU**2
        self.hydro.parameters.courant=0.2    
        self.hydro.parameters.artificial_viscosity_alpha = 0.1 
    
        self.hydro.gas_particles.add_particles(self.disk)
        self.hydro.dm_particles.add_particle(self.star)
        self.hydro_to_disk = self.hydro.gas_particles.new_channel_to(self.disk)
        self.hydro_to_star = self.hydro.dm_particles.new_channel_to(self.star.as_set())
        self.disk_to_hydro = disk.new_channel_to(self.hydro.gas_particles)
        self.star_to_hydro = self.star.as_set().new_channel_to(self.hydro.dm_particles)
        self.hydro.evolve_model(1|units.hour)
        self.hydro_to_disk.copy()
        self.hydro_to_star.copy()

        self.rad = SPHRay(redirection="file")#, number_of_workers=24)#, debugger="gdb")
        self.rad.parameters.number_of_rays=self.Nray/(1|units.hour)
        self.rad.parameters.default_spectral_type=-3.
        self.rad.parameters.box_size=10000. | units.AU
        self.rad.parameters.ionization_temperature_solver=2
        print self.rad.parameters

        self.rad.gas_particles.add_particles(self.disk)

        self.gas_to_rad = self.disk.new_channel_to(self.rad.gas_particles)
        self.rad_to_gas = self.rad.gas_particles.new_channel_to(self.disk)

    def initialize_supernova(self, Rsn, inc, sn_name, ff=0.1):
        self.ff = ff
        delay_time = 1|units.day
        self.supernova_type = Supernova_IIp(sn_name, delay_time)

        self.Rsn = Rsn
        self.Msn = self.supernova_type.mass
        z = Rsn * numpy.cos(numpy.radians(inc))
        y = 0 | units.parsec
        x = numpy.sqrt(Rsn**2 - z**2)
        self.supernova=Particle()
        self.name = sn_name
        self.time = 0 | units.Myr
        self.supernova.time = self.time
        self.supernova.name = sn_name
        self.supernova.ff = ff
        self.supernova.position = (x.value_in(units.parsec),
                              y.value_in(units.parsec),
                              z.value_in(units.parsec)) |units.parsec
        self.supernova.position *= self.ff
        self.supernova.luminosity = self.ff**2 * self.supernova_type.luminosity_at_time(self.time)/(20.|units.eV)
        self.supernova.xion = 0.0
        self.supernova.u = (10**51 | units.erg)/self.Msn
        
        self.rad.src_particles.add_particle(self.supernova)

    def update_source_particle(self, time):
        self.supernova.time = time
        self.supernova.luminosity = self.ff**2 * self.supernova_type.luminosity_at_time(time)/(20.|units.eV)
        self.rad.src_particles[0].luminosity = self.supernova.luminosity
        
    def get_output_filename(self):
        nd = int(numpy.log10(len(self.disk)))
        nr = int(numpy.log10(self.Nray))
        rs = int(10*self.Rsn.value_in(units.parsec))
        nf = int(numpy.log10(self.ff))
        ms = self.Msn.value_in(units.MSun)
        print "n=", nd, nr, nf, rs, ms
        filename = "ID_Nd%1.1dNr%1.1dMs%2.2dRs0%1.1dpcff%1.1d_i%4.4d.amuse"%(nd, nr, ms, rs, nf, self.index)
        return filename

    def write_files(self):

        filename = self.get_output_filename()
        print "Writing file:", filename
        self.index += 1
        
        write_set_to_file(self.supernova.as_set(), filename, "amuse", append_to_file=False)
        write_set_to_file(self.star.as_set(), filename, "amuse")
        write_set_to_file(self.disk, filename, "amuse")
        
    def evolve_model(self, model_time):
        dt = model_time - self.time
        self.old_time = self.time
        self.time += dt/2

        self.rad.evolve_model(self.time)

        print "RT done at time:", self.time.in_(units.day)
        self.rad_to_gas.copy()

        self.time += dt/2

        self.disk_to_hydro.copy()
        self.star_to_hydro.copy()
        self.hydro.evolve_model(self.time)
        self.hydro_to_disk.copy()
        self.hydro_to_star.copy()

        self.update_source_particle(self.time)
        self.rad.evolve_model(self.time)
        print "RT done at time:", self.time.in_(units.day)
        self.rad_to_gas.copy()

    def print_diagnostics(self):
        umin = self.disk.u.min()
        umean = self.disk.u.mean()
        umax = self.disk.u.max()
        T = mu() / constants.kB * self.disk.u
        Tmin =  T.min()
        Tmean =  T.mean()
        Tmedian =  T.median()
        Tmax =  T.max()
        Tsorted = numpy.sort(T)
        Nhot = 0
        for ti in range(len(Tsorted)):
            if T[ti]>=1800|units.K:
                Nhot = ti
                break

        print "Time=", self.time.in_(units.day)
        print "Supernova luminosity:", (self.supernova.luminosity*(20.|units.eV)).in_(units.LSun)
        print "Ionization:", self.disk.xion.min(), self.disk.xion.mean(), self.disk.xion.max()
        print "Intenal energy:", umin, umean, umax
        print "Temperature:", Tmin, Tmean, Tmax, Tmedian
        print "Fraction hot gas:", len(self.disk), Nhot, "f=", 1.0*Nhot/len(self.disk)
        print "Density:", self.disk.density.min().in_(units.amu/units.cm**3), self.disk.density.mean().in_(units.amu/units.cm**3), self.disk.density.max().in_(units.amu/units.cm**3)
        print "scaleheight:", abs(self.disk.z.value_in(units.AU)).mean()
        
    def stop(self):
        self.rad.stop()
        self.hydro.stop()

        
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


def main(Ndisk, Mstar, Mdisk, Rin, Rout, dt, dt_diag, t_end, Nray, Rsn, inc, ff, sn_name):

    if dt_diag<dt:
        print "dt_diag<dt."
        print "STOP"
        return
    
    print "M=", Mdisk/Mstar    
    converter=nbody_system.nbody_to_si(Mstar, 1 | units.AU)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter,
                              Rmin=Rin.value_in(units.AU), 
                              Rmax=Rout.value_in(units.AU),
                              q_out=25.0, discfraction=Mdisk/Mstar).result
    print disk.x.max().in_(units.AU)
    print disk.mass.sum().in_(units.MSun)
    print disk.u.max().in_(units.kms**2)
    print disk.mass.min().in_(units.MSun)
    disk.name = "gas"
    disk.flux = 0. | units.s**-1
    disk.xion = 0.0

    radhydro = RadHydro(disk, Nray)
    radhydro.initialize_supernova(Rsn, inc, sn_name, ff)
#    vsn = 30000. | units.kms
#    t_end = Rsn/vsn

    radhydro.write_files()

    t_diag = dt_diag
    time = 0 | units.Myr
    while time<t_end:
        time += dt

        radhydro.evolve_model(time)
        radhydro.print_diagnostics()

        if os.path.isfile("STOP"):
            os.remove("STOP")
            outfile = radhydro.get_output_filename()
            print "Terminate run at output:", outfile
            return
        
        if time>t_diag:
            t_diag += dt_diag
            radhydro.write_files()

    radhydro.print_diagnostics()
    radhydro.write_files()
    radhydro.stop()

def plot_temperature(t, tmin, tmean, tmax):

    x_label = "t [day]"
    y_label = 'T [K]'
    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=14, ysize=8)
    pyplot.plot(t.value_in(units.day), tmean.value_in(units.K), c='k')
    pyplot.plot(t.value_in(units.day), tmin.value_in(units.K), c='r')
    pyplot.plot(t.value_in(units.day), tmax.value_in(units.K), c='b')
    pyplot.show()

def plot_ionization_fraction(pos, xion):
    r = []
    x = []
    for pi, xi in zip(pos, xion):
        #r.append(pi.length())
        r.append(numpy.log10(pi.value_in(units.AU)+0.000001))
        r.append(pi.value_in(units.AU))
        x.append(numpy.log10(xi+0.000001))
    r, x = zip(*sorted(zip(r, x)))
    
    from matplotlib import pyplot
    x_label = "r [pc]"
    y_label = r'$\xi_{\rm ion}$'
    figure = single_frame(x_label, y_label, logx=False, logy=False, xsize=14, ysize=8)
    pyplot.scatter(r, x, c=get_distinct(1), lw=0, s=100)
#    pyplot.xlim(-1, 1)
#    pyplot.ylim(-0.04, 1.19)
    #pyplot.savefig("fig_ionization_of_GMC")
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Ndisk", dest="Ndisk", type="int", default = 10000,
                      help="number of disk particles [%default]")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 2|units.yr,
                      help="radiation time [%default]")
    result.add_option("--dt_diag", unit=units.day,
                      dest="dt_diag", type="float", default = 3|units.day,
                      help="diagnostics output timestep [%default]")
    result.add_option("--dt", unit=units.day,
                      dest="dt", type="float", default = 1|units.day,
                      help="output timestep [%default]")
    result.add_option("--ff", 
                      dest="ff", type="float", default = 0.01,
                      help="distance adjustment factor [%default]")
    result.add_option("-r", unit=units.AU, type="float",
                      dest="Rin", default = 1.|units.AU,
                      help="inner disk radius [%default]")
    result.add_option("-R", unit=units.AU, type="float",
                      dest="Rout", default = 100.|units.AU,
                      help="outer disk radius [%default]")
    result.add_option("--Mstar", unit=units.MSun, type="float",
                      dest="Mstar", default = 1|units.MSun,
                      help="stellar mass")
    result.add_option("--Mdisk", unit=units.MSun, type="float",
                      dest="Mdisk", default = 0.01|units.MSun,
                      help="disk mass")
    result.add_option("--Rsn", unit=units.parsec, type="float",
                      dest="Rsn", default = 0.1|units.parsec,
                      help="supernova x-position")
    result.add_option("--inc", type="float",
                      dest="inc", default = 10,
                      help="supnova inclination wrt disk")
    result.add_option("--Nray", type="int",
                      dest="Nray", default = 6,
                      help="number of rays [%default]")
    result.add_option("--SNname", 
                      dest="sn_name", default = "PS1-12bku",
                                     #default = "PS1-11aof",
                      help="identify the supernova [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


