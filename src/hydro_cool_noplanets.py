from amuse.lab import *
from cooling_class import Cooling, SimplifiedThermalModelEvolver

class Hydro:
    def __init__(self, hydro_code, particles, Sun=Particles(0), impact_shell=Particles(0), converter=None, dt=1|units.hour):

        self.escapers=Particles(0)
        self.Sun = Sun
        self.particles = particles
        self.impact_shell = impact_shell

        self.box_size = 20000|units.AU
        #self.converter = nbody_system.nbody_to_si(1|units.MSun, system_size)
        self.converter = converter
        self.code = hydro_code(self.converter, mode="openmp")
        self.code.parameters.use_hydro_flag=True
        self.code.parameters.radiation_flag=False
        self.code.parameters.self_gravity_flag=True
        self.code.parameters.gamma=1.
        self.code.parameters.isothermal_flag=True
        self.code.parameters.integrate_entropy_flag=False
        self.code.parameters.timestep= dt
        #self.code.parameters.periodic_box_size = 100.*system_size
        #self.openboxsize=0.9*self.code.parameters.periodic_box_size/2
        self.code.parameters.courant=0.2    
        self.code.parameters.artificial_viscosity_alpha = 0.1 
        self.code.commit_parameters()

        self.code.dm_particles.add_particles(self.Sun)
        self.code.gas_particles.add_particles(self.particles)
        self.code.gas_particles.add_particles(self.impact_shell)

        self.attributes=["x", "y", "z", "vx", "vy", "vz", "mass", "u", "rho", "h_smooth", "pressure"]
        self.channel_from_hydro = self.code.gas_particles.new_channel_to(self.particles,
                                                                     attributes=self.attributes)
        self.channel_from_hydro_to_shell = self.code.gas_particles.new_channel_to(self.impact_shell,
                                                                     attributes=self.attributes)
        self.channel_to_hydro = self.particles.new_channel_to(self.code.gas_particles,
                                                              attributes=self.attributes)
        self.channel_from_shell_to_hydro = self.impact_shell.new_channel_to(self.code.gas_particles,
                                                              attributes=self.attributes)

        self.channel_from_Sun_to_hydro = self.Sun.new_channel_to(self.code.dm_particles)
        self.channel_from_hydro_to_Sun = self.code.dm_particles.new_channel_to(self.Sun)

        self.channel_from_hydro.copy()
        self.channel_from_hydro_to_shell.copy()

        # In order to be using the bridge
        self.get_gravity_at_point = self.code.get_gravity_at_point
        self.get_potential_at_point = self.code.get_potential_at_point
        self.get_hydro_state_at_point = self.code.get_hydro_state_at_point

        self.cooling = SimplifiedThermalModelEvolver(self.code.gas_particles)
        #self.cooling = Cooling(self.gas_particles)
        self.cooling.model_time=self.code.model_time
        
    @property
    def kinetic_energy(self):
        return self.code.get_kinetic_energy()
    @property
    def potential_energy(self):
        return self.code.get_potential_energy()

    @property
    def model_time(self):
        return self.code.model_time
        
    @property
    def gas_particles(self):
        return self.code.gas_particles

    def add_particles(self, impact_shell):
        #self.code.gas_particles.add_particles(impact_shell)
        self.particles.add_particles(impact_shell)
        self.particles.synchronize_to(self.code.gas_particles)
        #self.cooling.particles.add_particles(impact_shell)

    def remove_far_away_and_interacting_particles(self):
        escapers=self.particles.select_array(lambda x,y,z: (x**2+y**2+z**2 > self.box_size**2), ["x","y","z"])
        if len(escapers)>0:
            print "Remove escaping particles:", len(escapers)
            self.particles.remove_particles(escapers)
            self.particles.synchronize_to(self.code.gas_particles)

        amin2 = (0.1 |units.AU)**2    
        mergers=self.particles.select_array(lambda x,y,z: (x**2+y**2+z**2 < amin2), ["x","y","z"])
        if len(mergers)>0:
            print "Remove merging particles:", len(mergers)
            self.particles.remove_particles(mergers)
            self.particles.synchronize_to(self.code.gas_particles)
    
    @property
    def stop(self):
        return self.code.stop
        
    def evolve_model(self, model_time):
        model_time_old = self.model_time
        dt=model_time - model_time_old
        print "Evolve Hydrodynamics:", dt.in_(units.day)
        print "Cool gas for dt=", (dt/2).in_(units.day)

        self.remove_far_away_and_interacting_particles()
        
        print "cool gas"
        self.cooling.evolve_for(dt/2)

        print "run_hydro for ", dt.in_(units.yr), "to t=", (self.code.model_time + dt).in_(units.yr)
        self.code.evolve_model(self.code.model_time + dt)
        print "gas evolved to t=", self.code.model_time.in_(units.yr)

        print "Cool gas for another dt=", (dt/2).in_(units.day)
        self.cooling.evolve_for(dt/2)
        print "...done."

        self.channel_from_hydro.copy()
        self.channel_from_hydro_to_shell.copy()
        self.channel_from_hydro_to_Sun.copy()
        print "Hydro arrived at:", self.model_time.in_(units.day)

        
