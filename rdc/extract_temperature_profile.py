from amuse.lab import *
from matplotlib import pyplot
import pickle
import os

def read_planetary_system(filename): #lim, snapshot_id):
    gas = Particles(0)
    shell = Particles(0)
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
      if hasattr(bi, "name"):
        if "gas" in bi.name:
            gas.add_particles(bi.copy())
        elif "shell" in bi.name:
            shell.add_particles(bi.copy())
    print "N=", len(gas), len(shell)
    return gas, shell

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

def process_data(filename):
    gas, shell = read_planetary_system(filename)
    t0 = 244 | units.day
    t00 = 244*3 | units.day
    if "ID_" in filename:
        t = 3*int(filename.split("_i")[1].split(".amuse")[0]) | units.day
#    elif "t_i" in filename:
#        t = t00 + (((int(filename.split("_i")[2].split(".amuse")[0])-t0) | units.day)
    elif "s_i" in filename:
        t = t00 + (((int(filename.split("_i")[1].split(".amuse")[0])-t0.value_in(units.day))) | units.day)
    else:
        print "filename cannot be processed."
        return
    print "time=", t

    search_bound = False
    if search_bound: 
        bound_gas = gas[gas.ecc<1.0]
    else:
        bound_gas = gas

    m_kB = mu() / constants.kB
    if not hasattr(gas, "T"):
        bound_gas.T =  (m_kB * bound_gas.u)
    bound_gas = bound_gas.sorted_by_attributes("T")

    Tt = [] | units.K
    it = []
    i = 0
    for gi in bound_gas:
        Tt.append(gi.T)
        it.append(i)
        i+=1

    f = open("Temperatures.pkl", 'awb')
    pickle.dump(t.value_in(units.day), f)
    pickle.dump(Tt.value_in(units.K), f)
    f.close()
        
def main(f):

    if f is None:
        f = os.listdir("./")
        print f
        for fi in f:
            if ".amuse" in fi:
                process_data(fi)
    else:
        process_data(f)

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="f", default = None,
                      help="output filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
