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

def main(filename):
    source, planets, debris, gas, shell = read_planetary_system(filename)
    print "mean Shell velocity:", shell.vx.mean().in_(units.kms), shell.vy.mean().in_(units.kms), shell.vz.mean().in_(units.kms)
    print "max Shell velocity:", shell.vx.max().in_(units.kms), shell.vy.max().in_(units.kms), shell.vz.max().in_(units.kms)
    
def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = None,
                      help="input filename [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
