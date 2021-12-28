import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from os import path, mkdir
from sys import argv, exit
from astropy import units as u
from astropy_healpix import HEALPix
from astropy.coordinates import Galactic, TETE, SkyCoord
import csv

# Configuration
OUTPUT_DIR = './footprint_maps'
NSIDE = 64
#open_clusters_data_file = './Kharchenko_catalog_OCs.dat'
open_clusters_data_file = './Kharchenko_selected.csv'
globular_clusters_data_file = './baumgardt_harris_GCs.csv'

class CelestialRegion:
    """Class to describe a region on sky, including its position and
    extend in an on-sky visualization"""

    def __init__(self,params={}):
        self.l_center = None
        self.b_center = None
        self.l_width = None
        self.b_height = None

        for key, value in params.items():
            if key in dir(self):
                setattr(self,key,value)

    def calc_healpixels_for_region(self,ahp):
        n_points = 500
        cosfactor = np.cos(self.b_center*3.141592/180.0)
        cosfactor2 = cosfactor*cosfactor
        halfwidth_l = (self.l_width / 2.0) / cosfactor
        halfheight_b = self.b_height / 2.0

        l_min = max( (self.l_center-halfwidth_l), 0 )
        l_max = min( (self.l_center+halfwidth_l), 360.0 )
        b_min = max( (self.b_center-halfheight_b), -90.0 )
        b_max = min( (self.b_center+halfheight_b), 90.0 )

        l = np.linspace(l_min, l_max, n_points) * u.deg
        b = np.linspace(b_min, b_max, n_points) * u.deg

        LL,BB = np.meshgrid(l, b)
        pointsincircle = np.where(np.sqrt(
            (LL.value-self.l_center)*(LL.value-self.l_center)/cosfactor2 +
            (BB.value-self.b_center)*(BB.value-self.b_center)
            ) <= halfheight_b)
        #coords = SkyCoord(LL, BB, frame=Galactic())
        coords = SkyCoord(LL[pointsincircle], BB[pointsincircle], frame=Galactic())
        self.pixels = ahp.skycoord_to_healpix(coords)

def generate_map():
    """Function to plot given pixel positions on a HEALpix map"""

    # Get user optional parameters
    options = get_args()

    # Initialze the HEALpix map
    ahp = HEALPix(nside=NSIDE, order='ring', frame=TETE())

    # Load data on the locations of regions of interest depending on user selections:
    if options['object_list'] == 'O':
        regions = load_open_cluster_data()
    elif options['object_list'] == 'G':
        regions = load_globular_cluster_data()

    # Use the HEALpix map together with the regions of interest to calculate the
    # HEALpixels included within those regions:
    print('Generating sky map...')
    for r in regions:
        r.calc_healpixels_for_region(ahp)

    # Build a map combining the pixel regions of all regions of interest:
    map = build_sky_map(regions)

    # Output the map data:
    output_sky_map(map,options)

def output_sky_map(map,options):

    if not path.isdir(OUTPUT_DIR):
        mkdir(OUTPUT_DIR)

    if options['object_list'] == 'O':
        title = 'Open Clusters'
        file_name = 'open_clusters_map'
    elif options['object_list'] == 'G':
        title = 'Globular Clusters'
        file_name = 'globular_clusters_map'

    fig = plt.figure(3,(10,10))
    hp.mollview(map, title=title)
    hp.graticule()
    #plt.tight_layout()
    plt.savefig(path.join(OUTPUT_DIR,file_name+'.png'))
    plt.close(3)

    hp.write_map(path.join(OUTPUT_DIR,file_name+'.fits'), map, overwrite=True)

    print('Output sky map data to '+OUTPUT_DIR+', files '+file_name+'.png & .fits')

def build_sky_map(regions):

    NPIX = hp.nside2npix(NSIDE)
    map = np.zeros(NPIX)

    for r in regions:
        map[r.pixels] += 1.0

    return map

def load_open_cluster_data():

    regions = []
    if path.isfile(open_clusters_data_file):
        if '.dat' in open_clusters_data_file:
            file_lines = open(open_clusters_data_file,'r').readlines()
            for line in file_lines[1:]:
                entries = line.replace('\n','').split()
                params = {'l_center': float(entries[5]), 'b_center': float(entries[6]),
                            'l_width': float(entries[9])*2.0, 'b_height': float(entries[9])*2.0}
                r = CelestialRegion(params)
                regions.append(r)

        elif '.csv' in open_clusters_data_file:
            with open(open_clusters_data_file, newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for i,row in enumerate(reader):
                    if i >= 1:
                        line = scan_for_quotes(row[0])
                        entries = line.split(',')
                        try:
                            params = {'l_center': float(entries[6]), 'b_center': float(entries[7]),
                                        'l_width': float(entries[10])*2.0, 'b_height': float(entries[10])*2.0}
                            r = CelestialRegion(params)
                            regions.append(r)
                        except ValueError:
                            pass
        print('Loaded data on Open Clusters')

    else:
        raise IOError('Cannot find data file for Open Clusters at '+open_clusters_data_file)

    return regions

def scan_for_quotes(line):
    if '"' in line:
        i0 = line.index('"')
        i1 = line[i0+1:].index('"')
        new_line = line[0:i0]+line[i0:i1].replace(',','_')+line[i1:]
    else:
        new_line = line 
    return new_line

def load_globular_cluster_data():

    regions = []
    if path.isfile(globular_clusters_data_file):
        with open(globular_clusters_data_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
            for i,row in enumerate(reader):
                if i >= 1:
                    entries = row[0].split(',')
                    try:
                        # What's the cluster radius we want to sample? Here we set it to 1.5*half-light radius
                        # The half-light radius is r_hl(pc) / R_sun(kpc) converted to deg
                        radius = 1.5*(0.001*float(entries[16])/float(entries[5]))*(180.0/3.141692)
                        # Another alternative: use the tidal radius r_t = r_c^c (with max at 30 arcmin, minimum at 2)
                        radius = float(entries[74])*np.power(10.0,float(entries[72])) / 60.0
                        if (radius>0.5): radius = 0.5
                        if (radius<2.0/60.0): radius = 2.0/60.0
                        
                        params = {'l_center': float(entries[3]), 'b_center': float(entries[4]),
                                  'l_width': 2.0*radius, 'b_height': 2.0*radius}
                        # Another alternative: use r_t = r_c^c (with max at 30 arcmin, 5 arcmin if unknown)
                        params = {'l_center': float(entries[3]), 'b_center': float(entries[4]),
                                  'l_width': float(entries[74])*np.power(10.0,float(entries[74])),
                                  'b_height': 1.5*(0.001*float(entries[16])/float(entries[5]))*2.0*(180.0/3.141692)}
                        r = CelestialRegion(params)
                        regions.append(r)
                    except ValueError:
                        pass
        print('Loaded data on Globular Clusters')
    else:
        raise IOError('Cannot find data file for Globular Clusters at '+globular_clusters_data_file)

    return regions

def get_args():
    """Present menu options and gather user selections"""

    options = {}
    if len(argv) == 1:
        menu = """Plot Open Cluster locations....................O
                  Plot Globular Cluster locations................G
                  Cancel.........................................C"""
        options['object_list'] = str(input('Please select an option to plot: ')).upper()
    else:
        options['object_list'] = str(argv[1]).upper()

    if options['object_list'] == 'C':
        exit()

    return options


if __name__ == '__main__':
    generate_map()
