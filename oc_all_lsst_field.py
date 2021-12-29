# Dataset describing the locations in galactic coordinates of all known
# Open Clusters across the sky
# catalog produced by Leo Girardi and Simone Zaggia starting from Karchenko's
# code partially copied from Rachel Street generate_star_maps

import csv

def scan_for_quotes(line):
    if '"' in line:
        i0 = line.index('"')
        i1 = line[i0+1:].index('"')
        new_line = line[0:i0]+line[i0:i1].replace(',','_')+line[i1:]
    else:
        new_line = line 
    return new_line

def fetch_OpenClusters_in_LSST_footprint():
    oc_list = []
    open_clusters_data_file = './Kharchenko_selected.csv'

    with open(open_clusters_data_file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for i,row in enumerate(reader):
            if i >= 1:
                line = scan_for_quotes(row[0])
                entries = line.split(',')
                try:
                    params = {'name': entries[1], 'l': float(entries[6]), 'b': float(entries[7]),
                              'rad': float(entries[10]) }
                    oc_list.append(params)
                except ValueError:
                    pass
    print('Loaded data on '+str(len(oc_list))+' Open Clusters')

    return oc_list
