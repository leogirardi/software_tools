
# Summary of Leo's changes on late dec. 2021: 

1. Using extended density maps: 

In the configuration of `generate_galplane_priority_maps.py`:
```
#changed to use extended area map: STAR_MAP_FILE = 'TRIstarDensity_r_nside_64.npz'
STAR_MAP_FILE = 'TRIstarDensity_r_nside_64_ext.npz'
```
This means that the entire area of dec<40 deg will be used, and not the limited area used in pre-2018 discussions. Note: the `TRIstarDensity_*_nside_64_ext.npz` files are already included in MAF's `rubin_sim_data/maps/` folder.


2. Adding open clusters to `science_priority_regions.py`:
* Note, this is optional, indeed open clusters should be inserted here only if they are not included in some other metrics package! 

This involves adding a simple block to `science_priority_regions.py`, soon after reading the GCs:
```
    ##### ADDING open clusters:
    # ee module oc_all_lsst_field.py
    filterset_oc = { 'u': 0.0, 'g': 0.2, 'r': 0.3, 'i': 0.3, 'z': 0.2, 'y': 0.0 }
    oc_list = oc_all_lsst_field.fetch_OpenClusters_in_LSST_footprint()
    cluster0_pix = calc_hp_pixels_for_roundregion(oc_list[0]['l'], oc_list[0]['b'], oc_list[0]['rad'], 20, ahp)
    cluster1_pix = calc_hp_pixels_for_roundregion(oc_list[1]['l'], oc_list[1]['b'], oc_list[1]['rad'], 20, ahp)
    oc_regions = np.concatenate((cluster0_pix.flatten(), cluster1_pix.flatten()))
    for cluster in oc_list[2:]:
        cluster0_pix = calc_hp_pixels_for_roundregion(cluster['l'], cluster['b'], cluster['rad'], 20, ahp)
        oc_regions = np.concatenate((oc_regions, cluster0_pix.flatten()))
```
In addition to just adding open clusters, there are two other novelties here: 
- The open cluster list is read from a catalog, in the new module `oc_all_lsst_field.py`. It also reads the cluster radii. Eventually the same code could be used to read the globular cluster catalog (see the code just before the open clusters).
- The size of every cluster is given by the cluster radius, and not by the 3.5 degrees assumed for globular clusters.
- Pixel coordinates around the cluster centre are computed with the new function `calc_hp_pixels_for_roundregion`. Eventually, also this function could be adapted to be used on the globular clusters.

3. Then there are changes in `generate_sky_maps.py`, that were simply used to plot the maps of open and globular clusters, separately from the other priority regions.
