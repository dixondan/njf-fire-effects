
import ee
ee.Initialize()
import pandas as pd
import numpy as np
import geopandas as gpd
import os
import random
import rasterio
import numpy.ma as ma

""" This script uses Google Earth Engine to sample surrounding proportional flowering within target and adjacent control fires """

# this function maps over year, nid (unique identifier, and an inside / outside buffer to find adjacent flowering counterfactuals

def return_fcover_by_nid(year, nid, inbuffer, outbuffer):

    # fire history composites (these are the same geotiffs as found in the directory ~/data/fire
    # it has public access
    fire_comp = ee.Image(
        'users/danieljdixon1991/NJF/PlanetScope_SWWA/v3/current_2022/fire_composites_2000_2019/fire_composite_{}'.format(year))

    fire_comp = fire_comp.unmask()
    # if years_since == 0, then it's unburnt, convert band 1 (ds, days since fire) to years since fire
    ys = fire_comp.select('ds').divide(365).ceil().rename('ys')
    fire_comp = fire_comp.addBands(ys)

    # masking out all 6, 7 and 8s because their severity scores are either unavailable or on the edge of two fires
    sev = fire_comp.select('sev')
    im8mask = sev.lte(5)
    fire_comp = fire_comp.mask(im8mask)

    # fire polygon
    # this is a fire perimeter dataset (DBCA 060) with all the metadata we need include NID (numerical identifiers)
    fires = ee.FeatureCollection(
        'users/danieljdixon1991/NJF/PlanetScope_SWWA/v3/current_2022/DBCAFireHistory_2000_2021_all_fires_explode_wID_edited_geometry_wNID_wModisDate_1kmgrid')

    # because we are mapping over nid values, we go fetch that nid to use as the target (treatment) fire
    fire = ee.Feature(fires.filterMetadata('nid', 'equals', nid).first())
    # make a 1000 m buffer aroud that perimeter (or however large is needed)
    fire_poly_1 = fire.buffer(inbuffer) #1000
    fire_poly_2 = fire.buffer(outbuffer) #2000 we want all fire information wihtin 1000 and 2000 m from the perimter
    inner_buffer = fire.buffer(-30) # remove edge effects
    outer_buffer = fire_poly_2.symmetricDifference(fire_poly_1, 1) # essentially do an inside clip with buffer 1 and 2

    flower = ee.Image('projects/euc-fire-sw-2022/assets/stoten2023_repo/flowering_{}'.format(year))\
                        .rename('flower_binary') # go grab the flowering data -- this is the same as ~/data/flowering/...
    # clip that data to the new buffer
    flower = flower.clip(fire_poly_2)

    # bring in tpi data publically available on ee
    tpi = ee.Image('users/danieljdixon1991/NJF/PlanetScope_SWWA/v3/current_2022/tpi_6m_1234at10m')

    # we have to clip fire_composite within the outer buffer -- have to do these inner and outer calculations
    # separately
    fire_comp_OUT = fire_comp.clip(outer_buffer)

    ''' function that vectorizes by nid for outer calculation of flower counts '''
    def vectorize_nid(n):
        classImage = fire_comp_OUT.select('nid').eq(ee.Number(n))
        vectors = classImage.updateMask(classImage) \
            .reduceToVectors(reducer=ee.Reducer.countEvery(),
                             geometry=outer_buffer.geometry(),
                             scale=30,
                             maxPixels=1e8).geometry()
        return ee.Feature(vectors, {"nid": n})

    classes = ee.List([0] + [int(i) for i in fireNIDs])
    classes = classes.map(vectorize_nid) # here we map over nid values and vectorize them, so we can use them to extract more information
    nid_vectors = ee.FeatureCollection(classes)

    # organize the flowering data to be two bands with 1 for flowering or not flowering into yesno_flower
    yes_flower_mask = flower.eq(1)
    yes_flower = flower.updateMask(yes_flower_mask).rename('yesflower')

    no_flower_mask = flower.eq(0)
    no_flower = flower.updateMask(no_flower_mask).rename('noflower')

    yesno_flower = yes_flower.addBands(no_flower) # we want the count of 1 and 0 flowering

    outpath = 'data' # or create an out location
    outname = os.path.join(outpath, "outer_{}_{}.csv".format(year, nid))

    l1 = [] # list1 takes outside flower counts for each unique nid, fire severity and tpi zone
    l2 = [] # list2 takes inside flower couns for each unique nid, fire severity and tpi zone
    for tpi_zone in [1, 2, 3, 4]:
        for hsev in [1, 2, 3, 4, 5]:
            print(tpi_zone, hsev)
            ''' need to mask out areas with sev 8: these were fires but don't have a severity '''
            ''' also need to mask out depending on high/low severity and tpi zone '''
            # sev 8 is already masked out
            # # mask tpi_zone first
            mask1_tpizone = tpi.eq(tpi_zone)
            yesno_flower_masked1 = yesno_flower.updateMask(mask1_tpizone)

            hsev_img = fire_comp.select('sev')

            mask2_hsev = hsev_img.eq(hsev)#.focal_max(radius = 15, units = 'meters') #.buffer(-15) # including a 15 m buffer

            yesno_flower_masked2 = yesno_flower_masked1.updateMask(mask2_hsev)

            yesno_flower_masked2_inner = yesno_flower_masked2.clip(inner_buffer)
            yesno_flower_masked2_outside = yesno_flower_masked2.clip(outer_buffer)

            counts_inside = yesno_flower_masked2_inner.reduceRegions(
                collection=inner_buffer,
                reducer=ee.Reducer.count()).getInfo() # this may cause memory issues calling .getInfo()

            for c in counts_inside['features']:
                tmp = c['properties']
                tmp = pd.DataFrame(tmp, index=[0])
                tmp = tmp[['nid', 'noflower', 'yesflower']]
                tmp['tpi_zone'] = tpi_zone
                tmp['hsev'] = hsev
                tmp['year'] = year
                l2.append(tmp)

            if hsev <= 5:
                counts_outside = yesno_flower_masked2_outside.reduceRegions(
                    collection=nid_vectors,
                    reducer=ee.Reducer.count()).getInfo() # again potential memory issues
                for c in counts_outside['features']:
                    tmp = c['properties']
                    tmp = pd.DataFrame(tmp, index = [0])
                    tmp = tmp[['nid', 'noflower', 'yesflower']]
                    tmp['tpi_zone'] = tpi_zone
                    tmp['hsev'] = hsev
                    tmp['year'] = year
                    l1.append(tmp)
            else:
                pass

    outpath = "data\fire_data_from_ee_wZero_wNID_wBuffer_{}_{}".format(
        inbuffer, outbuffer)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    res_inner = pd.concat(l2)
    res_inner['inout'] = 'inner'
    res_inner.to_csv(os.path.join(outpath, "inner_{}_{}.csv".format(year, nid)), index=False)

    res_outer = pd.concat(l1)
    res_outer['inout'] = 'outer'
    res_outer.to_csv(os.path.join(outpath, "outer_{}_{}.csv".format(year, nid)), index=False)

    print("completed", year, nid)




# gather a list of possible NID (unique numerical identifiers linking each fire ID)
l=[]
for year in [2018, 2019, 2020]:
    raster = "data\fire\fire_{}.tif".format(year)
    with rasterio.open(raster) as src:
        arr = src.read()
        sev = arr[1, :, :]
        nid = arr[2, :, :]
        nid = ma.masked_where(sev > 5, nid)#.data
        unique, counts = np.unique(nid, return_counts=True)
        s = pd.DataFrame(list(unique), columns = ['nid'])
        s['counts'] = pd.Series(list(counts))
        s = s.dropna()
        s = s[s['nid'].apply(lambda x: str(x).isdigit())]
        s['nid'] = s['nid'].astype(int)
        s = s[s.nid != 0]
        s['year'] = year
        s = s[s.counts > 50]
        l.append(s)
df = pd.concat(l)
fireNIDs = list(df['nid'].unique())

# now need to start mapping that function
for year in [2018, 2019, 2020]:
    for nid in fireNIDs:
        print(year, nid)
        for innerb, outerb in [(100, 1100)]:
            print(innerb, outerb)
            return_fcover_by_nid(int(year), int(nid), inbuffer=innerb, outbuffer=outerb)

