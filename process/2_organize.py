

import pandas as pd
import numpy as np
import geopandas as gpd
import os
import random
import rasterio
import numpy.ma as ma

# this function returns the levels and sfd csv files used to run the models in 3_make_fig3.R

def generate_out_and_fit():

    innerb, outerb = 100, 1100 # inside and outside buffer size surrounding target fire perimeter
    nobs_count = 50  # minimum number of pixel level observations

    # these csvs went through cleaning based on outputs from 1_sample_from_ee.py
    dfin = pd.read_csv("data/samples/inner.csv")
    dfout = pd.read_csv("data/samples/outer.csv")

    # step 2: put together the inner and outer observations only where the number of observations in > 50
    dfin = dfin[dfin.nobs >= nobs_count]
    dfout = dfout[dfout.nobs >= nobs_count]

    l = []
    for flower_year in [2018, 2019, 2020]:
        dfin_sub = dfin[dfin.flower_year == flower_year]
        nids = list(dfin_sub['nid'].unique())
        for target_nid in nids:
            dft_in = dfin[(dfin.target == target_nid) & (dfin.flower_year == flower_year)]
            dft_in['ys'] = dft_in['flower_year'] - dft_in['yof']
            dft_out = dfout[(dfout.target == target_nid) & (dfout.flower_year == flower_year)]
            dft_out = dft_out.dropna()
            # remove fires that occured after the target of interest, but we may be able to add this back in
            target_yof = dft_in['yof'].unique()[0]
            dft_out = dft_out[dft_out.yof <= target_yof]

            if len(dft_out) > 0 and len(dft_in) > 0:
                for tpi in list(dft_in['tpi_zone'].unique()):
                    dft_in_tpi = dft_in[dft_in.tpi_zone == tpi]
                    dft_out_tpi = dft_out[dft_out.tpi_zone == tpi]
                    m1 = pd.merge(dft_in_tpi, dft_out_tpi, on = 'tpi_zone', suffixes = ('_in', '_out'))
                    m1['f_prop_in'] = m1['yesflower_in'] / m1['nobs_in']
                    m1['f_prop_out'] = m1['yesflower_out'] / m1['nobs_out']
                    m1 = m1[['f_prop_in', 'f_prop_out', 'ys_in', 'ys_out', 'hsev_in', 'hsev_out',
                             'tpi_zone', 'yof_in', 'yof_out', 'flower_year_in', 'nid_in', 'nid_out', 'b_ub']]
                    l.append(m1)
    merged_inout = pd.concat(l)

    merged_inout['ys_diff'] = merged_inout['ys_in'] - merged_inout['ys_out']

    print('selecting in out fires')
    modis_checked = pd.read_csv('data/misc/fire_w_modis_check.csv')
    modis_checked = modis_checked['nid'].to_list()

    # selecting only fires which have a modis date check (accurate date information)
    merged_inout = merged_inout[merged_inout.nid_in.isin(modis_checked)]
    merged_inout = merged_inout[merged_inout.nid_out.isin(modis_checked)]

    # selecting only areas with a history of burning
    merged_inout = merged_inout[merged_inout.b_ub == 'burnt']

    ''' part 2: randomly sampling the dataset to find balance between ys, sev, and tpi '''

    print('sampling across ys, sev and tpi')
    l2 = []
    for ys in range(-15, 0 + 1):
        for sev in [1, 2, 3, 4, 5]:
            for tpi in [1, 2, 3, 4]:
                #print(ys, sev, tpi)
                tmp = merged_inout[merged_inout.ys_diff == ys]
                tmp = tmp[tmp.hsev_in == sev]
                tmp = tmp[tmp.tpi_zone == tpi]

                for sev_out in [1, 2, 3, 4, 5]:
                    tmp_out = tmp[tmp.hsev_out == sev_out]
                    if len(tmp_out) >= nobs_count:
                        samp = tmp_out.sample(nobs_count)
                    else:
                        samp = tmp_out

                    l2.append(samp)

    res_sampled = pd.concat(l2)

    # add on the clustered SEs here
    df_se = pd.read_csv("data/misc/SE_clustering_grids.csv")
    # dict grid2
    df_se2 = df_se[df_se.se == 'gridse2']
    dict2 = dict(zip(df_se2['nid'], df_se2['gridcell']))
    # dict grid3
    df_se3 = df_se[df_se.se == 'gridse3']
    dict3 = dict(zip(df_se3['nid'], df_se3['gridcell']))
    # dict grid4
    df_se4 = df_se[df_se.se == 'gridse4']
    dict4 = dict(zip(df_se4['nid'], df_se4['gridcell']))
    # dict grid5
    df_se5 = df_se[df_se.se == 'gridse5']
    dict5 = dict(zip(df_se5['nid'], df_se5['gridcell']))

    res_sampled['se2'] = res_sampled['nid_in'].map(dict2)
    res_sampled['se3'] = res_sampled['nid_in'].map(dict3)
    res_sampled['se4'] = res_sampled['nid_in'].map(dict4)
    res_sampled['se5'] = res_sampled['nid_in'].map(dict5)

    res_sampled.loc[res_sampled['hsev_in'] == 1, 'sev_binary_in'] = 1
    res_sampled.loc[res_sampled['hsev_in'] > 1, 'sev_binary_in'] = 2345
    res_sampled.loc[res_sampled['hsev_out'] == 1, 'sev_binary_out'] = 1
    res_sampled.loc[res_sampled['hsev_out'] > 1, 'sev_binary_out'] = 2345
    res_sampled['sev_binary_in'] = res_sampled['sev_binary_in'].astype(int)
    res_sampled['sev_binary_out'] = res_sampled['sev_binary_out'].astype(int)

    ''' part 3: packing it up and creating df_out and df_fit '''
    res_sampled['index'] = res_sampled.index
    res_samp_tmp = res_sampled
    l = []
    print("creating df out")
    for i, grp in res_samp_tmp.groupby('index'):
        obs = grp
        tmp_in = obs[['f_prop_in', 'hsev_in', 'sev_binary_in', 'ys_in', 'nid_in', 'tpi_zone', 'se2', 'se3', 'se4', 'se5', 'flower_year_in']]
        tmp_in.columns = ["f_prop", "hsev", 'sev_binary', "ys", "nid", "tpi_zone",  'se2', 'se3', 'se4', 'se5', 'flower_year']
        tmp_in['in_out'] = 1

        tmp_out = obs[['f_prop_out', 'hsev_out', 'sev_binary_out', 'ys_out', 'nid_out', 'tpi_zone',  'se2', 'se3', 'se4', 'se5', 'flower_year_in']]
        tmp_out.columns = ["f_prop", "hsev", 'sev_binary', "ys", "nid", "tpi_zone", 'se2', 'se3', 'se4', 'se5', "flower_year"]
        tmp_out['in_out'] = 2

        tmp_df = pd.concat([tmp_in, tmp_out])
        tmp_df['idx'] = i
        l.append(tmp_df)

    df_out = pd.concat(l)
    df_out = df_out.reset_index()

    dummies_binary = pd.get_dummies(df_out['sev_binary'].astype(str), prefix='cat')
    dummies_each = pd.get_dummies(df_out['hsev'].astype(str), prefix = 'cat')

    df_out = pd.merge(df_out, dummies_binary, left_index=True, right_index=True)
    df_out = pd.merge(df_out, dummies_each, left_index=True, right_index=True)

    df_out['cat_1'] = df_out['cat_1_x']
    df_out['cat_1_int'] = df_out['cat_1'] * df_out['ys']
    df_out['cat_2_int'] = df_out['cat_2'] * df_out['ys']
    df_out['cat_3_int'] = df_out['cat_3'] * df_out['ys']
    df_out['cat_4_int'] = df_out['cat_4'] * df_out['ys']
    df_out['cat_5_int'] = df_out['cat_5'] * df_out['ys']
    df_out['cat_2345_int'] = df_out['cat_2345'] * df_out['ys']

    # this will be the dfout csv
    out_fit_loc = "data/samples"
    df_out_outname = os.path.join(out_fit_loc,  "dfout_a50_g100-200.csv")
    df_out.to_csv(df_out_outname, index=False)

    print("creating df fit")
    l2 = []
    for i, grp in df_out.groupby('idx'):
        pass
        obs = grp
        tmp_in = obs[obs.in_out == 1]
        tmp_in = tmp_in[['f_prop', 'hsev', 'ys', 'cat_1', 'cat_2', 'cat_3', 'cat_4', 'cat_5', 'cat_2345',
                         'cat_1_int', 'cat_2_int', 'cat_3_int', 'cat_4_int', 'cat_5_int', 'cat_2345_int']]
        tmp_out = obs[obs.in_out == 2]
        tmp_out = tmp_out[['f_prop', 'hsev', 'ys', 'cat_1', 'cat_2', 'cat_3', 'cat_4', 'cat_5', 'cat_2345',
                         'cat_1_int', 'cat_2_int', 'cat_3_int', 'cat_4_int', 'cat_5_int', 'cat_2345_int']]
        stack = pd.concat([tmp_out, tmp_in])
        stack = stack.astype(float)
        stack = stack.diff()
        stack = stack.dropna()
        stack['se2'] = obs.iloc[0]['se2'].astype(int).astype(str)
        stack['se3'] = obs.iloc[0]['se3'].astype(int).astype(str)
        stack['se4'] = obs.iloc[0]['se4'].astype(int).astype(str)
        stack['se5'] = obs.iloc[0]['se5'].astype(int).astype(str)
        stack['flower_year'] = obs.iloc[0]['flower_year'].astype(int)
        stack['idx'] = i
        l2.append(stack)

    df_fit = pd.concat(l2)
    df_fit['ys'] = df_fit['ys'].astype(int)
    df_fit['cat_1'] = df_fit['cat_1'].astype(int)
    df_fit['cat_2'] = df_fit['cat_2'].astype(int)
    df_fit['cat_3'] = df_fit['cat_3'].astype(int)
    df_fit['cat_4'] = df_fit['cat_4'].astype(int)
    df_fit['cat_5'] = df_fit['cat_5'].astype(int)
    df_fit['cat_2345'] = df_fit['cat_2345'].astype(int)

    df_fit['cat_1_int'] = df_fit['cat_1_int'].astype(int)
    df_fit['cat_2_int'] = df_fit['cat_2_int'].astype(int)
    df_fit['cat_3_int'] = df_fit['cat_3_int'].astype(int)
    df_fit['cat_4_int'] = df_fit['cat_4_int'].astype(int)
    df_fit['cat_5_int'] = df_fit['cat_5_int'].astype(int)
    df_fit['cat_2345_int'] = df_fit['cat_2345_int'].astype(int)

    df_fit_outname = os.path.join(out_fit_loc, "dffit_a50_g100-200.csv")
    df_fit.to_csv(df_fit_outname, index=False)


if __name__ == "__main__":
    generate_out_and_fit()
