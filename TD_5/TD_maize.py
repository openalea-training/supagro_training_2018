import pandas
import multiprocessing
import sys

import openalea.plantgl.all as pgl
from alinea.astk.sun_and_sky import sun_sky_sources, sky_sources
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.light import light_sources
from simple_maize import parametric_leaf, simple_maize
from generator import cereals
from display import create_scene_from_mtg

import matplotlib as mpl
from matplotlib import pyplot


def colorScale(minval=0, maxval=1, label='relative irradiance'):
    '''    Produce an plot of the colorscale used by ViewMapOnCan when gamma == 1
    '''
    # Make a figure and axes with dimensions as desired.
    fig = pyplot.figure(figsize=(8, 1.5))
    ax1 = fig.add_axes([0.05, 0.4, 0.9, .5])

    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    cmap = mpl.cm.jet
    norm = mpl.colors.Normalize(vmin=minval, vmax=maxval)

    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
    cb1.set_label(label)

    return fig,


# default location and dates (Montpellier)
_daydate = '2000-06-21'
_timezone = 'Europe/Paris'
_longitude = 3.52
_latitude = 43.36
_altitude = 56

def maize(plant_area=7000,
          phytomer=16,
          plant_height=200,
          rmax=0.7,
          skew=0.01,
          wl_int=0.08,
          wl_slp=0.003,
          w0_int=0.5,
          w0_slp=0.01,
          lm_int=0.5,
          lm_slp=-0.02,
          incli_base=75,
          incli_top=15,
          infl=30,
          pos_l=0.6,
          plant_orientation=0,
          phyllotactic_angle=180,
          phyllotactic_deviation=10,
          nb_leaf_segment=10,
          seed=1,
          stage=None):
    """
    generate parameters of a maize plant

    Args:
        plant_area: plant area
        plant_height: plant height
        rmax: relative position of maximal area (from base to top)
        skew: skewness of the leaf area profile
        wl: width / length ratio
        incli_base: inclination angle at thae base of the plant (deg).
        angle at the top is set to 10
        infl: inflexion point of the logistic relation between insertion angle and leaf tip angle
        pos_l: relative position on the midrib where the rigidity decrease
        plant_orientation: azimuth (deg, positive clockwise from X+) of the first leaf
        phyllotactic_angle: phyllotactic angle between successive leaves (deg)
        phyllotactic_deviation: deviation amplitude from phyllotactic angle
        seed: the seed of the random number generator

    Returns:

    """

    # phytomer = 16
    phytomer = phytomer
    ranks = range(1, phytomer + 1)
    dinc = float(incli_top - incli_base) / (phytomer - 1)
    incli = [incli_base + i * dinc for i in range(phytomer)]
    # ddel = float(delta_angle_top - delta_angle_base) / (phytomer - 1)
    # delta_angle = [delta_angle_base + i * ddel for i in range(phytomer)]
    leaves = {
    # rank: parametric_leaf(nb_segment=nb_leaf_segment, insertion_angle=inc,
    #                       delta_angle=delta)
    # for rank, inc, delta in zip(ranks, incli, delta_angle)}
    rank: parametric_leaf(nb_segment=nb_leaf_segment, insertion_angle=inc, infl=infl, pos_l=pos_l, w0=w0_int + rank * w0_slp,
                          lm=lm_int + rank * lm_slp)
    for rank, inc in zip(ranks, incli)}

    # return simple_maize_raph(plant_area=plant_area, phytomer=phytomer, plant_height=plant_height,
    #                          rmax=rmax, leaves=leaves,
    #                          phyllotactic_angle=phyllotactic_angle,
    #                          phyllotactic_deviation=phyllotactic_deviation,
    #                          plant_orientation=plant_orientation, wl=wl, skew=skew,
    #                          seed=seed, stage=stage)

    return simple_maize(plant_area=plant_area, phytomer=phytomer, plant_height=plant_height,
                             rmax=rmax, leaves=leaves,
                             phyllotactic_angle=phyllotactic_angle,
                             phyllotactic_deviation=phyllotactic_deviation,
                             plant_orientation=plant_orientation,
                             wl_int=wl_int,
                             wl_slp=wl_slp,
                             w0_int=w0_int,
                             w0_slp=w0_slp,
                             lm_int=lm_int,
                             lm_slp=lm_slp,
                             skew=skew,
                             seed=seed, stage=stage)

def reader(data_file='rayostpierre2002.csv'):
    """ reader for mango meteo files """



def generate_mtg(**kwds):
    parameters = maize(**kwds)
    return cereals(leaf_volume=0, inclination=1, relative=True,plant=parameters)


def illuminate(g, isolated=True, density=9, clear_sky=False, daydate=_daydate, longitude=_longitude, latitude=_latitude,
               altitude=_altitude, timezone=_timezone):
    """ Illuminate a plant
    Args:
        isolated : is the plant isolated or within a canopy ?
        clear_sky: use clear_sky (homogeneous sky is used otherwise)
        irradiance: (float) sum of horizontal irradiance of all sources. If None
         diffuse horizontal clear_sky irradiance are used for clear_sky type and
          20% attenuated clear_sky global horizontal irradiances are used for
          soc and uoc types.
        dates: A pandas datetime index (as generated by pandas.date_range). If
            None, hourly values for daydate are used.
        daydate: (str) yyyy-mm-dd (not used if dates is not None).
        longitude: (float) in degrees
        latitude: (float) in degrees
        altitude: (float) in meter
        timezone:(str) the time zone (not used if dates are already localised)

    Returns:
        elevation (degrees), azimuth (degrees, from North positive clockwise),
        and horizontal irradiance of sources
    """
    if not clear_sky:
        light = light_sources(*sky_sources())
    else:
        sun, sky = sun_sky_sources(daydate=daydate, longitude=longitude, latitude=latitude, altitude=altitude,
                                   timezone=timezone, normalisation=1)
        light = light_sources(*sun) + light_sources(*sky)
    inter_row = 80
    inter_plant = 1. / density / (inter_row / 100.) * 100
    pattern = (-0.5 * inter_row, -0.5 * inter_plant,
               0.5 * inter_row, 0.5 * inter_plant)
    cs = CaribuScene(g, light=light, pattern=pattern, scene_unit='cm')
    raw, agg = cs.run(direct=True, simplify=True, infinite=not isolated)
    return cs, raw, agg


def display(g, light=False, illumination=None, **kwds):
    if not light:
        scene = create_scene_from_mtg(g)
        pgl.Viewer.display(scene)
    else:
        if illumination is None:
            cs, raw, agg = illuminate(g, **kwds)
        else:
            cs, raw, agg = illumination
        cs.plot(raw['Ei'], minval=0, maxval=1)


def plant_irradiance(g, illumination=None, **kwds):
    if illumination is None:
        _, _, agg = illuminate(g, **kwds)
    else:
        _, _, agg = illumination
    df = pandas.DataFrame(agg)
    leaves = [k for k,v in g.property('label').iteritems() if v.startswith('Leaf')]
    plant_area = df.area.sum()
    plant_leaf_area = df.loc[leaves,'area'].sum()
    dfl = df.loc[leaves,]
    return {'Ei': (df.Ei * df.area).sum() / plant_area,
            'Area': plant_area,
            'Area_leaf': plant_leaf_area,
            'Ei_leaf': (dfl.Ei * dfl.area).sum() / plant_leaf_area}

def leaf_irradiance(g, illumination=None, **kwds):
    if illumination is None:
        _, _, agg = illuminate(g, **kwds)
    else:
        _, _, agg = illumination
    df = pandas.DataFrame(agg)
    leaves = [k for k, v in g.property('label').iteritems() if v.startswith('Leaf')]
    dfl = df.loc[leaves,]
    return {'Area_leaf': dfl.area,
            'Ei_leaf': dfl.Ei}


def display_res(row, light=False):
    g_pars = ['plant_area', 'phytomer', 'seed', 'plant_height', 'rmax', 'skew', 'wl_slp', 'w0_slp', 'lm_slp',
              'incli_top', 'incli_base', 'infl', 'phyllotactic_angle', 'phyllotactic_deviation', 'stage']
    output = ['Ei', 'Ei_leaf', 'Area', 'Area_leaf']
    g = generate_mtg(**row[row.index[row.index.isin(g_pars)]].to_dict())
    if not light:
        scene = create_scene_from_mtg(g)
        pgl.Viewer.display(scene)
    else:
        kwds = row[row.index[~row.index.isin(g_pars + output)]].to_dict()
        cs, raw, agg = illuminate(g, **kwds)
        cs.plot(raw['Ei'], minval=0, maxval=1)


def run_sim(row, **kwds):
    g = generate_mtg(**row.to_dict())
    res = plant_irradiance(g, **kwds)
    for k in res:
        row[k] = res[k]
    for arg in kwds:
        row[arg] = kwds[arg]
    return row


def run_sim_xrun(xargs):
    row, kwds = xargs
    return run_sim(row, **kwds)

# ==============================================================================
# ==============================================================================


def process(path_input=None, path_output=None, nb_process=1,
            start=None, end=None, df_input=None, **kwds):

    if df_input is None:
        if path_input is not None:
            df_input = pandas.read_csv(path_input)
        else:
            df_input = pandas.DataFrame({'plant_height':[100,200], 'wl':[0.1,0.5]})

    rows = [row for index, row in df_input.iterrows()]

    if start is None:
        start = 0
    if end is None:
        end = len(rows)
    rows = rows[start:end]

    df_output = pool_function(rows, nb_process=nb_process, **kwds)

    if path_output is not None:
        df_output.to_csv(path_output)

    return df_output


def pool_function(rows, nb_process=2, verbose=True, **kwds):

    if nb_process <= 1:
        function = run_sim
        return run_function(function, rows, verbose=verbose, **kwds)

    function = run_sim_xrun
    pool = multiprocessing.Pool(nb_process)

    nb_rows = len(rows)
    df = pandas.DataFrame()
    args = [(row, kwds) for row in rows]
    for i, row in enumerate(pool.imap(function, args)):

        if verbose:
            print("%s : %d / %d" % (function.__name__, i, nb_rows))
            sys.stdout.flush()

        if row is not None:
            df = df.append(row)

    pool.close()
    pool.join()

    if verbose:
        print("%s : %d / %d" % (function.__name__, nb_rows, nb_rows))
        sys.stdout.flush()

    return df


def run_function(function, rows, verbose=True, **kwds):
    nb_rows = len(rows)
    df = pandas.DataFrame()
    for i, row in enumerate(rows):

        if verbose:
            print("%s : %d / %d" % (function.__name__, i, nb_rows))
            sys.stdout.flush()

        row = function(row, **kwds)
        if row is not None:
            df = df.append(row)

    if verbose:
        print("%s : %d / %d" % (function.__name__, nb_rows, nb_rows))
        sys.stdout.flush()

    return df


# ==============================================================================
# ==============================================================================

if __name__ == '__main__':
    # exp='ZA16'
    if len(sys.argv) > 1:
        # modulor config
        _, input, output, isolated, nbproc , density, latitude, daydate = sys.argv
        nbproc = int(nbproc)
        process(path_input=input, path_output=output, nb_process=nbproc, isolated=eval(isolated), density=density, latitude=latitude, daydate=daydate)
