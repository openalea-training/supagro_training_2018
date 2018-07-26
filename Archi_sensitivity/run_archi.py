from simple_maize import *
from matplotlib.pylab import plot, show, xlim, ylim, axes
import numpy

from TD_maize import maize,display,generate_mtg, illuminate,plant_irradiance, process, leaf_irradiance

# %pylab inline
# pylab.rcParams['figure.figsize'] = (5, 5)



g5=generate_mtg(plant_area=6000,
               plant_height=184.0921,
               rmax=0.2815643,
               skew=0.1,
               wl_int=0.08,
               wl_slp=0.007,
               w0_int=0.05,
               w0_slp=0.01,
               lm_int=0.5,
               lm_slp=-0.02,
               incli_base=80.24322,
               incli_top=20,
               infl=43.87615,
               pos_l=0.6,
               plant_orientation=0,
               phyllotactic_angle=182.023 ,
               phyllotactic_deviation=30,
               nb_leaf_segment=10,
               phytomer=16,
               stage=None,
               seed=1)


g12=generate_mtg(plant_area=6000,
               plant_height=167.0215 ,
               rmax=0.2812295,
               skew=0.1,
               wl_int=0.08,
               wl_slp=0.007,
               w0_int=0.05,
               w0_slp=0.01,
               lm_int=0.5,
               lm_slp=-0.02,
               incli_base=78.36513,
               incli_top=20,
               infl=42.44517,
               pos_l=0.6,
               plant_orientation=0,
               phyllotactic_angle=187.4795,
               phyllotactic_deviation=30,
               nb_leaf_segment=10,
               phytomer=16,
               stage=None,
               seed=1)
# display(g,light=True, isolated=True, density=9, clear_sky=False, daydate='2000-06-21', longitude=3.52, latitude=43.36,
#             altitude=56, timezone="Europe/Paris")

# plant_irradiance(g, isolated=False, density=9, clear_sky=True, daydate='2000-06-21', longitude=3.52, latitude=0,
# altitude=56, timezone="Europe/Paris")

plant_irradiance(g, isolated=False, density=9, clear_sky=True)

plant_irradiance(g)
# display(g,light=True,isolated=True)


# t=leaf_irradiance(g)
# display(g12,light=True,isolated=False)
# display(w,light=True,isolated=False)


# t=plant_irradiance(g, isolated=True, illuminated=None)

# t.to_csv('maxk_density6.csv')

# display(g,light=True,isolated=True)
