from simple_maize import *
from matplotlib.pylab import plot, show, xlim, ylim, axes
import numpy

from TD_maize import maize,display,generate_mtg, illuminate,plant_irradiance, process

# %pylab inline
# pylab.rcParams['figure.figsize'] = (5, 5)



g=generate_mtg(plant_area=8000,
               plant_height=156,
               rmax=0.355,
               skew=0.091,
               wl_int=0.08,
               wl_slp=0.01,
               w0_int=0.05,
               w0_slp=0.01,
               lm_int=0.5,
               lm_slp=-0.02,
               incli_base=72,
               incli_top=19.6,
               infl=41.6,
               pos_l=0.6,
               plant_orientation=0,
               phyllotactic_angle=172.704,
               phyllotactic_deviation=24,
               nb_leaf_segment=10,
               phytomer=16,
               stage=None,
               seed=1)

plant_irradiance(g)
# display(g,light=False,isolated=True)


# t=leaf_irradiance(g, isolated=False, illuminated=None)
# display(g,light=True,isolated=False)



# t=plant_irradiance(g, isolated=True, illuminated=None)

# t.to_csv('maxk_density6.csv')

# display(g,light=True,isolated=True)
