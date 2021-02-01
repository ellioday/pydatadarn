import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
import pydatadarn

vectors = pydatadarn.classes.griddata.GridData()
vectors.add_station("bks", 37.1, -77.95, "W")
#vectors.add_station("ade", 51.89, -176.63, "E")
#vectors.add_station("adw", 51.89, -176.63, "W")
#vectors.add_station("fhe", 38.859, -99.389, "E")
#vectors.add_station("fhw", 38.859, -99.389, "W")
#vectors.add_station("cve", 43.271, -120.358, "E")
#vectors.add_station("cvw", 43.271, -120.358, "W")
vectors.add_station("wal", 37.93, -75.47, "E")
vectors.add_data("2013/10/02 00:00:00", "2013/10/02 09:00:00")

time_i = "2013/10/02 06:58:00"

time_index = np.where(vectors.times == time_i)
mcolats = vectors.mcolats[time_index]
mlons = vectors.mlons[time_index]
kvecs = vectors.kvecs[time_index]
los_vs = vectors.los_vs[time_index]
pydatadarn.plotting.vector_plot(mcolats, mlons, kvecs, los_vs, time=time_i)
