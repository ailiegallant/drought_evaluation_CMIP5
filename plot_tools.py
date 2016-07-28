#THIS MODULE CONTAINS TOOLS AND CODE FOR ASSISTING PLOTTING
#EACH SUBROUTINE IS SEPARATED BY A LINE --------
#-----------------------------------------------------------------------------

'''
This code sets specific levels according to a series of thresholds

INPUT
data - the data from which levels are derived
inc - the major increment from which levels should be derived
p - the percentiles below/above which data are clipped

'''
def levels(data, inc, p) :
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator