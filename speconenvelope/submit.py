import os
import numpy as np


bexposure =       1500         ##  base exposure time.                                                                                                       
survey    =      'pfs'
run       =     'four'
nsig      =         1

for EW in [52.63, 11.00, -1.1, -14.92]: 
  ##  Scaling of base exposure time.
  for nexp in np.arange(38, 39, 1):
    for lw in [2400, 600]:
      os.system('python sn.py %f %d %d %d %s %s %d' % (EW, bexposure, nexp, nsig, survey, run, lw))
    
