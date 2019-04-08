import  os 
import  numpy as np


RUN   =   'four'
runs  =      []

##   survey, program, target, repeat
#runs.append(['desi',  'dark', 'lbg', 1])
#runs.append(['beast', 'dark', 'lbg', 1])
runs.append(['pfs',   'dark', 'lbg', 1])
#runs.append(['desi',  'dark', 'BC03', 1])                                                                                                            
#runs.append(['beast', 'dark', 'BC03', 1]) 
#runs.append(['pfs',   'dark', 'BC03', 1])
#runs.append(['desi',  'dark', 'elg', 1])
#runs.append(['beast', 'dark', 'elg', 1])
#runs.append(['desi',  'dark', 'qso', 1])
#runs.append(['beast', 'dark', 'qso', 1])
#runs.append(['desi',  'dark', 'lrg', 1])
#runs.append(['beast', 'dark', 'lrg', 1])
#runs.append(['desi',  'dark', 'bgs', 1])                                                                                                                  
#runs.append(['beast', 'dark', 'bgs', 1]) 
#runs.append(['desi',  'dark', 'VANDELS', 1])  
#runs.append(['beast', 'dark', 'VANDELS', 1])  
#runs.append(['beast', 'dark', 'VANDELS', 1])

print('\n\nWelcome.\n\n')

sbatch  = False

##  input=$DESIHUB/desimodel/data/spectra/$SPECTRA.dat                                                                                                    
##  INPUT=/global/cscratch1/sd/mjwilson/desi/simspec/$TARGET-input-spectra-nir.fits                                                                       
##  
##  output=./dat/$SPECTRA-exp-$EXPOSURE.fits                                                                                                              
##  OUTPUT=/global/cscratch1/sd/mjwilson/desi/simspec/$TARGET-$SURVEY-spectra-exp-$EXPOSURE.fits                                                          

for [SURVEY, PROGRAM, TARGET, REPEAT] in runs:
 for nexp in np.arange(1, 9, 1):
  EXPOSURE = nexp * 1500. 
  EXPOSURE = '%.0lf' % EXPOSURE

  if SURVEY != 'pfs':
    INPUT='/global/cscratch1/sd/mjwilson/desi/runs/%s/%s-input-spectra.fits'                % (RUN, TARGET)
    OUTPUT='/global/cscratch1/sd/mjwilson/desi/runs/%s/%s-%s-spectra-exp-%s.fits'           % (RUN, TARGET, SURVEY, EXPOSURE)

  else:
    INPUT='/global/cscratch1/sd/mjwilson/desi/runs/%s/%s-input-spectra-nearir.fits'         % (RUN, TARGET)
    OUTPUT='/global/cscratch1/sd/mjwilson/desi/runs/%s/%s-%s-spectra-nearir-exp-%s.fits'    % (RUN, TARGET, SURVEY, EXPOSURE)

  if sbatch:
     args = (RUN, SURVEY, TARGET, PROGRAM, REPEAT, INPUT, OUTPUT, EXPOSURE, SURVEY + '-' + TARGET)
     cmd  = "sbatch  --export=All,RUN=%s,SURVEY=%s,TARGET=%s,PROGRAM=%s,REPEAT=%d,INPUT=%s,OUTPUT=%s, EXPOSURE=%s --job-name %s " % args 
     cmd += "quickspectra.sh"

  else: 
     args = (RUN, SURVEY, TARGET, PROGRAM, REPEAT, INPUT, OUTPUT, EXPOSURE)
     cmd  = 'export RUN=%s;export SURVEY=%s;export TARGET=%s;export PROGRAM=%s;export REPEAT=%d;export INPUT=%s;export OUTPUT=%s; export EXPOSURE=%s;' % args
     cmd += 'source ./quickspectra.sh'

  ## os.system('echo ' + cmd)
  os.system(cmd)

print('\n\nDone.\n\n')
