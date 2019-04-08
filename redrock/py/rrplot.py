import os


rrplot = None

for p in os.environ['PATH'].split(':'):
    rrplot = os.path.join(p, 'rrplot')

    if os.path.exists(rrplot):
        break

if rrplot is None:
    print('ERROR: unable to find rrplot in your $PATH')

else:
    print('Using ' + rrplot)

run             = 'three'
target          =   'lbg'
survey          =   'pfs'
nexps           =      1
base_exp        =   1500.        ##  [Seconds].                                                                                                                                                                                       
repeat          =      1

exposure        =  nexps
exposure       *=  base_exp

if survey == 'pfs':
  nearir  = '-nearir'
  TEMPDIR = os.environ['BEAST'] + '/desihub/redrock/py/redrock/templates-nearir/'

else:
  nearir  = ''
  TEMPDIR = os.environ['BEAST'] + '/desihub/redrock/py/redrock/templates/'

##  
root            = '/global/cscratch1/sd/mjwilson/desi/runs/%s/' % run

##  Input files.
specfile        = root + '%s-%s-spectra%s-exp-%s.fits'  % (target, survey, nearir, '%.0lf' % exposure)
rrfile          = root + '%s-%s-spectra%s-exp-%s-rr.h5' % (target, survey, nearir, '%.0lf' % exposure) 

##  Now actually run it
os.system('%s --specfile  %s  --rrfile  %s  --templates  %s' % (rrplot, specfile, rrfile, TEMPDIR))
