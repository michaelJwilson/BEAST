import  os 


RUN   =   'four'
runs  =      []

##  survey, program, target, repeat
#runs.append(['desi',  'dark',   'lbg'])
#runs.append(['beast', 'dark',   'lbg'])
runs.append(['pfs', 'dark',   'lbg'])
#runs.append(['desi',  'dark',   'BC03'])                                                                                                       
#runs.append(['beast', 'dark',   'BC03']) 
#runs.append(['desi',  'dark',   'elg'])
#runs.append(['beast', 'dark',   'elg'])
#runs.append(['desi',  'dark',   'qso'])
#runs.append(['beast', 'dark',   'qso'])
#runs.append(['desi',  'dark',   'lrg'])                                                                                                      
#runs.append(['beast', 'dark',   'lrg'])
#runs.append(['desi',  'bright', 'bgs'])                                                                                                         
#runs.append(['beast', 'bright', 'bgs']) 
#runs.append(['beast', 'dark', 'VANDELS'])

store  = os.environ['NERSC_HOST']
sbatch = True

print('\n\nWelcome.\n\n')

for [SURVEY, PROGRAM, TARGET] in runs:
    if sbatch:
      os.system("unset NERSC_HOST")
       
      ##  --export=ALL; env -i
      cmd = 'sbatch  --job-name %s -C haswell --export=ALL,SURVEY=%s,TARGET=%s,SBATCH=%s,RUN=%s redrock.sh' % (SURVEY + '-' + TARGET, SURVEY,\
                                                                                                               TARGET, sbatch, RUN)  
            
    else:
      cmd = 'export SURVEY=%s; export TARGET=%s; export SBATCH=%s; export RUN=%s; source ./redrock.sh' % (SURVEY, TARGET, sbatch, RUN)  

    ##  os.system('echo ' + cmd)
    os.system(cmd)

os.system("export NERSC_HOST=%s" % store)

print('\n\nDone.\n\n')
