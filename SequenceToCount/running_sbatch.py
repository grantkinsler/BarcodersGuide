import os

for sample in ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']:
    os.system(f'sbatch techRep_{sample}_demultiplex.sbatch')




##################################

import os

prefixes = ['070418_diagonal','080618_Nextseq','082818_Hiseq']
n_list = ['N' + str(n) for n in [716,718,719,720,721,722,723,724,726,727,728,729]] # nextera n primers
s_list = ['S' + str(s) for s in [513,515,516,517,518,520,521,522]] # nextera s primers


for prefix in prefixes:
    for n in n_list:
        for s in s_list:
            os.system(f'sbatch iHop_{prefix}_{n}{s}_demultiplex.sbatch')

