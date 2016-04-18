import os

from methods_plotter import *

from overall_formating import *

from overall_best import *


# plot for each receptor
for i in ['RMSD','IC50_SPM','IC50_R2']:
    os.system('mkdir %s/figure/'%i)
    l_names = [x.split('_')[0] for x in os.listdir('./%s/'%i) if x.endswith('csv')]
    for j in l_names:
        methods_plotter('./%s/%s_%s.csv'%(i,j,i),'./%s/figure/%s_%s.png'%(i,j,i))

# generate data for whole analysis
for i in ['RMSD','IC50_SPM','IC50_R2']:
    l_names = [x.split('_')[0] for x in os.listdir('./%s/'%i) if x.endswith('csv')]
    overall_formating(l_names,i,'./%s/figure/overall.csv'%i, dire='./%s/'%i)
    for j in ['min-cross','align-cross','dock-cross','align-close','dock-close']:
        overall_formating(l_names,i,'./%s/figure/%s.csv'%(i,j),False, j, dire='./%s/'%i)


for i in ['RMSD','IC50_SPM','IC50_R2']:
    overall_best('./%s/figure/overall.csv'%i, './%s/figure/overall.png'%i)
    for j in ['min-cross','align-cross','dock-cross','align-close','dock-close']:
        overall_best('./%s/figure/%s.csv'%(i,j), './%s/figure/%s.png'%(i,j))


for i in ['RMSD','IC50_SPM','IC50_R2']:
    all_in_one('%s_methods_all_compared.png'%i,'./%s/figure/'%i)
