# this plotter is used for plotting a series of analysis
# created by AvrinJoker 2016
# e-mail:avrinjoker@hotmail.com

import os


from methods_plotter import *

from overall_best import *

from RMSD_formating import *

from IC50_formating import *


def target_analysis(target,dire):
    os.system('mkdir %s/%s/'%(dire,target))
    RMSD_formating('%s/'%target,'%s/%s/%s_RMSD.csv'%(dire,target,target))
    methods_plotter('%s/%s/%s_RMSD.csv'%(dire,target,target),'%s/figures/%s_RMSD.png'%(dire,target))
    #IC50_formating('%s/'%target,'%s/%s/%s_SPM.csv'%(dire,target,target),'%s/%s/%s_R2.csv'%(dire,target,target))
    #methods_plotter('%s/%s/%s_SPM.csv'%(dire,target,target),'%s/figures/%s_SPM.png'%(dire,target))
    #methods_plotter('%s/%s/%s_R2.csv'%(dire,target,target),'%s/figures/%s_R2.png'%(dire,target))


#os.system('mkdir rst/')
#os.system('mkdir rst/figures/')
#target_analysis('MAP4K4','./rst')
overall_best('example_methods_RMSD.csv','orange.png')

