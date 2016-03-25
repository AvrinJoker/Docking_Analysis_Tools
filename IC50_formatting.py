import numpy as np

import os

import random

from scipy import stats

from scipy.stats import spearmanr as spm

import itertools

def IC50_formating(home_dire,outfile1='IC50_SPM.csv',outfile2='IC50_R2.csv',iter=1000):
    '''calculate the Spearman Correlation and R2 between score and exp data across different methods
        
        output file includes calculation results for different methods, which shows mean
        and stderr by randomly reduced 10% of the structures over 1000 iteration. p-value
        is also calculated.
        
        e.g. of output file:
        ------------------------------
        Receptor,MAP4K4
        Measure,RMSD
        min-cross,0.499229,0.023409
        align-cross,0.499714,0.023593
        dock-cross,0.395514,0.023077
        align-close,0.448143,0.023901
        dock-close,0.263457,0.020882
        p-value,0.000000,1,3
        ------------------------------
    '''
    def _genpair(l):
        '''generate rankings (small to large) for all combination'''
        lib = {}
        # histogram the list
        for i in l:
            if i[1] in lib:
                lib[i[1]].append(i[0])
            else:
                lib[i[1]] = [i[0],]
        # sort from samll to large
        l_rst = [(x,lib[x]) for x in lib]
        l_rst = sorted(l_rst)
        # generate all ranking combinations
        l_rst = [x[1] for x in l_rst]
        l_rst = [_genComb(x) for x in l_rst]
        l_rst = list(itertools.product(*l_rst))
        l_rst = [''.join(x) for x in l_rst]
        l_rst = [list(x) for x in l_rst]
        return l_rst

    def _genComb(l):
    '''generate all combinations of ranking for same value ones'''
        if len(l) == 1:
            return l[0]
        else:
            l = list(itertools.permutations(l,len(l)))
            l = [''.join(x) for x in l]
            return l

    def _name2rank(l1,l2):
        '''transfer the name ranking to the number ranking'''
        l1 = [(l1[x],x) for x in range(len(l1))]
        lib = dict(l1)
        l2 = [lib[x] for x in l2]
        l1 = [x[1] for x in l1]
        return l1,l2


    def _SPM(l1,l2):
    '''calculate the spearman correlation between two lists'''
        # check wether l1 and l2 are comparable
        if len(l1) != len(l2):
            print 'input lists have different length'
            return None
        # generate pairs for both list:
        l1 = _genpair(l1)
        l2 = _genpair(l2)
        #compute SPM
        l_rst = []
        for i in l1:
            for j in l2:
                a,b = _name2rank(i,j)
                SPM_value = spm(a,b)
                l_rst.append(SPM_value[0])
        return max(l_rst)

    
    def _IC50(l_lig,l_rec, dire, home_dire, close=False, ref='similarity.csv'):
        l_rst = []
        if close==False:
            for i in l_rec:
                num = 0
                for j in l_lig:
                    RMSD = [x.split(',') for x in open('%s%s/%s_%s.csv'%(dire,i,i,j),'r')][0][1]
                    if float(RMSD) <= 2.0:
                        num+=1
                l_rst.append(num/float(len(l_lig)))
        else:
            l_ref = [x.split(',') for x in open('%s%s'%(home_dire,ref),'r') if x!='']
            l_ref = [(x[0],x[1]) for x in l_ref]
            l_ref = dict(l_ref)
            num = 0
            for i in l_lig:
                rec = l_ref[i][:-3]+'PRO'
                RMSD =[x.split(',') for x in open('%s%s/%s_%s.csv'%(dire,rec,rec,i),'r')][0][1]
                if float(RMSD) <= 2.0:
                    num+=1
            l_rst.append(num/float(len(l_lig)))
        
        return max(l_rst)
    
    def _generate_outfile(lib):
        lib_new = {}
        for i in lib:
            l = lib[i]
            l = sorted(l, reverse=True)
            mean = np.mean(l)
            sterr = np.std(l)
            lib_new[i]=[mean,sterr]
        
        # label the first p
        l_rank = [(lib_new['min_cross'][0],0),(lib_new['align_cross'][0],1),(lib_new['dock_cross'][0],2),(lib_new['align_close'][0],3),(lib_new['dock_close'][0],4)]
        l_rank = sorted(l_rank,reverse=True)
        l_rank = [x[1] for x in l_rank]
        l_rank = [(l_rank[0],x) for x in l_rank[1:]]
        
        # t-test for all pairs
        l_tt = [(0,lib['min_cross']),(1,lib['align_cross']),(2,lib['dock_cross']),(3,lib['align_close']),(4,lib['dock_close'])]
        l_tt_r = []
        for i,j in l_rank:
            t,p = stats.ttest_ind(l_tt[i][1],l_tt[j][1])
            if p <= 0.05:
                l_tt_r.append('p-value,%f,%d,%d'%(p,i,j))
        return lib_new, l_tt_r[0]
    
    
    # get names for ligands and receptor
    l_lig = [x[:-4] for x in os.listdir('%sligand/'%home_dire) if x.endswith('.sdf')]
    l_rec = [x[:-4] for x in os.listdir('%sreceptor/'%home_dire) if x.endswith('.pdb')]
    
    # compute RMSD percentage for 1000 iteration
    l_rst = {'min_cross':[],'align_cross':[],'dock_cross':[],'align_close':[],'dock_close':[]}
    
    counter=0
    for i in range(1000):
        # randomly reduce 10% of the list, minimal 1
        num = max(len(l_lig)/10,1)
        l_lig_r = random.sample(l_lig, len(l_lig)-num)
        # compute RMSD among different methods
        l_rst['min_cross'].append(_RMSD(l_lig_r, l_rec, '%smin_cross/min_score/'%home_dire, home_dire))
        l_rst['align_cross'].append(_RMSD(l_lig_r, l_rec, '%salign_cross/align_score/'%home_dire, home_dire))
        l_rst['dock_cross'].append(_RMSD(l_lig_r, l_rec, '%sdock_cross/dock_score/'%home_dire, home_dire))
        l_rst['align_close'].append(_RMSD(l_lig_r, l_rec, '%smin_cross/min_score/'%home_dire, home_dire, close=True))
        l_rst['dock_close'].append(_RMSD(l_lig_r, l_rec, '%sdock_cross/dock_score/'%home_dire, home_dire, close=True))
        # set the counter to viualize process
        counter+=1
        if counter%100==0:
            print str(counter/10)+'%'

    #
    l_rst,p = _generate_outfile(l_rst)
    
    # write outfile
    f = open(outfile,'w')
    f.write('Receptor,%s\n'%home_dire[:-1])
    f.write('Measure,IC50\n')
    f.write('min-cross,%f,%f\n'%(l_rst['min_cross'][0],l_rst['min_cross'][1]))
    f.write('align-cross,%f,%f\n'%(l_rst['align_cross'][0],l_rst['align_cross'][1]))
    f.write('dock-cross,%f,%f\n'%(l_rst['dock_cross'][0],l_rst['dock_cross'][1]))
    f.write('align-close,%f,%f\n'%(l_rst['align_close'][0],l_rst['align_close'][1]))
    f.write('dock-close,%f,%f\n'%(l_rst['dock_close'][0],l_rst['dock_close'][1]))
    f.write('%s'%p)
    f.close()
