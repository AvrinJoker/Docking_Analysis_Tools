import numpy as np

import os

import random

from scipy import stats

from scipy.stats import spearmanr as spm

from scipy.stats import linregress as lrg

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
        l_value = [x[1] for x in l]
        if len(l_value) == len(set(l_value)):
            l =[(x[1],x[0]) for x in l]
            l = sorted(l)
            l = [x[1] for x in l]
            return [l,]
        else:
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
    
    def _R2(l1,l2):
        '''calculate the r2 between two lists'''
        # check wether l1 and l2 are comparable
        if len(l1) != len(l2):
            print 'input lists have different length'
            return None
        # compute R2
        l1 = sorted(l1)
        l1 = [x[1] for x in l1]
        l2 = sorted(l2)
        l2 = [x[1] for x in l2]
        slope, intercept, r, p_value, std_err = lrg(l1,l2)
        return r**2


    def _commonPart(l1,l2):
        '''find the common part between two lists'''
        l2_names = [x[0] for x in l2]
        l_common = [x for x in l1 if x in l2_names]
        l1 = [x for x in l1 if x in l_common]
        l2 = [x for x in l2 if x[0] in l_common]
        return l1,l2

    
    def _Measure(l_lig,l_rec, dire, home_dire, func, close=False, ref='similarity.csv', exp='IC50.csv'):
        exp = '%s%s'%(home_dire,exp)
        l_rst = []
        
        l_exp = [x.split(',') for x in open(exp,'r')]
        l_exp = [(x[0],float(x[1])) for x in l_exp]

        # format l_exp and l_lig to get the common part
        l_lig,l_exp = _commonPart(l_lig,l_exp)
        
        # cross methods
        if close==False:
            for i in l_rec:
                l_small = []
                for j in l_lig:
                    small_v = [x.split(',') for x in open('%s%s/%s_%s.csv'%(dire,i,i,j),'r')][0][0]
                    l_small.append((j,float(small_v)))
                small_value = func(l_small,l_exp)
                l_rst.append(small_value)
        # close methods
        else:
            l_ref = [x.split(',') for x in open('%s%s'%(home_dire,ref),'r') if x!='']
            l_ref = [(x[0],x[1]) for x in l_ref]
            l_ref = dict(l_ref)
            l_small = []
            for i in l_lig:
                rec = l_ref[i][:-3]+'PRO'
                small_v =[x.split(',') for x in open('%s%s/%s_%s.csv'%(dire,rec,rec,i),'r')][0][0]
                l_small.append((i,float(small_v)))
            small_value = func(l_small,l_exp)
            l_rst.append(small_value)
        # best value
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
    
    
    def _iterst(iter, l_rec,l_lig,home_dire, func):
        '''iterate for many times'''
        # compute for 1000 iteration
        l_rst = {'min_cross':[],'align_cross':[],'dock_cross':[],'align_close':[],'dock_close':[]}
    
        counter=0
        for i in range(iter):
            # randomly reduce 10% of the list, minimal 1
            num = max(len(l_lig)/10,1)
            l_lig_r = random.sample(l_lig, len(l_lig)-num)
            # compute RMSD among different methods
            l_rst['min_cross'].append(_Measure(l_lig_r, l_rec, '%smin_cross/min_score/'%home_dire, home_dire, func))
            l_rst['align_cross'].append(_Measure(l_lig_r, l_rec, '%salign_cross/align_score/'%home_dire, home_dire, func))
            l_rst['dock_cross'].append(_Measure(l_lig_r, l_rec, '%sdock_cross/dock_score/'%home_dire, home_dire, func))
            l_rst['align_close'].append(_Measure(l_lig_r, l_rec, '%smin_cross/min_score/'%home_dire, home_dire, func, close=True))
            l_rst['dock_close'].append(_Measure(l_lig_r, l_rec, '%sdock_cross/dock_score/'%home_dire, home_dire, func, close=True))
            # set the counter to viualize process
            counter+=1
            if counter%100==0:
                print str(counter/10)+'%'
        return l_rst
    
    def _writeoutfile(outfile,method, l_rst):
        # write outfile
        f = open(outfile,'w')
        f.write('Receptor,%s\n'%home_dire[:-1])
        f.write('Measure,%s\n'%method)
        f.write('min-cross,%f,%f\n'%(l_rst['min_cross'][0],l_rst['min_cross'][1]))
        f.write('align-cross,%f,%f\n'%(l_rst['align_cross'][0],l_rst['align_cross'][1]))
        f.write('dock-cross,%f,%f\n'%(l_rst['dock_cross'][0],l_rst['dock_cross'][1]))
        f.write('align-close,%f,%f\n'%(l_rst['align_close'][0],l_rst['align_close'][1]))
        f.write('dock-close,%f,%f\n'%(l_rst['dock_close'][0],l_rst['dock_close'][1]))
        f.write('%s'%p)
        f.close()


    
    l_lig = [x[:-4] for x in os.listdir('%sligand/'%home_dire) if x.endswith('.sdf')]
    l_rec = [x[:-4] for x in os.listdir('%sreceptor/'%home_dire) if x.endswith('.pdb')]
    # generate IC50 summary
    l_rst_1 = _iterst(iter, l_rec,l_lig,home_dire, _SPM)
    l_rst_1,p = _generate_outfile(l_rst_1)
    _writeoutfile(outfile1,'SPM_IC50',l_rst_1)
    # generate R2 summary
    l_rst_2 = _iterst(iter, l_rec,l_lig,home_dire, _R2)
    l_rst_2,p = _generate_outfile(l_rst_2)
    _writeoutfile(outfile2,'R2_IC50',l_rst_2)






