import matplotlib.pyplot as plt

def methods_plotter(infile,figsize=(4,3)):
    '''plotting methods comparison
    
    methods include (1)min_cross, (2)align_cross, (3)dock_cross, (4)align_close, (5)dock_close
    data includes RMSD, Spearman (IC50 & Kd) and R2 (IC50 & Kd)
    '''
    def _infile_parser(infile):
        lib = {'RMSD':'Percentage of RMSD <= 2.0 A','SPM_Kd':'Kd Spearman Correlation', 'SPM_IC50':'IC50 Spearman Correlation', 'R2_Kd':'Kd R2', 'R2_IC50':'IC50 R2'}
        l = [x.split(',') for x in open(infile,'r')]
        rec_name = l[0][1][:-1]
        measure_name = lib[l[1][1][:-1]]
        value = l[2:]
        x_value = [float(x[1]) for x in value]
        x_label = [x[0] for x in value]
        err_value = [float(x[2]) for x in value]
        return rec_name, measure_name, x_value, x_label, err_value

    rec_name, measure_name, x_value, x_label, err_value = _infile_parser(infile)

    # set figure basic info
    fig = plt.figure(figsize=figsize,dpi=300)
    ax = fig.add_subplot(111, ylim=(0,1))
    x = [x for x in range(len(x_value))]
    # add error bar
    ax.errorbar(x, x_value,yerr=err_value, ecolor='g', capthick=2, fmt=None, capsize=3, elinewidth=2, markeredgewidth=1)
    # plot the methods
    ax.scatter(x, x_value,facecolors='g', s=100, linewidths=0)
    # format the fig
    ax.set_ylabel(measure_name,fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(x_label)
    plt.xticks(rotation='vertical')
    plt.title(rec_name,fontweight='bold')
    #plt.show()
    plt.savefig('apple.png',dpi=300,bbox_inches='tight')