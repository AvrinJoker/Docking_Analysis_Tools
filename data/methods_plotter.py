import matplotlib.pyplot as plt

def methods_plotter(infile,outfile='apple.png',figsize=(4,3)):
    '''plotting methods comparison
    
    methods include (1)min_cross, (2)align_cross, (3)dock_cross, (4)align_close, (5)dock_close
    data includes RMSD, Spearman (IC50 & Kd) and R2 (IC50 & Kd)
    '''
    def _infile_parser(infile):
        lib = {'RMSD':'Percentage of RMSD <= 2.0 A','Kd_SPM':'Kd Spearman Correlation', 'IC50_SPM':'IC50 Spearman Correlation', 'Kd_R2':'Kd R2', 'IC50_R2':'IC50 R2'}
        l = [x.split(',') for x in open(infile,'r')]
        rec_name = l[0][1][:-1]
        rec_name = rec_name.split('/')[-1]
        measure_name = lib[l[1][1][:-1]]
        value = l[2:-1]
        x_value = [float(x[1]) for x in value]
        x_label = [x[0] for x in value]
        err_value = [float(x[2]) for x in value]
        p = l[-1]
        p = [float(p[1]), int(p[2]),int(p[3])]
        return rec_name, measure_name, x_value, x_label, err_value, p
    
    # get plotting data from file
    rec_name, measure_name, x_value, x_label, err_value, p = _infile_parser(infile)

    # set figure basic info
    fig = plt.figure(figsize=figsize,dpi=600)
    ax = fig.add_subplot(111, ylim=(0,1), xlim=(-1,len(x_value)))
    x = [x for x in range(len(x_value))]
    
    # add error bar & plot the methods
    ax.errorbar(x, x_value,yerr=err_value, ecolor='g', capthick=2, fmt=None, capsize=3, elinewidth=2, markeredgewidth=1)
    ax.scatter(x, x_value,facecolors='g', s=100, linewidths=0, alpha=0.8)
    
    # set the axises & format the fig
    plt.title(rec_name,fontweight='bold')
    ax.set_ylabel(measure_name,fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(x_label)
    plt.xticks(rotation='vertical')
    
    # mark the x value
    x_value = [float('%.3f'%x) for x in x_value]
    for i in range(len(x_value)):
        ax.annotate(x_value[i], xy=(i,x_value[i]-err_value[i]-0.07),ha='center',fontsize=8, color='g')
    
    # draw the p value label & connection line
    if p[0] < 0.001:
        p_label = 'p < 0.001'
    else:
        p_label = 'p = %.4f'%p[0]
    fraction = 0.4*(p[2]-p[1])/abs(p[1]-p[2])
    ax.annotate(p_label, xy=(0.5*p[1]+0.5*p[2],0.14+max(x_value[p[1]]+err_value[p[1]],x_value[p[2]]+err_value[p[2]])), ha='center', fontsize=8, alpha=0.8)
    ax.annotate('',xy=(p[1],err_value[p[1]]+x_value[p[1]]), xytext=(p[2],err_value[p[2]]+x_value[p[2]]), arrowprops=dict(arrowstyle='-',connectionstyle='bar,angle=-180,fraction=%f'%fraction))
    # save figure
    #plt.show()
    plt.savefig(outfile,dpi=600,bbox_inches='tight')