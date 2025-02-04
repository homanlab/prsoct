"""
@author: Finn Rabe
@email: finn.rabe@bli.uzh.ch
@terms: CC-BY-NC-ND
"""

# set working dir
import os
#os.chdir('/retpsy/src')
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# importing  all the functions defined in retpsy_load.py
from retinflam_load import *


def load_data_type(dis, disorder, format_type, meas_gr):
    """load data depending on different demands
        dis: define disorder abbreviation, e.g. SZ,BD,PD
        disorder: ex-/include diagnosed participants, choose 'diagnosed', 'none' or 'all'
        format_type: format oct data either to 'wide' or 'long' format (is needed for mixed linear model)
        meas_gr: select if you want to analyse general retinal phenotypes or macular subfields
    """
    # load and format data
    df_data = general_formatting(dis,disorder,format_type,meas_gr)
    y_meas,cov,cp = sel_oct_group(meas_gr, format_type)
    x_meas = ['PRS'+dis]
    df_data_trans, y_meas_lat = oct_formatting(df_data,meas_gr,y_meas,format_type, disorder)
    #df_data_trans.iloc[0].to_csv('/Users/frabe/Desktop/col.csv')
    return df_data_trans,y_meas,y_meas_lat,x_meas,cov,cp
    
def parreg_oct(df_zcorr,y_meas,x_meas,cov,figname):
    # Individualal retina measures partial regression

    cpfig5 = sns.color_palette("Paired")[2:]
    ax_dim = [1,len(y_meas)]
    fig, subfig = plt.subplots(ax_dim[0], ax_dim[1], figsize=fig_size, sharey=True,layout='constrained')
    #plt.title('Retinal thickness '+'('+side+')')
    #fig.suptitle('Retinal thickness '+'('+side+')')
    if x_meas[0] == 'PRSSZ':
        xlbl = 'Polygenic risk score for schizophrenia (z)'
    
    if ax_dim[1] > 1:
        if meas_gr == 'overall':
            side = ['LEFT EYE','RIGHT EYE']
            ylbl = 'Macular thickness (µm)'
        else:
            side = ['INNER RETINA','OUTER RETINA']
            ylbl = 'Thickness (µm)'
    else:
        if y_meas[0] == 'Macula_mean':
            side = 'Overall Macula'
        elif y_meas[0] == 'GCIPL_mean':
            side = 'GC-IPL'
        elif y_meas[0] == 'INL_RPE_mean':
            side = 'Outer retina (INL-RPE)'
        ylbl = 'Mean Thickness (µm)'
    fig.supxlabel(xlbl)
    fig.supylabel(ylbl)
    count = 0 
    # Regress out nuissance regressors 
    #xarr = regress_out(df_corr[pc_lbl[yv]], df_corr[cov])
    xarr = regress_out(df_zcorr[x_meas[0]], df_zcorr[cov])

    # Compute partial regression incl. covariates as nuissance reg
    pstats = pg.partial_corr(data=df_zcorr, y=y_meas[count], x=x_meas[0], method='pearson')#,covar=cov)
    r_partial = np.round(pstats.r.values[0],2)
    p_partial = np.round(pstats['p-val'].values[0],3)

    # compute bivariate correlation
    #r,p = r2(df_corr,x_meas[0],y_meas[count])
    #print(r)

    if ax_dim[1] > 1:
        for xv in range(ax_dim[0]):
            for yv in range(ax_dim[1]):
                axid = subfig[yv]
                g = sns.regplot(ax=axid, y=df_zcorr[y_meas[count]], x=xarr, scatter=False, line_kws={'linewidth':line_size}, color=cpfig5[yv], ci=95, robust=True)
                if meas_gr == 'overall':
                    xlm = [-4.5,4.5]
                    ylm = [275,281]
                else:
                    xlm = [-3.5,3.5]
                    ylm = [133,143.5]
                # set subplot title 
                g.title.set_text(side[count])
                g.set(xlabel=None)
                g.set(ylabel=None)
                g.set(xlim=xlm)
                g.set(ylim=ylm)

                #xticks = ['{:,.2f}'.format(x) for x in g.get_xticks()]
                #yticks = ['{:,.2f}'.format(x) for x in g.get_yticks()]

                #g.set_xticklabels(xticks, fontsize=tick_size)
                #g.set_yticklabels(yticks,fontsize=tick_size-1)
                count += 1
                sns.despine()
    else:
        g = sns.regplot(ax=subfig, y=df_zcorr[y_meas[count]], x=xarr, scatter=False, line_kws={'linewidth':line_size}, color=cpfig5[0], ci=95, robust=True)
        xlm = [-4.5,4.5]
        if y_meas[0] == 'Macula_mean':
            ylm = [276,280]
        elif y_meas[0] == 'GCIPL_mean':
            ylm = [73,74]
        elif y_meas[0] == 'INL_RPE_mean':
            ylm = [141,143]
            yticks = [round(x) for x in g.get_yticks()]
            g.set_yticklabels(yticks)
        # set subplot title 
        g.title.set_text(side)
        g.set(xlabel=None)
        g.set(ylabel=None)
        g.set(xlim=xlm)
        g.set(ylim=ylm)

        #xticks = ['{:,.2f}'.format(x) for x in g.get_xticks()]
        #yticks = ['{:,.2f}'.format(x) for x in g.get_yticks()]

        #g.set_xticklabels(xticks, fontsize=tick_size)
        #g.set_yticklabels(yticks,fontsize=tick_size-1)
        sns.despine()

    plt.savefig('../output/figures/'+figname+'_'+meas_gr+'_retinap_partialreg.png', format='png', bbox_inches='tight')

def parreg_pc(df_zcorr,pc_lbl,x_meas,cov,figname='Fig3'):
    # Partial regression with first principle component 
    figlbl = ['a','b']
    ax_dim = [1,1]
    pc_len = 1

    count = 0 
    for yv in range(pc_len):
        fig, subfig = plt.subplots(ax_dim[0], ax_dim[1], sharex=True, figsize=fig_size)
        fig.tight_layout() 
        # Regress out nuissance regressors 
        #xarr = regress_out(df_corr[pc_lbl[yv]], df_corr[cov])
        xarr = regress_out(df_zcorr[x_meas[0]], df_zcorr[cov])

        # Compute partial regression incl. covariates as nuissance reg
        pstats = pg.partial_corr(data=df_zcorr, y=pc_lbl[yv], x=x_meas[0], covar=cov, method='pearson')
        r_partial = np.round(pstats.r.values[0],2)
        p_partial = np.round(pstats['p-val'].values[0],3)
        df_results[pc_lbl[yv]+'partialr'] = [r_partial]
        df_results[pc_lbl[yv]+'partialp'] = [p_partial]

        axid = subfig

        g = sns.regplot(ax=axid, y=xarr, x=df_zcorr[pc_lbl[yv]], scatter_kws={"color": "white"}, line_kws={'linewidth':line_size}, color=cp[yv], ci=95, robust=True)
        #add regression equation to plot
        #axx = g.axes[xv,yv]
        #g.text(0, 0.1,'r\u00b2= {}'.format(r_partial) +pval_asterisks(p_partial,0.05), fontsize=val_size)
        # g.text(0, 0.1,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
        # g.set(ylim=(-0.3,0.3))
        # #g.xaxis.set_ticks(g.get_xticks())
        # #g.yaxis.set_ticks(ylabels)
        # xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
        # print('xticks:',xlabels)
        # g.set_xticklabels(xlabels, fontsize=tick_size)
        # ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
        # g.set_yticklabels(ylabels,fontsize=tick_size)

        
        if figname == 'Fig1':
            xlbl = 'PRS SZ (z)'
            g.text(0, 0.05,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.1,0.1))
            g.set(xlim=(-20,20))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        elif figname == 'Fig3A':
            xlbl = 'Polygenic risk score for neuroinflammatory pathway'
            #g.text(0, 0.05,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.1,0.1))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            ylbl = 'C-reactive protein log(mg/L)'
            count = 1
            g.set_ylabel(ylbl,fontsize=label_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        elif figname == 'Fig3B':
            xlbl = 'Polygenic risk score for WnT signaling pathway'
            g.text(0, 0.05,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.1,0.1))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        elif figname == 'Fig3C':
            xlbl = 'C-reactive protein log(mg/L)'
            g.text(0, 0.4,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.2,0.7))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        if count == 0:
            ylbl = 'Principal component 1'
            g.set_ylabel(ylbl,fontsize=label_size)
            #g.set_ylabel(g.get_ylabel(),fontsize=label_size)
        #else:
        #    g.set_ylabel('',fontsize=label_size)
        count += 1
        sns.despine()

        plt.savefig('../output/figures/'+figname+figlbl[yv]+'_'+meas_gr+'_'+pc_lbl[yv]+'_partialreg.png', format='png', bbox_inches='tight')

def mixed_lm(df_corr, y_var):
    cov_ed = '+'.join(cov)
    mlmodel = sm.MixedLM.from_formula(y_var+' ~ '+x_meas[0]+'+'+cov_ed, groups=df_corr['Participant ID'],data=df_corr)
    rel_mlm = mlmodel.fit(method=["lbfgs"])
    #method optimisation: 
    #https://docs.scipy.org/doc/scipy/tutorial/optimize.html#broyden-fletcher-goldfarb-shanno-algorithm-method-bfgs
    #rel_mlm_table = rel_mlm.summary().tables[1]
    #rel_mlm_table.to_csv('../output/figures/'+figname+'_'+meas_gr+'_'+y_var+'_mixedlm_results.csv')

    # create dataframe for summary statistics table
    #df_mlmresults = _to_dataframe(rel_mlm, y_var)
    return rel_mlm

def pca_loadscores(df_corr, y_meas, pc_lbl, pca, figname = 'Fig2B'):
    # weights are nothing but the eigenvectors of X
    y_df = df_corr[y_meas]
    fig, g = plt.subplots(figsize=fig_size)

    #cp = sns.color_palette("colorblind")
    cp = sns.color_palette("Paired")
    pca_weights = pca.components_.flatten()

    df_pca_weights = pd.DataFrame(data={'Principal Component': np.repeat(pc_lbl,len(y_df.columns)),\
                        'Variables': np.tile(y_df.columns,len(pc_lbl)),\
                        'PCA loading scores': pca_weights})

    for lc in range(len(pca_weights)):
        pcn = df_pca_weights['Principal Component'][lc]
        retp = df_pca_weights['Variables'][lc]
        ldc = df_pca_weights['PCA loading scores'][lc]
        df_results[meas_gr+'_pca'+pcn+'_'+retp+'_loadingscore'] = [np.round(ldc,2)]
    #sns.set(rc={'figure.figsize':(10,8)})
    g = sns.barplot(x=df_pca_weights.columns[0], y=df_pca_weights.columns[2], hue=df_pca_weights.columns[1], data=df_pca_weights, palette=cp)
    #axes.legend(loc='lower left', bbox_to_anchor=(1, 0.5),fontsize = tick_size)
    if meas_gr == 'overall':
        ti = 'Retinal phenotypes'
    else:
        ti = 'Macular bilateral subfields'
    sns.move_legend(g,'lower center', bbox_to_anchor=(.5, 1), ncol=3, fontsize = tick_size, title=ti, frameon=False)
    sns.despine()
    fig = g.get_figure()
    fig.savefig('../output/figures/'+figname+'_'+meas_gr+'_pca_loadingscores.png', format='png', bbox_inches='tight')
    return df_pca_weights

def pca_cumvar(df_corr,y_meas,var_cutoff=0.8, figname = 'Fig2A'):
    # Compute Principle component analysis
    
    y_df = df_corr[y_meas]
    pca = PCA(n_components=var_cutoff)
    pcs = pca.fit_transform(y_df)

    pc_lbl = ["PC%01d" %i for i in range(1,pca.n_components_+1)]
    pca = PCA(n_components=len(pc_lbl))

    # weight oct measures
    pcs = pca.fit_transform(y_df)

    # cumulative variance 
    exp_var = pca.explained_variance_ratio_
    cum_sum_eigenvalues = np.cumsum(pca.explained_variance_ratio_)
    print('shape before pcs:',df_corr.shape)
    for pc in range(len(pc_lbl)):
        df_corr[pc_lbl[pc]] = pcs[:,pc]
    
    # save to results df
    if figname == 'Fig2A':
        df_results[meas_gr+'_pca_ncomponents'] = [pca.n_components_]
        for ev in range(len(exp_var)):
            print('varr:',np.round(exp_var[ev],2))
            df_results[meas_gr+'_pca'+str(ev+1)+'_expvariance'] = [np.round(exp_var[ev],2)]
    
    # Display Scree plot
    df_pca = pd.DataFrame(data={'Principal Component': pc_lbl,\
                        'Variance explained': pca.explained_variance_ratio_,\
                        'Cumulative Variance': cum_sum_eigenvalues
                        })
    plt.figure(figsize=fig_size)
    #sns.barplot(x=df_pca.columns[0], y=df_pca.columns[1], data=df_pca, palette=cp[:len(pc_lbl)])
    #plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')

    sns.barplot(x=df_pca.columns[0], y=df_pca.columns[2], data=df_pca, palette=cp[:len(pc_lbl)])
    sns.lineplot(x=df_pca.columns[0], y=df_pca.columns[2], data=df_pca, marker='o', sort = False)
    plt.axhline(var_cutoff, ls='--',color=thr_line)
    #plt.text(0,0.9, '90% threshold')
    plt.ylim((0,1))
    sns.despine()
    print('shape after pcs:',df_corr.shape)

    plt.savefig('../output/figures/'+figname+'_'+meas_gr+'_pca.png', format='png',bbox_inches='tight')
    return pc_lbl, pca, df_pca, df_corr, exp_var, cum_sum_eigenvalues

def corr_matrix(df_corr,y_meas,figname='Appx_Fig2'):
    covm_all = df_corr[y_meas].corr()
    ma = covm_all.round(2)

    fig, g = plt.subplots(figsize=fig_size)
    cp = sns.color_palette("YlOrRd", as_cmap=True)
    g = sns.heatmap(ma, annot=False,cmap=cp)
    #g = sns.heatmap(ma, annot=True,annot_kws={"fontsize":val_size},cmap=cp)

    yl = [t.get_text()  for t in g.get_yticklabels()]
    y_lbl = [s.strip('_') for s in yl]
    #y_lbl = ['GC-IPL left','GC-IPL right','Macula left','Macula right','RNFL left','RNFL right','INL left','INL right']#,'RPE left','RPE right']

    g.set_xticklabels(y_lbl, fontsize=tick_size)
    g.set_yticklabels(y_lbl, fontsize=tick_size)

    fig = g.get_figure()
    fig.savefig('../output/figures/'+figname+'_'+meas_gr+'_corrmatrix.png', format='png', bbox_inches='tight')

def population(raw,filt_hc,filt_scz,figname='Appx_Table1'):
    # number of participants
    n_total = len(raw)
    n_filt = len(filt_hc)
    n_scz = len(filt_scz)

    # sex ratio
    sex_vcounts = filt_hc.Sex.value_counts()

    # Compute female, male ratio
    df_sex = pd.DataFrame()
    df_sex['Sex'] = ['Male','Female']
    df_sex['N'] = [sex_vcounts[1],sex_vcounts[2]]

    t,p = stats.ttest_ind(filt_hc.Sex==1, filt_hc.Sex==2)
    p = pval_asterisks(p, alpha=0.05)

    fig =  ff.create_table(df_sex)
    fig.update_layout(
        autosize=False,
        width=700,
        height=400,
    )
    fig.write_image('../output/figures/'+figname+'_'+meas_gr+'_sex_ratio.png', scale=2)
    return n_total, n_filt, n_scz, sex_vcounts, t, p

def collinarity_dist(df_corr, y_meas, figname = 'Appx_Fig1'):
    ## Compute significance of collinarity of OCT measures
    df_coll = df_corr[y_meas].rcorr(padjust = 'bonf')

    #replace nan values with empty string
    coll_val = df_coll.values.flatten()
    coll_val[coll_val == '-'] = ''
    #coll_val = coll_val + ' $^3$'

    # hide upper triangle 
    def hide_current_axis(*args, **kwds):
        plt.gca().set_visible(False)

    fig, g = plt.subplots(figsize=fig_size)

    g = sns.pairplot(df_corr[y_meas], kind="reg", corner=False,\
                        #palette=scatter_plots,
                        markers="o", \
                        diag_kws = {'color': scatter_plots},\
                        plot_kws={'line_kws':{'color':reg_line}, 'scatter_kws': {'alpha': 0.5, 'color': scatter_plots}})
    
    # iterate through each axes and insert collinerity value
    count = 0
    for xpl in range(len(y_meas)):
        for ypl in range(len(y_meas)):
            #add regression equation to plot
            axx = g.axes[xpl,ypl]

            if meas_gr == 'overall':
                #axx.set(ylim=(-20, 50))
                #axx.set(xlim=(0, 800))
                if coll_val[count] == '':
                    axx.text(15, 20, str(coll_val[count]), size=16, color='black')
                else:
                    axx.text(15, 20, str(coll_val[count]), size=16, color='black') #+r'$^{***}$'
            else:
                #axx.set(ylim=(0, 1000))
                #axx.set(xlim=(0, 1000))
                if coll_val[count] == '':
                    axx.text(0, 10, str(coll_val[count]), size=23, color='black')
                else:
                    #axx.text(0, 10, str(coll_val[count])+r'$^{***}$', size=16, color='black')
                    axx.text(0, 10, str(coll_val[count]), size=23, color='black')
            count += 1

    g.map_upper(hide_current_axis)
    plt.savefig('../output/figures/'+figname+'_'+meas_gr+'_regmatrix.png', format='png')

def vif(df_corr,y_meas,figname = 'Appx_Table2'):
    ## Compute variance inflation factor
    # VIF dataframe
    df_sex = pd.DataFrame()
    df_sex['Phenotype'] = y_meas
    
    # calculating VIF for each feature
    df_sex['VIF'] = [variance_inflation_factor(df_corr[y_meas].values, i)
                            for i in range(len(y_meas))]
    df_sex['VIF'] = df_sex['VIF'].round(2)
    
    fig =  ff.create_table(df_sex)
    fig.update_layout(
        autosize=False,
        width=700,
        height=400,
    )
    fig.write_image('../output/figures/'+figname+'_'+meas_gr+'_vif.png', scale=2)

def save_pcs(df_dlm_path,df_res,y_meas_lat):
    df_dlm = pd.read_csv(df_dlm_path+'.csv')
    #y_meas_lat = ['GC_IPL_left','GC_IPL_right','Macula_left','Macula_right','RNFL_left','RNFL_right']

    y_df = df_res[y_meas_lat]
    pca = PCA(n_components=0.7)
    pcs = pca.fit_transform(y_df)
    pc_lbl = ["PC%01d" %i for i in range(1,pca.n_components_+1)]
    for pc in range(len(pc_lbl)):
            df_dlm[pc_lbl[pc]] = pcs[:,pc]
    df_dlm.to_csv(df_dlm_path+'_pc.csv')
    return df_dlm

def count_dropped_rows(df, columns):
    n_before = len(df)
    
    # Create a mask for NaN and infinite values
    mask = df[columns].isin([np.nan, np.inf, -np.inf])
    
    # Replace infinite values with NaN
    df = df.replace([np.inf, -np.inf], np.nan)
    
    df_dropped = df.dropna(subset=columns)
    n_after = len(df_dropped)
    
    # Find columns with at least one missing or infinite value
    missing_cols =  df[columns].columns[df[columns].isnull().any()].tolist()
    df = df[missing_cols]
    
    # Count NaN and infinite values in each column
    dropped_counts = mask[missing_cols].sum().sort_values(ascending=False)
    dropped_counts['Total'] = n_before - n_after
    dropped_counts.index = dropped_counts.index.str.replace('_', ' ')
    
    return df_dropped[columns], dropped_counts.to_frame(name='Missing values')

def check_mcar(df):
    # Initialize the MCAR Test
    mcar_test = MCARTest(method="little")
    # Perform Little's MCAR test
    p_value = mcar_test.little_mcar_test(df)
    if p_value < 0.05:
        result = "Data is NOT MCAR"
        print(result)
    else:
        result = "Data is MCAR"
        print(result)
    return p_value

if __name__ == "__main__":
    # set figure layout
    fig_size, tick_size, label_size, title_size, val_size, line_size, dist_plots, scatter_plots, reg_line, thr_line = set_fig_style()
    
    #eye indicies 
    side = ['left','right'] 

    # Get macular thickness values
    meas_gr = 'overall'
    df_data_trans, y_meas, y_meas_lat, x_meas, cov, cp = load_data_type('SZ', 'none', 'wide', meas_gr)
    n_nondiag = len(df_data_trans)
    # Average across eyes
    df_data_trans['Macula_mean'] = df_data_trans[y_meas_lat].mean(axis=1)
    df_data_trans['Image_quality_mean'] = df_data_trans[['OCT_quality_left','OCT_quality_right']].mean(axis=1)
    # Add covariates to df
    cova = ['Age', 'Age_squared', 'Sex', 'Smoking_status', 'Alcohol_drinker_status', 'BMI', 'Hypertension','Diabetes_mellitus',
           'Genotype_array', 'Townsend_index', 
           'Genetic PC1', 'Genetic PC2', 'Genetic PC3', 'Genetic PC4', 'Genetic PC5', 
           'Genetic PC6', 'Genetic PC7', 'Genetic PC8', 'Genetic PC9', 'Genetic PC10', 'Image_quality_mean', 'OCT_quality_left', 'OCT_quality_right']
    prs_thres = ['PRSSZ_0_001', 'PRSSZ_0_05', 'PRSSZ_0_1', 'PRSSZ_0_2', 'PRSSZ_0_3', 'PRSSZ_0_4', 'PRSSZ_0_5', 'PRSSZ_1']
    col_incl = y_meas_lat+x_meas+cova+['Participant ID']+['Macula_mean']+prs_thres
    print('columns included:',col_incl)

    df_macula_miss = df_data_trans[col_incl]
    df_macula_miss.to_csv('../output/data/df_mlm_macula_miss.csv')
    df_data_corr, missing_vals = count_dropped_rows(df_data_trans, col_incl)
    print('len df for macula analysis:', len(df_data_corr))
    df_data_corr.to_csv('../output/data/df_mlm_macula.csv')
    
    # Plot partial regression between PRS and average macular thickness and left/right
    parreg_oct(df_data_corr,['Macula_mean'],x_meas,cova,figname='Fig1a')
    parreg_oct(df_data_corr,y_meas_lat,x_meas,cova,figname='Appx_FigS1')

    n_nondiag_nan = len(df_data_trans)-len(df_data_corr)
    print('n_nondiag_nan:',n_nondiag_nan)
    missing_vals.to_csv('../output/data/df_missing_vals.csv')

    # Create a new DataFrame including only the columns with missing values
    mis_col_inc = x_meas+cova+['Participant ID']+['Macula_mean']
    df_colfilt = df_data_trans[mis_col_inc]
    missing_col = df_colfilt.columns[df_colfilt.isna().any()].tolist()
    print('col missing data:', missing_col)
    mcar_pval = check_mcar(df_data_trans[missing_col])
    print("Missing data assessment:")
    print(mcar_pval)

    # Test for normal distributed phenotypes
    df_norm = pd.DataFrame({'Phenotype':[],'Statistic':[],'p':[]})
    for retsub in range(len(y_meas_lat)):
        res = stats.normaltest(df_data_corr[y_meas_lat[retsub]])
        new_row = pd.DataFrame({
            'Phenotype': [y_meas_lat[retsub]],
            'Statistic': [res.statistic],
            'p': [pval_asterisks(res.pvalue, alpha=0.05)]
        })
        df_norm = pd.concat([df_norm, new_row], ignore_index=True)
    df_norm.to_csv('../output/data/df_normtest_macula.csv')

    curr_smoker = df_data_corr.Smoking_status[df_data_corr.Smoking_status == 2].count()
    prev_smoker = df_data_corr.Smoking_status[df_data_corr.Smoking_status == 1].count()
    non_smoker = df_data_corr.Smoking_status[df_data_corr.Smoking_status == 0].count()
    curr_alcohol_drinker = df_data_corr.Alcohol_drinker_status[df_data_corr.Alcohol_drinker_status == 2].count()
    prev_alcohol_drinker = df_data_corr.Alcohol_drinker_status[df_data_corr.Alcohol_drinker_status == 1].count()
    non_alcohol_drinker = df_data_corr.Alcohol_drinker_status[df_data_corr.Alcohol_drinker_status == 0].count()
    mean_age = str(np.round(df_data_corr['Age'].mean(),2))+" ("+str(np.round(df_data_corr['Age'].std(),2))+")"
    mean_bmi = str(np.round(df_data_corr['BMI'].mean(),2))+" ("+str(np.round(df_data_corr['BMI'].std(),2))+")"
    mean_towsend = str(np.round(df_data_corr['Townsend_index'].mean(),2))+" ("+str(np.round(df_data_corr['Townsend_index'].std(),2))+")"

    df_data_trans, y_meas, y_meas_lat, x_meas, cov, cp = load_data_type('SZ', 'diagnosed', 'wide', meas_gr)
    cova = ['Age', 'Age_squared', 'Sex', 'Smoking_status', 'BMI', 
           'Eye_disorders', 'Genotype_array', 'Townsend_index', 
           'Genetic PC1', 'Genetic PC2', 'Genetic PC3', 'Genetic PC4', 'Genetic PC5', 
           'Genetic PC6', 'Genetic PC7', 'Genetic PC8', 'Genetic PC9', 'Genetic PC10', 
           'Macula_centered_left', 'Macula_centered_right', 'OCT_quality_left', 'OCT_quality_right']
    print('after diag scz:',len(df_data_trans))
    col_incl = y_meas_lat+x_meas+cova+['Participant ID']
    n_diag = len(df_data_trans)
    df_data_corr_scz = drop_nan(df_data_trans,col_incl)
    n_diag_nan = len(df_data_trans)-len(df_data_corr_scz)
    print('after diag scz nan:',n_diag_nan)
    df_data_corr_scz.to_csv('../output/data/df_mlm_scz.csv')

    _, n_filt, n_scz, sex_vcounts, t_sex, p_sex = population(df_data_trans,df_data_corr,df_data_corr_scz,'Appx_Table1')
    # read to writ to results df
    df_results = pd.read_csv('../output/data/retinflam_results.csv')
    df_results['n_nondiag'] = [n_nondiag]
    df_results['n_nondiag_nan'] = [n_nondiag_nan]
    df_results['n_diag'] = [n_diag]
    df_results['n_diag_nan'] = [n_diag_nan]
    df_results['mean_age'] = [mean_age]
    df_results['mean_bmi'] = [mean_bmi]
    df_results['mean_townsend_index'] = [mean_towsend]
    df_results['male'] = [sex_vcounts[1]]
    df_results['female'] = [sex_vcounts[2]]
    df_results['mcar_pval'] = [mcar_pval]
    df_results['tstat_sex'] = [np.round(t_sex,2)]
    df_results['tstat_sex_p'] = [p_sex]
    df_results['curr_smoker'] = [curr_smoker]
    df_results['prev_smoker'] = [prev_smoker]
    df_results['non_smoker'] = [non_smoker]
    df_results['curr_alcohol_drinker'] = [curr_alcohol_drinker]
    df_results['prev_alcohol_drinker'] = [prev_alcohol_drinker]
    df_results['non_alcohol_drinker'] = [non_alcohol_drinker]
    df_results.to_csv('../output/data/retinflam_results_m.csv')

    meas_gr = 'subfields'
    df_data_trans, y_meas, y_meas_lat, x_meas, cov, cp = load_data_type('SZ', 'none', 'wide', meas_gr)
    # Add covariates to df
    cova = ['Age', 'Age_squared', 'Sex', 'Smoking_status', 'Alcohol_drinker_status', 'BMI', 'Hypertension','Diabetes_mellitus',
           'Genotype_array', 'Townsend_index', 
           'Genetic PC1', 'Genetic PC2', 'Genetic PC3', 'Genetic PC4', 'Genetic PC5', 
           'Genetic PC6', 'Genetic PC7', 'Genetic PC8', 'Genetic PC9', 'Genetic PC10','OCT_quality_left', 'OCT_quality_right']
    col_incl = y_meas_lat+x_meas+cova+['Participant ID']
    df_data_corr, missing_vals = count_dropped_rows(df_data_trans, col_incl)
    df_data_corr.to_csv('../output/data/df_mlm.csv')

    # Collect data for inner (RNFL-INL) and outer (INL-RPE) retinal thickness
    meas_gr = 'inner_outer'
    df_data_trans, y_meas, y_meas_lat, x_meas, cov, cp = load_data_type('SZ', 'none', 'wide', meas_gr)
    df_data_trans['RNFL_mean'] = df_data_trans[['RNFL_left','RNFL_right']].mean(axis=1)
    df_data_trans['GC_IPL_mean'] = df_data_trans[['GC_IPL_left','GC_IPL_right']].mean(axis=1)
    df_data_trans['INL_mean'] = df_data_trans[['INL_left','INL_right']].mean(axis=1)
    df_data_trans['INL_RPE_mean'] = df_data_trans[['INL_RPE_left','INL_RPE_right']].mean(axis=1)
    df_data_trans['Image_quality_mean'] = df_data_trans[['OCT_quality_left','OCT_quality_right']].mean(axis=1)
    
    cova = ['Age', 'Age_squared', 'Sex', 'Smoking_status', 'Alcohol_drinker_status', 'BMI', 'Hypertension','Diabetes_mellitus','Fasting_time',
           'Genotype_array', 'Townsend_index', 
           'Genetic PC1', 'Genetic PC2', 'Genetic PC3', 'Genetic PC4', 'Genetic PC5', 
           'Genetic PC6', 'Genetic PC7', 'Genetic PC8', 'Genetic PC9', 'Genetic PC10', 
           'Macula_centered_left', 'Macula_centered_right', 'Image_quality_mean']
    
    ret_io_phen = ['RNFL_mean','GC_IPL_mean','INL_mean','INL_RPE_mean']
    pathprs = ['M43559', 'M36658', 'M6557', 'M24927', 'M18933', 'M15140', 'M25305', 'M24111', 'M17761']
    inflam_mark = ['SII','NLR','PLR','LMR','Monocytes','Neutrophils','CRP']
    col_incl = y_meas_lat+x_meas+cova+['Participant ID']+ret_io_phen+inflam_mark+pathprs
    print('columns included:',col_incl)

    #df_data_corr = drop_nan(df_data_trans,col_incl)
    df_data_corr, missing_vals = count_dropped_rows(df_data_trans, col_incl)

    # Collinarity between and distribution of retinal phenotypes
    collinarity_dist(df_data_corr,ret_io_phen,'Appx_FigS2')
    print('coll plot done!')

    # Partial regression average inner outer retina
    parreg_oct(df_data_corr,['INL_RPE_mean'],x_meas,cova,figname='Fig1c')

    df_data_corr.to_csv('../output/data/df_mlm_inner_outer.csv')
    

