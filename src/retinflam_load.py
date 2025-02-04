"""
@author: Finn Rabe
@email: finn.rabe@bli.uzh.ch
@terms: CC-BY-NC-ND
"""
import os
from pathlib import Path
from statsmodels.datasets import longley
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.decomposition import PCA
from scipy.stats import zscore
import pingouin as pg
from scipy.stats import shapiro, kstest
from statsmodels.compat import lzip
from statsmodels.stats.api import het_breuschpagan
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.graphics.api import abline_plot
import plotly.figure_factory as ff
from matplotlib import gridspec
from stargazer.stargazer import Stargazer
import warnings
import glob
from scipy.stats import chi2_contingency
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LinearRegression
from pyampute.exploration.mcar_statistical_tests import MCARTest
warnings.simplefilter(action='ignore', category=FutureWarning)

def set_fig_style():
    mpl.rcParams.update(mpl.rcParamsDefault)
    #seaborn.reset_orig
    #sns.set_context('poster')
    #sns.set(font_scale=3)
    #sns.set_theme(font='sans-serif') #set font type for all plots
    #sns.set_style("whitegrid", {'axes.grid' : True})

    # set resolution and figure sizes
    #%config InlineBackend.figure_format = 'retina'
    rcParams['figure.dpi'] = 300
    px = 1/rcParams['figure.dpi']  # pixel in inches
    fig_size = (1200*px,750*px)
    tick_size = 5
    label_size = 6
    title_size = 8
    val_size = 5
    line_size = 3

    # set color palettes
    dist_plots = 'grey'
    scatter_plots = 'grey'
    reg_line = 'darkgrey'
    thr_line = 'orange'
    return fig_size, tick_size, label_size, title_size, val_size, line_size, dist_plots, scatter_plots, reg_line, thr_line

def excl_antipsychotics(df_data_filt_wide):
    """exclude participants with prescribed antipsychotics"""
    antpych_lst = ['amisulpride', 'aripiprazole', 'asenapine', 'benperidol', 'brexpiprazole',\
        'cariprazine', 'chlorpromazine', 'clopenthixol', 'clozapine', 'flupenthixol',\
        'fluphenazine', 'haloperidol', 'iloperidone', 'levomepromazine', 'loxapine',\
        'lurasidone', 'molindone', 'olanzapine', 'paliperidone', 'quetiapine',\
        'penfluridol', 'perazine', 'perphenazine', 'pimozide', 'prochlorperazine', 'risperidone',\
        'sertindole', 'sulpiride', 'thioridazine', 'thiothixene', 'trifluoperazine',\
        'ziprasidone', 'zotepine', 'zuclopenthixiol']
    df_excl_ant = df_data_filt_wide[~df_data_filt_wide['Treatment/medication code | Instance 0'].str.contains('|'.join(antpych_lst),na=False)]
    return df_excl_ant

def general_formatting(prs_type,disorder,format_type,meas_gr):
     ## Create dataframe to append results values to
    df_results = pd.DataFrame()

    # plot formatting
    fig_size, tick_size, label_size, title_size, val_size, line_size, dist_plots, scatter_plots, reg_line, thr_line = set_fig_style()
    # Load data
    #df_data_filt = pd.read_csv('../data/retinflam_ukb_data.csv')
    df_data_filt = pd.read_csv('../data/retinflam_ukb_data_revised.csv')

    # Proposal: Filter by macular volume, because if it occurs, it usually also entails all other OCT measures
    df_data_filt = df_data_filt.dropna(subset=['Overall macular thickness (left) | Instance 0', 'Overall macular thickness (right) | Instance 0'])
    n_total = len(df_data_filt)
    print('available macula oct data:', n_total)
    df_results['n_total'] = [n_total]

    # add all general PRS
    #tip: use "for FILENAME in *; do mv $FILENAME $FILENAME.txt; done" for prs to txt conversion
    df_pathprs = pd.DataFrame()
    prs_lst = glob.glob('../data/generalPRS/*.txt')
    for f in range(len(prs_lst)):
        file_name = [os.path.basename(x) for x in glob.glob(prs_lst[f])]
        fsplit = file_name[0].split("_")
        #pname = fsplit[0]+'_PRS'
        pname = 'PRSSZ_2022'
        data_pathprs2022 = pd.read_csv(prs_lst[f], delim_whitespace=True)
        data_pathprs2022.columns = ["FieldID", "Participant ID", "Allele_CT", "Allele_Dosage", pname]
        df_pathprs[pname] = data_pathprs2022[pname]
    df_pathprs["Participant ID"] = data_pathprs2022["Participant ID"]
    df_data_filt = pd.merge(df_data_filt, df_pathprs, on='Participant ID', how='left')

    ## best-fit general PRS
    df_pathprs_bfit = pd.DataFrame()
    prs_lst = glob.glob('../data/generalPRS/best_fit/*best.txt')
    for f in range(len(prs_lst)):
        file_name = [os.path.basename(x) for x in glob.glob(prs_lst[f])]
        fsplit = file_name[0].split("_")
        pname = 'PRSSZ_bestfit'
        data_pathprs2022 = pd.read_csv(prs_lst[f], delim_whitespace=True)
        data_pathprs2022.columns = ["FieldID", "Participant ID", "In_Regression", pname]
        df_pathprs_bfit[pname] = zscore_data(data_pathprs2022[pname])
    df_pathprs_bfit["Participant ID"] = data_pathprs2022["Participant ID"]
    df_data_filt = pd.merge(df_data_filt, df_pathprs_bfit, on='Participant ID', how='left')

    # clumped based thresholded PRS
    df_prs_th = pd.DataFrame()
    #th_order = ['0_1','0_001','0_2','0_3','0_4','0_5','0_05','1']
    th_order = ['0_3','0_05','0_1','0_4','0_2','0_001','0_5','1']
    clump_set = 'clumped_1'
    prs_lst = glob.glob('../data/clumped_based_nscores/'+clump_set+'*.txt')
    for f in range(len(prs_lst)):
        file_name = [os.path.basename(x) for x in glob.glob(prs_lst[f])]
        #fsplit = file_name[0].split("_")
        pname = 'PRSSZ_'+th_order[f]
        #print('comp: filname | label:', file_name,pname)
        data_pathprs2022 = pd.read_csv(prs_lst[f], delim_whitespace=True)
        data_pathprs2022.columns = ["FieldID", "Participant ID", "Allele_CT", "Allele_Dosage", pname]
        df_prs_th[pname] = zscore_data(data_pathprs2022[pname])
    df_prs_th["Participant ID"] = data_pathprs2022["Participant ID"]
    df_data_filt = pd.merge(df_data_filt, df_prs_th, on='Participant ID', how='left')


    ## pathway-specific PRS results from PRSet
    df_pathprs_comp  = pd.DataFrame(
                             {'Pathway': [],
                            'Num_SNP': [],
                            'Phenotype': [],
                            'Estimate': [],
                            'SE': [],
                            'selfcontained.p': [],
                            'competitive.p': []
                             })

    pathprs_comp_lst = glob.glob('../data/pathwayPRS/competitivep/z_scored_revised/*.txt')
    for f in range(len(pathprs_comp_lst)):
        file_name = [os.path.basename(x) for x in glob.glob(pathprs_comp_lst[f])]

        fsplit = file_name[0].split("_")
        pname = fsplit[0]+'_PRS'
        pc = fsplit[1] #fsplit[-3].split(".")[0]
        #basename = fsplit[0]+'_Bestfit_PRS'
        #if pname in pf:
        #    data_pathprs2022 = pd.read_csv(pathprs_comp_lst[f], sep="\t")
        #else:
        data_pathprs2022 = pd.read_csv(pathprs_comp_lst[f], sep="\t")
        num_snp = data_pathprs2022.iloc[1]['Num_SNP']
        pathB = data_pathprs2022.iloc[1]['Coefficient']
        pathsd = data_pathprs2022.iloc[1]['Standard.Error']
        selfcp = data_pathprs2022.iloc[1]['P']
        comp = data_pathprs2022.iloc[1]['Competitive.P']

        new_row = pd.DataFrame({
            'Pathway': [pname],
            'Num_SNP': [num_snp],
            'Phenotype': [pc],
            'Estimate': [pathB],
            'SE': [pathsd],
            'selfcontained.p': [selfcp],
            'competitive.p': [comp]
        })

        df_pathprs_comp = pd.concat([df_pathprs_comp, new_row], ignore_index=True)
    
    df_pathprs_comp.to_csv('../output/data/pathway_comp.csv')
    
    # add pathway-specific PRS
    df_pathprs = pd.DataFrame()
    pathprs_lst = glob.glob('../data/pathwayPRS/new/*.txt')
    pathprs_lbl = {'CHROINFLAM': 'M15140',
        'MITO': 'M16257',
        'DOPPOSREG': 'M24111',
        'ABNOVAS':'M43559',
        'WNT': 'M25305',
        'NEUROINFLAM': 'M24927',
        'TGFB': 'M18933',
        'CATENIN': 'M17761',
        'ACUTEINFLAM': 'M6557',
        'CORART': 'M36658',
        'DIABET2': 'M36789'}
    for f in range(len(pathprs_lst)):
        file_name = [os.path.basename(x) for x in glob.glob(pathprs_lst[f])]
        fsplit = file_name[0].split("_")
        #pname = fsplit[0]+'_PRS'
        pname = pathprs_lbl[fsplit[0]]
        basename = fsplit[0]+'_Bestfit_PRS'
        data_pathprs2022 = pd.read_csv(pathprs_lst[f], sep=" ")
        data_pathprs2022.columns = ["FieldID", "Participant ID", basename, pname]
        df_pathprs[pname] = zscore_data(data_pathprs2022[pname])
        df_pathprs[basename] = zscore_data(data_pathprs2022[basename])
    df_pathprs["Participant ID"] = data_pathprs2022["Participant ID"]
    df_data_filt = pd.merge(df_data_filt, df_pathprs, on='Participant ID', how='left')

    # exclude individuals ID duplicates
    df_data_filt_wide = df_data_filt.drop_duplicates(subset=['Participant ID'])
    n_duplicat = len(df_data_filt) - len(df_data_filt_wide)
    print('number of duplicates:', n_duplicat)

    # exclude individuals based on genetic QC
    gen_qc = '../data/genetic_qc/sampleQC.id.txt'
    data_gen_qc= pd.read_csv(gen_qc, sep="\t")
    data_gen_qc.columns = ["FieldID", "Participant ID"]
    df_data_filt = df_data_filt_wide[df_data_filt_wide['Participant ID'].isin(data_gen_qc['Participant ID'])]
    n_qen_qc = len(df_data_filt_wide)-len(df_data_filt)
    df_results['n_qen_qc'] = [n_qen_qc]
    df_data_filt_wide = df_data_filt
    print('df after SNP QC:',len(df_data_filt_wide))

    print('df before array:',len(df_data_filt))
    # UKB utilizes two genotyping arrays, create binary column
    df_data_filt_wide['Genotype measurement batch'][df_data_filt_wide['Genotype measurement batch'].str.startswith('Batch_',na=False)] = 1
    df_data_filt_wide['Genotype measurement batch'][df_data_filt_wide['Genotype measurement batch'].str.startswith('UKBi',na=False)] = 2
    df_data_filt_wide['Genotype_array'] = df_data_filt_wide['Genotype measurement batch']#.astype(int)
    print('df after array:',len(df_data_filt_wide))

    #Add hypertension as covariate
    df_data_filt_wide['Hypertension'] = 0 
    df_data_filt_wide.loc[df_data_filt_wide['Diagnoses - ICD10'].str.contains('hypertension',na=False), 'Hypertension'] = 1
    n_hypertension = df_data_filt_wide['Hypertension'][df_data_filt_wide['Hypertension'] == 1].count()
    df_results['n_hypertension'] = [n_hypertension]

    #Add diabetes mellitus as covariate, consider only E11.3= with othamologic complications
    df_data_filt_wide['Diabetes_mellitus'] = 0
    df_data_filt_wide.loc[df_data_filt_wide['Diagnoses - ICD10'].str.contains('E11',na=False), 'Diabetes_mellitus'] = 1
    n_diabetestwo = df_data_filt_wide['Diabetes_mellitus'][df_data_filt_wide['Diabetes_mellitus'] == 1].count()
    df_results['n_diabetestwo'] = [n_diabetestwo]

    # Add disease that affected retina as covariate
    # based on https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6148
    #incl_crit = [float('nan'),'None of the above','Do not know','Prefer not to answer']
    print('df before eyedis:',len(df_data_filt_wide))
    eyedis_lst = ['Diabetes related eye disease','Glaucoma','Macular degeneration','Injury or trauma resulting in loss of vision']
    df_data_filt_wide['Eye problems/disorders | Instance 0'][~df_data_filt_wide['Eye problems/disorders | Instance 0'].isin(eyedis_lst)] = 0
    df_data_filt_wide['Eye problems/disorders | Instance 0'][df_data_filt_wide['Eye problems/disorders | Instance 0'].isin(eyedis_lst)] = 1
    n_eyedis = df_data_filt_wide['Eye problems/disorders | Instance 0'][df_data_filt_wide['Eye problems/disorders | Instance 0'] == 1].count()
    df_data_filt_wide = df_data_filt_wide[df_data_filt_wide['Eye problems/disorders | Instance 0']==0]
    df_results['n_eyedis'] = [n_eyedis]
    print('df after eyedis:',len(df_data_filt_wide))

    # Exclude individuals with extreme myopia and hypermetropia
    df_data_filt_wide.rename(columns={'Median spherical equivalent value (left) | Instance 2': 'MSE_left'
                                        }, inplace=True)
    df_data_filt_wide.rename(columns={'Median spherical equivalent value (right) | Instance 2': 'MSE_right'
                                        }, inplace=True)
    df_data_filt_wide['MSE_left'] = df_data_filt_wide['MSE_left'].fillna(0)
    df_data_filt_wide['MSE_right'] = df_data_filt_wide['MSE_right'].fillna(0)
    # Highly hyperopic SE >= 3
    n_hyperopic = len(df_data_filt_wide[(df_data_filt_wide['MSE_left'] >= 3) | (df_data_filt_wide['MSE_right'] >= 3)])
    print('n_hyperopic:', n_hyperopic)
    df_results['n_hyperopic'] = [n_hyperopic]
    df_data_filt_wide = df_data_filt_wide[(df_data_filt_wide['MSE_left'] <= 3) & (df_data_filt_wide['MSE_right'] <= 3)]
    print('df after hyperopic excl:', len(df_data_filt_wide))
    
    # Highly myopic SE <= -6, or "highly myopic"
    n_myopic = len(df_data_filt_wide[(df_data_filt_wide['MSE_left'] <= -6) | (df_data_filt_wide['MSE_right'] <= -6)])
    myopia_lst = ['highly myopic']
    n_myopia = df_data_filt_wide['Myopia diagnosis'][df_data_filt_wide['Myopia diagnosis'].isin(myopia_lst)].count()
    n_myopic = n_myopia+n_myopic
    df_results['n_myopic'] = [n_myopia+n_myopic]
    print('n_myopic:', n_myopic)
    df_data_filt_wide = df_data_filt_wide[(df_data_filt_wide['MSE_left'] >= -6) & (df_data_filt_wide['MSE_right'] >= -6) & (~df_data_filt_wide['Myopia diagnosis'].isin(myopia_lst))]
    print('df after excluding myopic:',len(df_data_filt_wide))
    df_results['n_my_hyperopic'] = [n_myopic+n_hyperopic]

    # Sex variables into numeric
    replace_sex= {'Male' : 1, 'Female' : 2}
    df_data_filt_wide['Sex'] = df_data_filt_wide['Genetic sex'].replace(replace_sex)
    #replace_bool = {False : 0, True : 1} 
    print('df after sex:',len(df_data_filt_wide))
    
    # Set Protein labels
    df_data_filt_wide.rename(columns={'C-reactive protein | Instance 0': 'CRP'
                                        }, inplace=True)
    
    df_data_filt_wide['SII'] = (df_data_filt_wide['Neutrophill count | Instance 0']*df_data_filt_wide['Platelet count | Instance 0'])/df_data_filt_wide['Lymphocyte count | Instance 0']
    df_data_filt_wide['NLR'] = df_data_filt_wide['Neutrophill count | Instance 0']/df_data_filt_wide['Lymphocyte count | Instance 0']
    df_data_filt_wide['PLR'] = df_data_filt_wide['Platelet count | Instance 0']/df_data_filt_wide['Lymphocyte count | Instance 0']
    df_data_filt_wide['LMR'] = df_data_filt_wide['Lymphocyte count | Instance 0']/df_data_filt_wide['Monocyte count | Instance 0']
    df_data_filt_wide['Monocytes'] = df_data_filt_wide['Monocyte count | Instance 0']
    df_data_filt_wide['Neutrophils'] = df_data_filt_wide['Neutrophill count | Instance 0']

    # Compute age
    df_data_filt_wide['Age'] = df_data_filt_wide['Age at recruitment']
    df_data_filt_wide['Age_squared'] = np.square(df_data_filt_wide['Age'])

    print('df before CRP:',len(df_data_filt_wide))

    #Log_transform plots of CRP data
    #df_data_filt_wide = df_data_filt_wide.rename({'C-reactive protein | Instance 0': 'CRP'},axis=1)
    pr_log = plot_log_transform(df_data_filt_wide,'CRP',fig_size,figname='Appx_FigureS1')
    df_data_filt_wide['CRP'] = pr_log

    # exclude non white british/irish participants
    print('df before ethn:',len(df_data_filt_wide))

    ethn_qc = '../data/genetic_qc/whitebrit.txt'
    data_ethn= pd.read_csv(ethn_qc, sep="\t")
    data_ethn.columns = ["Participant ID","Ethnicity"]
    data_ethn = data_ethn[(data_ethn['Ethnicity'] == 'British') | (data_ethn['Ethnicity'] == 'Irish')]
    df_data_filt = df_data_filt_wide[df_data_filt_wide['Participant ID'].isin(data_ethn['Participant ID'])]
    n_ethn = len(df_data_filt_wide)-len(df_data_filt)

    df_results['n_ethn'] = [n_ethn]
    #df_data_filt_wide = df_data_filt_wide[df_data_filt_wide['Ethnicity_combined'] == 'White British/Irish']
    df_data_filt_wide = df_data_filt
    print('df after ethn:',len(df_data_filt_wide))

    df_data_filt_wide = df_data_filt_wide.rename({'QC - Image quality (left) | Instance 0': 'OCT_quality_left'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'QC - Image quality (right) | Instance 0': 'OCT_quality_right'},axis=1)

    #exclude images with poorest 20% image quality
    pc_20_l = df_data_filt_wide['OCT_quality_left'].quantile(0.20)
    pc_20_r = df_data_filt_wide['OCT_quality_right'].quantile(0.20)
    n_octqc_left = df_data_filt_wide['OCT_quality_left'][df_data_filt_wide['OCT_quality_left'] <= pc_20_l].count()
    n_octqc_right = df_data_filt_wide['OCT_quality_right'][df_data_filt_wide['OCT_quality_right'] <= pc_20_r].count()
    #n_octqc = n_octqc_left + n_octqc_right
    n_octqc = len(df_data_filt_wide[(df_data_filt_wide['OCT_quality_left'] <= pc_20_l) | (df_data_filt_wide['OCT_quality_right'] <= pc_20_r)])
    df_results['n_octqc'] = n_octqc
    df_data_filt_wide = df_data_filt_wide[(df_data_filt_wide['OCT_quality_left'] > pc_20_l) & (df_data_filt_wide['OCT_quality_right'] > pc_20_r)]
    print('OCT qc excluded:', n_octqc)
    print('df after OCT qc :', len(df_data_filt_wide))

    df_data_filt_wide = df_data_filt_wide.rename({'QC - Macula center aline (left) | Instance 0': 'Macula_centered_left'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'QC - Macula center aline (right) | Instance 0': 'Macula_centered_right'},axis=1)

    # rename columns
    if prs_type == 'SZ':
        #df_data_filt_wide = df_data_filt_wide.rename({'Standard PRS for schizophrenia (SCZ)': 'PRSSZ'},axis=1) #use this for old PRS SZ from the year 2018
        df_data_filt_wide = df_data_filt_wide.rename({'PRSSZ_2022': 'PRSSZ'},axis=1)
        df_data_filt_wide['PRSSZ'] = zscore_data(df_data_filt_wide['PRSSZ'])
        df_data_filt_wide['Diagnoses - ICD10'][~df_data_filt_wide['Diagnoses - ICD10'].str.contains('F2',na=False)] = 0
        df_data_filt_wide['Diagnoses - ICD10'][df_data_filt_wide['Diagnoses - ICD10'].str.contains('F2',na=False)] = 1
        if disorder == 'diagnosed':
            df_data_filt_wide = df_data_filt_wide[df_data_filt_wide['Diagnoses - ICD10'] == 1]
        elif disorder == 'none':
            print('before diag filter:',len(df_data_filt_wide))
            df_data_filt_wide = df_data_filt_wide[df_data_filt_wide['Diagnoses - ICD10'] == 0]
        elif disorder == 'both':
            df_data_filt_wide = df_data_filt_wide

    # rename all other relevant columns
    df_data_filt_wide = df_data_filt_wide.rename({'Current eye infection | Instance 0': 'Eye_infection'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'Fasting time | Instance 0': 'Fasting_time'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'Smoking status | Instance 0': 'Smoking_status'},axis=1)
    df_data_filt_wide.Smoking_status[df_data_filt_wide.Smoking_status == 'Never'] = 0
    df_data_filt_wide.Smoking_status[df_data_filt_wide.Smoking_status == 'Previous'] = 1
    df_data_filt_wide.Smoking_status[df_data_filt_wide.Smoking_status == 'Current'] = 2
    df_data_filt_wide.Smoking_status[df_data_filt_wide.Smoking_status == 'Prefer not to answer'] = float('nan')

    df_data_filt_wide = df_data_filt_wide.rename({'Alcohol drinker status | Instance 0': 'Alcohol_drinker_status'},axis=1)
    df_data_filt_wide.Alcohol_drinker_status[df_data_filt_wide.Alcohol_drinker_status == 'Never'] = 0
    df_data_filt_wide.Alcohol_drinker_status[df_data_filt_wide.Alcohol_drinker_status == 'Previous'] = 1
    df_data_filt_wide.Alcohol_drinker_status[df_data_filt_wide.Alcohol_drinker_status == 'Current'] = 2
    df_data_filt_wide.Alcohol_drinker_status[df_data_filt_wide.Alcohol_drinker_status == 'Prefer not to answer'] = float('nan')

    # get first 10 genetic principle components 
    for gpc in range(10):
        df_data_filt_wide = df_data_filt_wide.rename({'Genetic principal components | Array '+str(gpc+1): 'Genetic PC'+str(gpc+1)},axis=1)
   
    df_data_filt_wide = df_data_filt_wide.rename({'Townsend deprivation index at recruitment': 'Townsend_index'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'Body mass index (BMI) | Instance 0': 'BMI'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'Eye problems/disorders | Instance 0': 'Eye_disorders'},axis=1)

    df_data_filt_wide = df_data_filt_wide.rename({'Standard PRS for cardiovascular disease (CVD)': 'CVDPRS'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({'Standard PRS for bipolar disorder (BD)': 'BDPRS'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({"Standard PRS for alzheimer's disease (AD)": 'ADPRS'},axis=1)
    df_data_filt_wide = df_data_filt_wide.rename({"Standard PRS for parkinson's disease (PD)": 'PDPRS'},axis=1)

    
    # exclude individuals using antipsychotics
    if disorder == 'none':
        df_data_filt_psy = excl_antipsychotics(df_data_filt_wide)
        n_antipsy = len(df_data_filt_wide) - len(df_data_filt_psy)
        df_data_filt_wide = df_data_filt_psy
        df_results['n_antipsy'] = [n_antipsy]
        print('n antipsychotics:', n_antipsy)
    
    #df_data_filt_wide = df_data_filt_wide.rename({"Treatment/medication code | Instance 0": 'Antipsychotics'},axis=1)

    ### Save dataframe with variables ###
    if disorder == 'none' and format_type == 'wide' and meas_gr == 'overall':
        df_results.to_csv('../output/data/retinflam_results.csv')

    return df_data_filt_wide

def oct_dict(meas_gr, retphen, side):
    if meas_gr == 'overall':
        ukb_lbl = [
            'Overall macular thickness '
        ]
    elif meas_gr == 'inner_outer':
        ukb_lbl = [
            'Average ganglion cell-inner plexiform layer thickness ',
            'Average retinal nerve fibre layer thickness ',
            'Average inner nuclear layer thickness ',
            'Average INL-RPE thickness ',
            'Overall average retinal pigment epithelium thickness '
        ]
    elif meas_gr == 'subfields':
        ukb_lbl = [
            'Macular thickness at the inner inferior subfield ',
            'Macular thickness at the outer inferior subfield ',
            'Macular thickness at the inner nasal subfield ',
            'Macular thickness at the outer nasal subfield ',
            'Macular thickness at the inner superior subfield ',
            'Macular thickness at the outer superior subfield ',
            'Macular thickness at the inner temporal subfield ',
            'Macular thickness at the outer temporal subfield ',
            'Macular thickness at the central subfield '
        ]
    side_br = '('+side+')'
    col_dict = dict()
    for lbl in range(len(retphen)):
        col_dict[retphen[lbl]+'_'+side] = ukb_lbl[lbl]+side_br+' | Instance 0'
    #col_dict['OCT_quality_'+side] =  'QC - Image quality '+side_br+' | Instance 0'
    return col_dict

def oct_formatting(df_data_filt_wide,meas_gr,retphen_unique,format_type, disorder):
    """Depending on the linear model, dataframe needs to be in long or wide fromat
    mixed linear model = long format (tile dataframe and add new column of each oct measure consisting of left and right and one binary column for eye side)
    robust linear model = wide format """
    
    retphen = retphen_unique #+ ['OCT_quality']
    side = ['left','right']
    # retphen_lat = [octm + '_' + eye for octm in retphen for eye in side]
    len_dataframe = len(df_data_filt_wide)
    
    if format_type == 'wide' or format_type == 'both':
        if disorder == 'none':
            ret_samp_dir = Path('../output/data/retphen_samp.csv')
            if meas_gr== 'overall' or meas_gr== 'inner_outer' and ret_samp_dir.is_file():
                df_samp_ret = pd.read_csv(ret_samp_dir)
            else:
                df_samp_ret = pd.DataFrame({'retinal_phenotype':[],'thickness':[]})
        y_meas_lat = [octm + '_' + eye for octm in retphen_unique for eye in side]
        df_data_filt = df_data_filt_wide.copy()
        for s in side:
            col_dict = oct_dict(meas_gr,retphen,s)
            for col in col_dict:
                if meas_gr == 'subfields':
                    if disorder == 'none':
                        df_samp_ret = df_samp_ret._append({
                                        'retinal_phenotype': col_dict[col].split('|')[0],
                                        'thickness': str(np.round(df_data_filt_wide[col_dict[col]].mean(),2))+" ("+str(np.round(df_data_filt_wide[col_dict[col]].std(),2))+")"
                                        }, ignore_index=True)
                    df_data_filt[col] = zscore_data(df_data_filt_wide[col_dict[col]])
                else:
                    if disorder == 'none':
                        df_samp_ret = df_samp_ret._append({
                                        'retinal_phenotype': col_dict[col].split('|')[0],
                                        'thickness': str(np.round(df_data_filt_wide[col_dict[col]].mean(),2))+" ("+str(np.round(df_data_filt_wide[col_dict[col]].std(),2))+")"
                                        }, ignore_index=True)
                    df_data_filt[col] = df_data_filt_wide[col_dict[col]]
        if disorder == 'none':
            print('length df_samp:', len(df_samp_ret))
            if len(df_samp_ret) < 5:
                df_samp_ret.to_csv('../output/data/retphen_samp.csv')
    elif format_type == 'long':
        y_meas_lat = retphen_unique #[octm + '_' + eye for octm in retphen[3:-1] for eye in side]
        # rearrange dataframe to long format to add each eye column into one
        df_data_filt = pd.DataFrame(np.tile(df_data_filt_wide, (2, 1)),columns=df_data_filt_wide.columns)
        # add column with binaries indicating which eye corresponds to which measure
        df_data_filt['Eye'] = [1]*len_dataframe + [0]*len_dataframe
        for col in range(len(retphen)):
            col_left = oct_dict(meas_gr,retphen,side[0])
            col_right = oct_dict(meas_gr,retphen,side[1])
            idx_left = retphen[col]+'_'+side[0]
            idx_right = retphen[col]+'_'+side[1]
            df_data_filt[retphen[col]] = zscore_data(df_data_filt_wide[col_left[idx_left]].append(df_data_filt_wide[col_right[idx_right]]).reset_index(drop=True))
    return df_data_filt, y_meas_lat

def sel_oct_group(meas_gr, format_type):
    """1. Examine overall measures 
       2. Examine macular subflieds"""#overall measures or macular subflieds
    
    if meas_gr == 'overall':
        y_meas = ['Macula']
        cp = sns.color_palette("crest")
    elif meas_gr == 'inner_outer':
        y_meas = ['GC_IPL','RNFL','INL','INL_RPE','RPE']
        cp = sns.color_palette("crest")
    elif meas_gr == 'subfields':
        y_meas = ['Inner_Inferior','Outer_Inferior','Inner_Nasal','Outer_Nasal','Inner_Superior','Outer_Superior','Inner_Temporal','Outer_Temporal','Central']
        cp = sns.color_palette("viridis")

    # Define covariates
    cov = ['Age',\
        'Age_squared',\
        'Sex',\
        'Fasting_time',\
        'Smoking_status',\
        'BMI',\
        'Eye_disorders',\
        'Genotype_array',\
        'Townsend_index',\
        'Genetic PC1',
        'Genetic PC2',
        'Genetic PC3',
        'Genetic PC4',
        'Genetic PC5',
        'Genetic PC6',
        'Genetic PC7',
        'Genetic PC8',
        'Genetic PC9',
        'Genetic PC10',
        'Macula_centered_left',
        'Macula_centered_right',
        'OCT_quality_left',
        'OCT_quality_right',
        'PRSSZ_0_001',
        'PRSSZ_0_05',
        'PRSSZ_0_1', 
        'PRSSZ_0_2',
        'PRSSZ_0_3',
        'PRSSZ_0_4', 
        'PRSSZ_0_5',
        'PRSSZ_1',
        'PRSSZ_bestfit',
        'CHROINFLAM_PRS',
        'MITO_PRS',
        'DOPPOSREG_PRS',
        'ABNOVAS_PRS',
        'WNT_PRS',
        'NEUROINFLAM_PRS',
        'TGFB_PRS',
        'CATENIN_PRS',
        'ACUTEINFLAM_PRS',
        'CORART_PRS',
        'DIABET2_PRS',
        'CHROINFLAM_Bestfit_PRS',
        'MITO_Bestfit_PRS',
        'DOPPOSREG_Bestfit_PRS',
        'ABNOVAS_Bestfit_PRS',
        'WNT_Bestfit_PRS',
        'NEUROINFLAM_Bestfit_PRS',
        'TGFB_Bestfit_PRS',
        'CATENIN_Bestfit_PRS',
        'ACUTEINFLAM_Bestfit_PRS',
        'CORART_Bestfit_PRS',
        'DIABET2_Bestfit_PRS',
        'SII',
        'NLR',
        'PLR',
        'LMR',
        'Monocyte_count',
        'Neutrophill_count',
        'CRP'
        ]
    if format_type == 'wide':
        cov = cov #+ ['OCT_quality_left'] + ['OCT_quality_right']
    else:
        cov = cov + ['Eye']
    
    return y_meas, cov, cp

def drop_nan(df_data_filt,columns):
    # filter dataframe by column names
    df_sel = df_data_filt[columns]

    # exclude all nan values
    df_nnan = df_sel.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    return df_nnan

def zscore_data(df_data_filt):
    # standarise/z-score data
    df_zcorr = zscore(df_data_filt,nan_policy='omit')
    #df_zcorr = df_data_filt.apply(zscore)
    return df_zcorr

def log_transform(df):
    data_log = np.log(df)
    return data_log

def plot_log_transform(df, col,fig_size,figname='FigureS1'):
    df_data = pd.DataFrame()
    log_data = log_transform(df[col])
    df_data[col+' levels'] = pd.concat([df[col],log_data], ignore_index=True)
    df_data['CRP data distribution'] = ['original']*len(df[col]) + ['log-transformed']*len(df[col])

    fig, g = plt.subplots(figsize=fig_size)
    g = sns.displot(df_data, x=col+' levels', hue='CRP data distribution', stat='density')
    sns.move_legend(obj = g, loc = 'lower left', bbox_to_anchor = (0.4, 0.7), frameon = True)
    plt.xlim((-5,30))

    #g.set_xticklabels(y_meas, fontsize=tick_size)
    #g.set_yticklabels(y_meas, fontsize=tick_size)

    #fig = g.get_figure()
    plt.savefig('../output/figures/'+figname+'_'+col+'_logtransform.png', format='png', bbox_inches='tight')
    return log_data

def regress_out(a, b):
    # use this function to regress out covariates from prs data
    """Regress b from a keeping a's original mean."""
    a = a.astype(float)
    b = b.astype(float)
    a_mean = a.mean()
    a = a - a_mean
    b = b - b.mean()
    b = np.c_[b]
    a_prime = a - b.dot(np.linalg.pinv(b).dot(a))
    return np.asarray(a_prime + a_mean).reshape(a.shape)

def pval_asterisks(pvalue, alpha=0.05):
    if pvalue < alpha/100:
        return "p < 0.0001****"
    elif pvalue < alpha/20:
        return "p < 0.001***"
    elif pvalue < alpha/10:
        return "p < 0.01**"
    elif pvalue < alpha:
        return "p < 0.05*"
    else:
        return "p = " +str(pvalue)
        
