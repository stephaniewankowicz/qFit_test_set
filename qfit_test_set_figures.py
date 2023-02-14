import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import glob
import os
import matplotlib.patches as mpatches
import argparse

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

def parse_args():
    parser = argparse.ArgumentParser(description='Compare the results of different water picking methods')
    parser.add_argument('--base_dir', type=str, help='The base directory of the water picking results')
    parser.add_argument('--comparison', type=str, help='The name of the comparison method')
    return parser.parse_args()

def rotamer_subset_compare(subset_rotamer, comparison):
    rotamer_ = []
    comparison_rotamer = comparison + '_rotamer'
    for i in subset_rotamer['original_PDB'].unique():
        tmp = subset_rotamer[subset_rotamer['original_PDB'] == i]
        for r in tmp['chain_resi'].unique():
            num_original = len(set(tmp[tmp['chain_resi'] == r][f'original_rotamer']))
            num_compairson = len(set(tmp[tmp['chain_resi'] == r][comparison_rotamer]))
            if (tmp[tmp['chain_resi']==r]['original_rotamer'] == tmp[tmp['chain_resi']==r][comparison_rotamer]).all() == True:
                rotamer = 'Same'
            elif bool(set(tmp[tmp['chain_resi']==r]['original_rotamer']) & set(tmp[tmp['chain_resi']==r][comparison_rotamer])) == True:
                if len(set(tmp[tmp['chain_resi']==r]['original_rotamer'])) > len(set(tmp[tmp['chain_resi']==r][comparison_rotamer])):
                    rotamer = 'Gain in Original'
                else:
                    rotamer = 'Gain in Comparison'
            else:
                rotamer = 'Different'
            rotamer_.append(tuple((i, r, rotamer, num_original, num_compairson)))
    rotamer_comp = pd.DataFrame(rotamer_, columns =['PDB', 'chain_resi', 'Rotamer', 'Num_Original', 'Num_Comp'])
    return rotamer_comp


def read_files(base_dir, comparison, file_start, file_end):
    #bfactors
    all_files = glob.glob(base_dir + '/' + comparison + "/*B_factors.csv")
    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, sep=',', header=0)
        df['PDB'] = filename[file_start:file_end] #file_start:file_end
        li.append(df)
    b_factor = pd.concat(li, axis=0, ignore_index=True)
    b_factor['category'] = comparison

    #rotamers
    all_files = glob.glob(base_dir + '/' + comparison + "/*_rotamer_output.txt")
    li = []
    for filename in all_files:
        try:
            df = pd.read_csv(filename, index_col=None, header=0, sep=':')
        except pd.errors.EmptyDataError:
            continue  
        df['PDB'] = filename[file_start:file_end]
        li.append(df)
    rotamer = pd.concat(li, axis=0, ignore_index=True)
    rotamer['category'] = comparison
    rotamer = rotamer[rotamer['residue']!= 'SUMMARY'].reset_index()
    split = rotamer['residue'].str.split(" ")
    for i in range(0, (len(rotamer.index)-1)):
        rotamer.loc[i,'chain'] = split[i][1]
        STUPID = str(rotamer.loc[i,'residue'])[3:6]
        tmp = []
        try:
            tmp = (int(''.join(i for i in STUPID if i.isdigit())))
        except:
            newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in STUPID)
            tmp = [float(i) for i in newstr.split()]
        rotamer.loc[i, 'resi'] = tmp

    #rvalues
    all_files = glob.glob(base_dir + '/' + comparison +  "/*_rvalues.csv")
    li = []

    for filename in all_files:
        df = pd.read_csv(filename, index_col = None, sep=',', header=0)
        df['PDB'] = filename[file_start:file_end]
        li.append(df)
    r_values = pd.concat(li, axis=0, ignore_index=True)
    return b_factor, r_values, rotamer


def compare_rvalues(comparison_rvalues, original_rvalues, comparison, base_dir):
    original_rvalues.columns = 'original_' + original_rvalues.columns.values
    comparison_rvalues.columns = comparison + '_' + comparison_rvalues.columns.values
    comparison_PDB = comparison + '_PDB'
    comparison_Rfree = comparison + '_Rfree_qFit'
    original_Rfree = 'original_Rfree_qFit'
    r_values = original_rvalues.merge(comparison_rvalues, left_on=['original_PDB'], right_on=[comparison_PDB])
    #calculating the difference in rvalues between each step
    r_values['Rfree_Differences'] = r_values[original_Rfree] - r_values[comparison_Rfree]

    #labeling structures that will be removed due to high rfree values
    r_values['Rfree_YN'] = np.where(((r_values[original_Rfree] - r_values[comparison_Rfree]) < -0.025), 1, 0)

    #plotting the difference in rvalues
    fig=plt.figure()
    sns.lmplot(original_Rfree, comparison_Rfree, data=r_values, hue='Rfree_YN', fit_reg=False)
    slope = 1
    y_intercept = 0
    x_fit = np.linspace(0.1, 0.3, 100)
    y_fit = slope * x_fit + y_intercept
    plt.plot(x_fit, y_fit)
    plt.xlabel('Original Rfree')
    plt.text(0.12, 0.25, 'Better Original qFit Rfree') 
    plt.text(0.24, 0.15, f'Better {comparison} Rfree')
    plt.ylabel('qFit Rfree')
    plt.savefig(base_dir + '/Refinement_scatterplot_comparison.png')

    fig = plt.figure()
    sns.kdeplot(r_values['Rfree_Differences'], shade=True)
    plt.xlabel('Rfree Difference')
    plt.savefig(base_dir + f'/Rfree_Differences_distribution_{comparison}.png')


def compare_rotamer(original_rotamer, comparison_rotamer, comparison, base_dir):    
    original_rotamer.columns = 'original_' + original_rotamer.columns.values
    comparison_rotamer.columns = comparison + '_' + comparison_rotamer.columns.values

    comparison_PDB = comparison + '_PDB'
    comparison_chain = comparison + '_chain'
    comparison_resi = comparison + '_resi'
    rotamer_all = original_rotamer.merge(comparison_rotamer, left_on=['original_PDB', 'original_chain', 'original_resi'], right_on=[comparison_PDB, comparison_chain, comparison_resi])
    rotamer_all['chain_resi'] = rotamer_all['original_chain'] + '_' + rotamer_all['original_resi'].astype(str)
    comp_rotamer= rotamer_subset_compare(rotamer_all, comparison)

    comp_rotamer.to_csv(f'{base_dir}/comp_rotamer.csv')

    comp_eval = comparison + '_evaluation'
    sns.countplot(x=comparison_PDB, data=comparison_rotamer[comparison_rotamer[comp_eval] == "OUTLIER"])
    plt.ylabel('PDB Id')
    plt.xlabel('Number of OUTLIER Rotamers')
    plt.savefig(base_dir + f'/Outlier_Rotamers_{comparison}.png')

    different = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Different'].index)
    same = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Same'].index)
    gain_comp = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Gain in Comparison'].index)
    gain_original = len(comp_rotamer[comp_rotamer['Rotamer'] == 'Gain in Original'].index)
    print(f'% Different:{different/len(comp_rotamer.index)}')
    print(f'% Same:{same/len(comp_rotamer.index)}')
    print(f'% Gain {comparison}:{gain_comp/len(comp_rotamer.index)}')
    print(f'% Gain Original: {gain_original/len(comp_rotamer.index)}')
    all_ = [same, different, gain_comp, gain_original]
    labels_ = ['Same', 'Different', f'Remodeled-Gain {comparison}', f'Remodeled-Loss {comparison}']

    fig = plt.figure()
    plt.bar([1,2,3,4], all_, color = ['r', 'b', 'g', 'c'])
    plt.xticks([1,2,3,4], labels_, rotation = 10) #, color=['black', 'red', 'blue', 'cyan']
    plt.savefig(f'{base_dir}/RotamerStatus_{comparison}.png')

def compare_bfactor(original_bfactor, comparison_bfactor, comparison, base_dir): 
    original_bfactor.columns = 'original_' + original_bfactor.columns.values
    comparison_bfactor.columns = comparison + '_' + comparison_bfactor.columns.values
    comparison_PDB = comparison + '_PDB'
    comparison_chain = comparison + '_chain'
    comparison_resi = comparison + '_resi'
    comparison_b = comparison + '_b_factor'
    comparison_altloc = comparison + '_num_altoc'

    b_factor_all = original_bfactor.merge(comparison_bfactor, left_on=['original_PDB', 'original_chain', 'original_resi'], right_on=[comparison_PDB, comparison_chain, comparison_resi])
    b_factor_all['bfactor_diff'] = b_factor_all['original_b_factor'] - b_factor_all[comparison_b]
    b_factor_all['altloc_diff'] = b_factor_all['original_num_altloc'] - b_factor_all[comparison_altloc]


    #compare altlocs
    fig = plt.figure()
    sns.histplot(b_factor_all['altloc_diff'])
    plt.text(-4.8, 6000, f'More Alt locs with {comparison}') 
    plt.text(2, 6000, 'More Altlocs with Original') 
    plt.savefig(f'{base_dir}/altloc_diff.png')

    #compare bfactor
    fig = plt.figure()
    sns.kdeplot(b_factor_all['bfactor_diff'])
    plt.axvline(x = 0, color = 'r', label = '')
    plt.text(-220, 0.03, f'Higher B-factor with {comparison}') 
    plt.text(100, 0.03, f'Higher B-factor Original')
    plt.savefig(f'{base_dir}/bfactor_diff.png')

    print(f'Average B-factor difference: {b_factor_all["bfactor_diff"].mean()}')
    print(f'Median Altloc difference: {b_factor_all["altloc_diff"].median()}')

    print('The following residues gain more than two altloc, you should look at them:')
    print(b_factor_all[b_factor_all['altloc_diff'] > 2][['original_PDB', 'original_chain', 'original_resi', 'original_num_altloc', comparison_altloc]])

    print('The following residues lose more than two altloc, you should look at them:')
    print(b_factor_all[b_factor_all['altloc_diff'] < -2][['original_PDB', 'original_chain', 'original_resi', 'original_num_altloc', comparison_altloc]])


    print('The following residues have a bfactor difference greater than 100, you should look at them:')
    print(b_factor_all[b_factor_all['bfactor_diff'] > 100][['original_PDB', 'original_chain', 'original_resi', 'original_b_factor', comparison_b]])
    print(b_factor_all[b_factor_all['bfactor_diff'] < -100][['original_PDB', 'original_chain', 'original_resi', 'original_b_factor', comparison_b]])


def main():
    args = parse_args()
    #read in original
    original_bfactor, original_rvalues, original_rotamer = read_files(args.base_dir, 'single_conf', 60, 64)  #51, 55
    #read in comparison
    comparison_bfactor, comparison_rvalues, comparison_rotamer = read_files(args.base_dir, args.comparison, 63, 67) #read in cctbx: 65-69 b-factor:60, 64
    compare_rvalues(comparison_rvalues, original_rvalues, args.comparison, args.base_dir)
    compare_rotamer(original_rotamer, comparison_rotamer, args.comparison, args.base_dir)
    compare_bfactor(original_bfactor, comparison_bfactor, args.comparison, args.base_dir)
if __name__ == '__main__':
    main()
