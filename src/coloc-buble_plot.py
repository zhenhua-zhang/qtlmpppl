#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

input_file1 = '/home/umcg-zzhang/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/coloc/comparison_cohort-per_cell_type-with_sig_hits.txt'
dtfm1 = pd.read_csv(input_file1, sep=',', index_col=False)

input_file2 = '/home/umcg-zzhang/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/coloc/comparison_measurements-per_cell_type-with_sig_hits.txt'
dtfm2 = pd.read_csv(input_file2, sep=',', index_col=False)

input_file3 = '/home/umcg-zzhang/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/coloc/ccr5-vs-viral_load-per_cell_type.txt'
dtfm3 = pd.read_csv(input_file3)

dtfm = pd.concat([dtfm1, dtfm2, dtfm3])
dtfm = dtfm.loc[dtfm['Chrom'] == 3, :]

exclued_cell_type = ['RApR7p', 'RAnR7p', 'nTreg', 'Naive_CD8', 'M']
dtfm = dtfm.query('CellType not in @exclued_cell_type')

dtfm.loc[:, 'BothWithSigSNP'] = (dtfm.loc[:, 'WithSigHits'].str.find('0') == -1)

sb.set()
axe = sb.scatterplot(data=dtfm, x='Measurement', y='CellType',
                     size='PP4', sizes=(50, 200),
                     hue='BothWithSigSNP', edgecolor='none')

axe.legend(bbox_to_anchor=(1, 1))
axe.set_xlim(-1, 8)

xticklabels = ['PCP(BCG vs HIV)', 'MFI(BCG vs HIV)', 'BCG(PCP vs MFI)', 'HIV(PCP vs MFI)', 'PCP(BCGCCR5 vs viralLoad)', 'MFI(BCGCCR5 vs viralLoad)',
               'PCP(HIVCCR5 vs viralLoad)', 'MFI(HIVCCR5 vs viralLoad)']
axe.set_xticklabels(xticklabels)

axe.set_ylabel('Cell type')
axe.set_xlabel('Comparison')

plt.setp(axe.get_xticklabels(), rotation=45, ha='right')

fig = axe.get_figure()
fig.set_tight_layout(True)
fig.set_figwidth(7)
fig.set_figheight(6)

fig.savefig('/home/umcg-zzhang/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/coloc/coloc-per_cell_type-heatmap.svg')
fig.savefig('/home/umcg-zzhang/Documents/projects/wp_hiv_ccr5/outputs/age_gender/integration/coloc/coloc-per_cell_type-heatmap.png')
