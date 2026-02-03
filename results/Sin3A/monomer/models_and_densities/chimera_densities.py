import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = [
        'SIN3A_PAH3',
        'SIN3A_HID',
        'SAP30_N-term',
        'SAP30_ZFD',
        'SAP30_SID',
        'SUDS3_N-term',
        'SUDS3_dimerization',
        'SUDS3_SID',
        'SUDS3_C-term',
        'HDAC1_1-376',
        'HDAC1_C-term_377-482'        
        ]

# Set visualization thresholds
threshold = {
        'SIN3A_PAH3': 0.0134,               # 10%
        'SIN3A_HID': 0.0478,                # 20%
        'SAP30_N-term': 0.00694,            # 20%
        'SAP30_ZFD': 0.0167,                # 10%
        'SAP30_SID': 0.0136,                # 10%
        'SUDS3_N-term': 0.006608,           # 28%
        'SUDS3_dimerization': 0.00775,      # 10%
        'SUDS3_SID': 0.0112,                # 10%
        'SUDS3_C-term': 0.0111,             # 10%
        'HDAC1_1-376': 0.0583,              # 10%
        'HDAC1_C-term_377-482': 0.012712  } # 28%

# Color of each protein/domain
col = {
        'SIN3A_PAH3': '#fdbe85',
        'SIN3A_HID': '#fd8d3c',
        'SAP30_N-term': '#fcae91',
        'SAP30_ZFD': '#fb6a4a',
        'SAP30_SID': '#de2d26',
        'SUDS3_N-term': '#bae4b3',
        'SUDS3_dimerization': '#a6d96a',
        'SUDS3_SID': '#74c476',
        'SUDS3_C-term': '#31a354',
        'HDAC1_1-376': 'yellow',
        'HDAC1_C-term_377-482': '#c2a300'        
        }

runCommand('set bgcolor white')
i=0

for p in prots:
    print(p)
#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    runCommand('volume #'+str(i)+' level '+str(threshold[p]))
    # runCommand('volume #'+str(i)+' transparency '+str(transp[p]))
    runCommand('color '+col[p]+' #'+str(i))
    # runCommand('2dlabels create ' + str(p) + '_lab text ' + str(p) + ' color '+col[p]+' size 30 xpos .1 ypos ' + str(0.9 - i / 25.0))
    i += 1

# runCommand('open emd_10627.mrc')
# runCommand('volume #'+str(i)+' step 1 level 0.008 style surface color "#aeaeae" transparency 0.8')