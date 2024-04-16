import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = [
        'SIN3A.0_PAH3',
        'SIN3A.0_HID',
        'SIN3A.1_PAH3',
        'SIN3A.1_HID',
        'SAP30.0_N-term',
        'SAP30.0_ZFD',
        'SAP30.0_SID',
        'SAP30.1_N-term',
        'SAP30.1_ZFD',
        'SAP30.1_SID',
        'SUDS3.0_N-term',
        'SUDS3.0_dimerization',
        'SUDS3.0_SID',
        'SUDS3.0_C-term',
        'SUDS3.1_N-term',
        'SUDS3.1_dimerization',
        'SUDS3.1_SID',
        'SUDS3.1_C-term',
        'HDAC1.0_1-376',
        'HDAC1.0_C-term_377-482',
        'HDAC1.1_1-376',
        'HDAC1.1_C-term_377-482',
        ]

# Set visualization thresholds
threshold = {
        'SIN3A.0_PAH3': 0.00692,              # 20%
        'SIN3A.0_HID': 0.0282,                # 20%
        'SIN3A.1_PAH3': 0.00706,              # 20%
        'SIN3A.1_HID': 0.02,                  # 20%
        'SAP30.0_N-term': 0.004172,           # 28%
        'SAP30.0_ZFD': 0.0117,                # 20%
        'SAP30.0_SID': 0.00606,               # 20%
        'SAP30.1_N-term': 0.004228,           # 28%
        'SAP30.1_ZFD': 0.0126,                # 20%
        'SAP30.1_SID': 0.00966,               # 20%      
        'SUDS3.0_N-term': 0.0022736,          # 28%
        'SUDS3.0_dimerization': 0.00388,      # 20%
        'SUDS3.0_SID': 0.00762,               # 20%
        'SUDS3.0_C-term': 0.004228,           # 28%
        'SUDS3.1_N-term': 0.0022204,          # 28%
        'SUDS3.1_dimerization': 0.00383,      # 10%
        'SUDS3.1_SID': 0.0095,                # 20%
        'SUDS3.1_C-term': 0.003276,           # 28%        
        'HDAC1.0_1-376': 0.0324,              # 20%
        'HDAC1.0_C-term_377-482': 0.004424,   # 28%
        'HDAC1.1_1-376': 0.0404,              # 20%
        'HDAC1.1_C-term_377-482': 0.0049  }   # 28%

# Color of each protein/domain
col = {
        'SIN3A.0_PAH3': '#fdbe85',
        'SIN3A.0_HID': '#fd8d3c',
        'SIN3A.1_PAH3': '#fdbe85',
        'SIN3A.1_HID': '#fd8d3c',
        'SAP30.0_N-term': '#fcae91',
        'SAP30.0_ZFD': '#fb6a4a',
        'SAP30.0_SID': '#de2d26',
        'SAP30.1_N-term': '#fcae91',
        'SAP30.1_ZFD': '#fb6a4a',
        'SAP30.1_SID': '#de2d26',
        'SUDS3.0_N-term': '#bae4b3',
        'SUDS3.0_dimerization': '#a6d96a',
        'SUDS3.0_SID': '#74c476',
        'SUDS3.0_C-term': '#31a354',
        'SUDS3.1_N-term': '#bae4b3',
        'SUDS3.1_dimerization': '#a6d96a',
        'SUDS3.1_SID': '#74c476',
        'SUDS3.1_C-term': '#31a354',
        'HDAC1.0_1-376': 'yellow',
        'HDAC1.0_C-term_377-482': '#c2a300',        
        'HDAC1.1_1-376': 'yellow',
        'HDAC1.1_C-term_377-482': '#c2a300'
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