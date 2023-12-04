import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = [
        'RCOR1_Nterm',
        'RCOR1_ELM-SANT1',
        'RCOR1_240-310',
        'RCOR1_LINKER',
        'RCOR1_SANT2',
        'RCOR1_Cterm',
        'KDM1A_1-170',
        'KDM1A_SWIRM',
        'KDM1A_AOD1',
        'KDM1A_TOWER',
        'KDM1A_AOD2',
        'KDM1A_Cterm',
        'HDAC1_1-376',
        'HDAC1_C-term_377-482'        
        ]

# Set visualization thresholds
threshold = {
        'RCOR1_Nterm': 0.0295,            # 10%
        'RCOR1_ELM-SANT1': 0.0646,        # 10%
        'RCOR1_240-310': 0.0991,          # 10%
        'RCOR1_LINKER': 0.0517,           # 10%
        'RCOR1_SANT2': 0.0656,            # 10%
        'RCOR1_Cterm': 0.0653,            # 15%
        'KDM1A_1-170': 0.0795,            # 10%
        'KDM1A_SWIRM': 0.0755,            # 10%
        'KDM1A_AOD1': 0.0548,             # 10%
        'KDM1A_TOWER': 0.0628,            # 10%
        'KDM1A_AOD2': 0.0812,             # 10%
        'KDM1A_Cterm': 0.0234,            # 10%
        'HDAC1_1-376': 0.0835,            # 10%
        'HDAC1_C-term_377-482': 0.0383  } # 10%

# Color of each protein/domain
col = {
        'RCOR1_Nterm': '#fdbe85',
        'RCOR1_ELM-SANT1': '#fed98e',
        'RCOR1_240-310': '#fc8d59',
        'RCOR1_LINKER': '#fe9929',
        'RCOR1_SANT2': '#d95f0e',
        'RCOR1_Cterm': '#e6550d',
        'KDM1A_1-170': '#bae4b3',
        'KDM1A_SWIRM': '#c2e699',
        'KDM1A_AOD1': '#78c679',
        'KDM1A_TOWER': '#31a354',
        'KDM1A_AOD2': '#006d2c',
        'KDM1A_Cterm': '#2ca25f',
        'HDAC1_1-376': 'yellow',
        'HDAC1_C-term_377-482': '#c2a300' }


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

runCommand('open emd_10627.mrc')
runCommand('volume #'+str(i)+' step 1 level 0.008 style surface color "#aeaeae" transparency 0.8')