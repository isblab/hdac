import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['HDAC1.0_C-term_377-482',
        'HDAC1.1_C-term_377-482',
        'MTA1_BAH_1-164',
        'MTA1.1_BAH_1-164',
        'HDAC1.0_1-376',
        'HDAC1.1_1-376',
        'MTA1_ELM2-SANT_165-333',
        'MTA1.1_ELM2-SANT_165-333',
        'MTA1_mid_334-467',
        'MTA1.1_mid_334-467',
        'MBD3_C-term',
        'MBD3_mid_uns',
        'MBD3_mid_struc',
        'MBD3_N-term',
        'GATAD2B_1-136',
        'GATAD2B_137-178',
        'GATAD2B_179-281',
        'MTA1_R1_468-546',
        'MTA1.1_R1_468-546',
        'MTA1_between_R1_and_R2',
        'MTA1.1_between_R1_and_R2',
        'MTA1_R2_670-691',
        'MTA1.1_R2_670-691',
        'MTA1_C-term_692-715',
        'MTA1.1_C-term_692-715',
        'RBBP4.0_1-425',
        'RBBP4.1_1-425',
        'RBBP4.2_1-425',
        'RBBP4.3_1-425']

# Set visualization thresholds
threshold = {'HDAC1.0_C-term_377-482':0.0115,           # ~29%
        'HDAC1.1_C-term_377-482':0.0115,                # ~29%
        'MTA1_BAH_1-164':0.05,                          # ~19%
        'MTA1.1_BAH_1-164':0.0292,                      # ~22%
        'HDAC1.0_1-376':0.0286,                         # ~10%
        'HDAC1.1_1-376':0.06,                           # ~16%
        'MTA1_ELM2-SANT_165-333':0.06,                  # ~23%
        'MTA1.1_ELM2-SANT_165-333':0.04,                # ~12%
        'MTA1_mid_334-467':0.0151,                      # ~29%
        'MTA1.1_mid_334-467':0.0102,                    # ~29%
        'MBD3_C-term':0.006,                            # ~21%
        'MBD3_mid_uns':0.025,                           # ~19%
        'MBD3_mid_struc':0.005,                         # ~17%
        'MBD3_N-term':0.012,                            # ~19%
        'GATAD2B_1-136':0.013,                          # ~20%
        'GATAD2B_137-178':0.00553,                      # ~17%
        'GATAD2B_179-281':0.015,                        # ~21%
        'MTA1_R1_468-546':0.011,                        # ~23%
        'MTA1.1_R1_468-546':0.009,                      # ~24%
        'MTA1_between_R1_and_R2':0.015,                 # ~24%
        'MTA1.1_between_R1_and_R2':0.017,               # ~20%
        'MTA1_R2_670-691':0.004,                        # ~22%
        'MTA1.1_R2_670-691':0.002,                      # ~15%
        'MTA1_C-term_692-715':0.0024,                   # ~26%
        'MTA1.1_C-term_692-715':0.0017,                 # ~18%
        'RBBP4.0_1-425':0.06,                           # ~21%
        'RBBP4.1_1-425':0.05,                           # ~21%
        'RBBP4.2_1-425':0.05,                           # ~22%
        'RBBP4.3_1-425':0.02, }                         # ~22%


# Color of each protein/domain
col = {'HDAC1.0_1-376':'yellow',\
    'HDAC1.1_1-376':'yellow',\
    'HDAC1.0_C-term_377-482':'#c2a300',\
    'HDAC1.1_C-term_377-482':'#c2a300',\
    'MTA1_BAH_1-164':'#FFA500',\
    'MTA1.1_BAH_1-164':'#FFA500',\
    'MTA1_ELM2-SANT_165-333':'#FF8C00',\
    'MTA1.1_ELM2-SANT_165-333':'#FF8C00',\
    'MTA1_mid_334-467':'#f27130',\
    'MTA1.1_mid_334-467':'#f27130',\
    'MTA1_R1_468-546':'#ff5f24',\
    'MTA1.1_R1_468-546':'#ff5f24',\
    'MTA1_between_R1_and_R2':'#FF6347',\
    'MTA1.1_between_R1_and_R2':'#FF6347',\
    'MTA1_R2_670-691':'#ff4f38',\
    'MTA1.1_R2_670-691':'#ff4f38',\
    'MTA1_C-term_692-715':'#ff4040',\
    'MTA1.1_C-term_692-715':'#ff4040',\
    'RBBP4.0_1-425':'#0984ad',\
    'RBBP4.1_1-425':'#0984ad',\
    'RBBP4.2_1-425':'#88d7f7',\
    'RBBP4.3_1-425':'#88d7f7',\
    'MBD3_C-term':'dark green',
    'MBD3_mid_uns':'#84ff6b',
    'MBD3_mid_struc':'#4ca84c',
    'MBD3_N-term':'#4affa4',
    'GATAD2B_1-136':'#cccccc',
    'GATAD2B_137-178':'#969696',
    'GATAD2B_179-281':'#636363' }


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

runCommand('open EMD_22895.mrc')
runCommand('volume #'+str(i)+' step 1 level 0.05 style surface color "#aeaeae" transparency 0.8')
