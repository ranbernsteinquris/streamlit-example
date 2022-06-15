import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import simulation.quris_config as cfg
import seaborn as sns
#sns.set_theme()

#pd.set_option('precision', 2)

import streamlit as st
#st.set_page_config(layout="wide")

import os
import pandas_gbq

def validate_well_allocation_is_mutualy_exclusive(plate_measurements, group_col_name):
    wells_groups = plate_measurements.groupby(group_col_name).well.apply(set).values
    for i1, g1 in enumerate(wells_groups):
        for i2, g2 in enumerate(wells_groups):
            if i1 == i2:
                continue
            if len(g1.intersection(g2)) > 0:
                return False
    return True

os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = 'braided-course-bigqueryai.json'
project_id = "braided-course-312015"
table_id = 'sensor.data'
df = pandas_gbq.read_gbq('select * from sensor.data', project_id=project_id)

############################################################################################################################################
st.title('Description')
'Total rows in DB', len(df)

ms8 = df[df['study'] == 'MS8']
'Total rows in ms8', len(ms8)

ms8 = ms8[ms8['treatment'] != 'empty']
'Total rows without empty', len(ms8)

'Rows per analyte:'
tmp = ms8.groupby('analyte').value.count()
tmp

lactate_days = sorted(ms8[ms8['analyte'] == 'lactate'].day.unique().tolist())
'Lactate Days', lactate_days

'Number of measurements per (plate, day, analyte):'
tmp = ms8.groupby(['day', 'plate', 'analyte']).value
tmp = tmp.count().sort_index().reset_index()
st.dataframe(tmp)
#tmp
ms8['group'] = ms8.apply(lambda row: row.compound.replace('merck_compound_', 'C') + '_' + str(row.concentration_uM), axis=1)
ms8['well_id'] = ms8.apply(lambda row: str(row.plate) + '_' + str(row.well), axis=1)

############################################################################################################################################
st.title('Validations')

def check_unique(key, analyte):
    tmp = ms8[ms8['analyte'] == analyte]
    tmp = tmp.set_index(key)
    tmp['has_dups'] = ms8.groupby(key).value.count() > 1
    tmp = tmp[tmp['has_dups'] == True].sort_index()
    #'len(tmp)', len(tmp)
    if len(tmp) == 0:
        'Passed: Single value check passed for', analyte
    else:
        'Failed: found multiple assignments for', analyte
        tmp


ctg_key = ['plate', 'row', 'col', 'analyte']
check_unique(ctg_key, 'atp')

lactate_key = ['plate', 'row', 'col', 'analyte', 'day']
check_unique(lactate_key, 'lactate')

plates = ms8.plate.unique()

if len(set(ms8[(ms8['day'] == 0) & (ms8['compound'] == 'medium')].well)) == 96*len(plates):
    'Passed: All wells were measured in medium only on day 0'
else:
    'Failed: Not all wells were measured in medium only on day 0'

failed = False
for plate in plates:
    plate_measurements = ms8[ms8['plate'] == plate]
    plate_measurements_without_day_0 = plate_measurements[plate_measurements['day'] != 0]
    if not validate_well_allocation_is_mutualy_exclusive(plate_measurements_without_day_0, 'group'):
        failed = True
        'Failed: Wells allocation is not mutually exclusive for plate', plate
if not failed:
    'Passed: Wells allocation is mutually exclusive for plate', plate

found_missing_lac = False
for _, ctg_row in ms8[ms8['analyte'] == 'atp'].iterrows():
    expected_lac_days = set([day for day in lactate_days if day <= ctg_row.day])
    actual_lac_days = set(ms8[(ms8['analyte'] == 'lactate') & (ms8['well_id'] == ctg_row.well_id)].day.values)
    if not expected_lac_days.issubset(actual_lac_days):
        found_missing_lac = True
        missing_lac_well_id = ctg_row.well_id
        missing_days = str(expected_lac_days.difference(actual_lac_days))
        break
    if found_missing_lac:
        break
if found_missing_lac:
    'Failed: Found missing lactate value for well_id and day', ctg_row.well_id, missing_days
else:
    'Passed: Each ctg measurement has all preceding lactate values'
############################################################################################################################################
st.title('Lactate')

def plot_analyte_over_time(analyte):
    fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(5,7))
    fig.supxlabel('Days')
    fig.supylabel('Lactate')
    #plot_lactate_over_time(1, axes[0])
    ms8_without_day_0 = ms8[(ms8['day'] != 0) & (ms8['analyte'] == analyte)]
    #groups = ms8_without_day_0['group'].unique()
    days = ms8_without_day_0['day'].unique()
    for plate in plates:
        axe = axes[plate-1]
        plate_df = ms8_without_day_0[ms8_without_day_0['plate'] == plate]
        for group in plate_df['group'].unique():
            means = []
            for day in days:
                means.append(plate_df[(plate_df['day'] == day) & (plate_df['group'] == group)].value.mean())
            axe.scatter(days, means, label=group)
        axe.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        axe.set_title('plate {}'.format(plate))
    st.pyplot(fig)
'In plate one the toxic treatments (150 and 300 uM) can be obtained by the trend change on day 3'

st.header('Production over days')
plot_analyte_over_time('lactate')

#ms8_without_day_0 = ms8[(ms8['day'] != 0) & (ms8['analyte'] == analyte)]



############################################################################################################################################
st.title('CTG')

st.header('ATP count over days')
fig, axe = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10, 5))
ctg_data = ms8[(ms8['day'] != 0) & (ms8['analyte'] == 'atp')]
days = [1, 4, 7]
for plate in [1,2]:
    centroids = []
    for day in days:
        centroid = [ctg_data[(ctg_data['plate'] == plate) & (ctg_data['day'] == day)].value.mean()]
        centroids.append(centroid)
    axe.scatter(days, centroids, label='P{}'.format(plate))
axe.legend()
fig

st.header('ATP count over concentrations')
fig, axe = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10, 5))
concentrations = ctg_data['concentration_uM'].unique()
for plate, color in zip([1,2], ['g','b']):
    for day, shape in zip([4, 7], ["o" , "v"]):
        centroids = []
        for concentration in concentrations:
            if concentration != np.nan:
                centroid = [ctg_data[(ctg_data['plate'] == plate) & (ctg_data['day'] == day) & (ctg_data['concentration_uM'] == concentration)].value.mean()]

            centroids.append(centroid)
        axe.scatter(concentrations, centroids, label='P{}_D{}'.format(plate, day), marker=shape, color=color)

for plate, color in zip([1,2], ['g','b']):
    for day in days:
        centroids = []
        centroid = [ctg_data[(ctg_data['plate'] == plate) & (ctg_data['day'] == day) & (ctg_data['compound'] == 'Medium')].value.mean()]
axe.legend()
axe.set_xlabel('Compound Concentration')
axe.set_ylabel('ATP')
axe.set_title('CTG comparison')
fig

'Insight: Is there a linear relationship between atp and concentration? We need more mild concentrations, and more repeations ' \
'per data point.'

def scatter_ctg(day):
    d7_atp = ms8[(ms8['day'] == day) & (ms8['analyte'] == 'atp')]
    groups = sorted(d7_atp.group.unique(), key=lambda x: float(x.split('_')[-1]))
    fig, axe = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(13, 5))
    for group in groups:
        group_values = d7_atp[d7_atp.group == group].value
        axe.scatter([group]*len(group_values), group_values)
    fig

st.header('ATP at the middle (D4)')
scatter_ctg(4)

st.header('ATP at the end point (D7)')
scatter_ctg(7)

'Insight: Medium and DMSO are currently useless baselines. The idea of allocating the in the first and last row, while ' \
'the treatments in columns should be reconsidered.'


############################################################################################################################################
st.title('CTG and Lactate correlation')

import scipy.stats as stats
def correlate_lac_and_ctg(lac_day, ctg_day, axe):
    ctg_rows = ms8[(ms8.day == ctg_day) & (ms8.analyte == 'atp')]
    wells = ctg_rows.well_id.values
    ctg_means = ctg_rows.groupby(['group']).value.mean()

    lac = []
    lac_means = ms8[(ms8.day == lac_day) & (ms8.analyte == 'lactate') & (ms8.well_id.isin(wells))].groupby('group').value.median()
    lac = lac_means.values
    ctg = ctg_means.values
    r, p =stats.pearsonr(lac, ctg)
    groups = ctg_means.index.values
    for l, c, g in zip(lac, ctg, groups):
        axe.scatter([l], [c], label=g)
    #axe.set_xlabel('lactate')
    #axe.set_ylabel('atp')
    axe.set_title('lactate day {}, ctg day {}, r:{:.2f}, p:{:.3f}'.format(lac_day, ctg_day, r, p))
    axe.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig, axes = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(5, 15))
correlate_lac_and_ctg(3,4, axes[0])
correlate_lac_and_ctg(5,7, axes[1])
correlate_lac_and_ctg(7,7, axes[2])
fig
#for plate in [1,2]:
#    lac = lac_d5[lac_d5.plate == plate].set_index(['well']).sort_index().stack()
#    ctg = ctg_d4[ctg_d4.plate == plate].set_index(['well']).sort_index().stack()

"""
############################################################################################################################################
st.title('Predict CTG from Lactate')

tmp = ms8[(ms8['analyte'] == 'atp') & (ms8['day'] == 7)][['well_id', 'value']]
wells = tmp['well_id'].values
y = tmp['value'].values
X = []
days_for_train = [0,1,3,5,7]

for well_id in wells:
    feature_vec = []
    for day in days_for_train:
        curr_val = ms8[(ms8['well_id'] == well_id) & (ms8['analyte'] == 'lactate') & (ms8['day'] == day)].value.values[0]
        if len(feature_vec) > 0:
            feature_vec.append(curr_val-feature_vec[-1])
        feature_vec.append(curr_val)
    X.append(feature_vec)

#from sklearn.linear_model import LinearRegression
#reg = LinearRegression().fit(X, y)

y_binary = [atp > 400 for atp in y]
from sklearn.ensemble import GradientBoostingClassifier
clf = GradientBoostingClassifier()
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y_binary, test_size=12/96, random_state=42)
clf.fit(X_train, y_train)
print(clf.score(X_test, y_test))
"""
