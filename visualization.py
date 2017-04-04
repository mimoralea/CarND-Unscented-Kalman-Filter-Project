import plotly.offline as py
from plotly.graph_objs import *
import pandas as pd
import math
#py.init_notebook_mode()

my_cols=['px_est','py_est','vx_est','vy_est','px_meas','py_meas','px_gt','py_gt','vx_gt','vy_gt']
with open('./output-1.txt') as f:
    table_ekf_output1 = pd.read_table(f, sep='\t', header=None, names=my_cols, lineterminator='\n')

with open('./output-2.txt') as f:
    table_ekf_output2 = pd.read_table(f, sep='\t', header=None, names=my_cols, lineterminator='\n')

#estimations
trace11 = Scatter(
    x=table_ekf_output1['px_est'],
    y=table_ekf_output1['py_est'],
    xaxis='x2',
    yaxis='y2',
    name='EKF-Estimate'
)
#estimations
trace12 = Scatter(
    x=table_ekf_output2['px_est'],
    y=table_ekf_output2['py_est'],
    xaxis='x2',
    yaxis='y2',
    name='EKF-Estimate'
)

#Measurements
trace21 = Scatter(
    x=table_ekf_output1['px_meas'],
    y=table_ekf_output1['py_meas'],
    xaxis='x2',
    yaxis='y2',
    name = 'Measurements',
    mode = 'markers'
)
#Measurements
trace22 = Scatter(
    x=table_ekf_output2['px_meas'],
    y=table_ekf_output2['py_meas'],
    xaxis='x2',
    yaxis='y2',
    name = 'Measurements',
    mode = 'markers'
)

#Measurements
trace31 = Scatter(
    x=table_ekf_output1['px_gt'],
    y=table_ekf_output1['py_gt'],
    xaxis='x2',
    yaxis='y2',
    name = 'Ground Truth'
)
#Measurements
trace32 = Scatter(
    x=table_ekf_output2['px_gt'],
    y=table_ekf_output2['py_gt'],
    xaxis='x2',
    yaxis='y2',
    name = 'Ground Truth'
)

data1 = [trace11, trace21, trace31]
data2 = [trace12, trace22, trace32]

layout = Layout(
    xaxis2=dict(

        anchor='x2',
        title='px'
    ),
    yaxis2=dict(

        anchor='y2',
        title='py'
    )
)

fig1 = Figure(data=data1, layout=layout)
fig2 = Figure(data=data2, layout=layout)
py.plot(fig1, filename= 'EKF_data1.html')
py.plot(fig2, filename= 'EKF_data2.html')
