#
#
#
# Under Construction
#
#
#
import plotly
from plotly.offline import init_notebook_mode
import plotly.graph_objs as go
plotly.offline.init_notebook_mode(connected=True)

intensity = [0, 0.142857142857143, 0.285714285714286, 0.428571428571429, 0.571428571428571, 0.714285714285714, 0.85714257142857, 1] 

i = [7, 0, 0, 0, 4, 4, 2, 6, 4, 0, 3, 7]
j = [3, 4, 1, 2, 5, 6, 5, 5, 0, 1, 2, 2]
k = [0, 7, 2, 3, 6, 7, 1, 2, 5, 5, 7, 6]


cavitycolor = [[0, 'rgb(0, 100, 0)']*6]

coppercolor = [[0, 'rgb(139, 69, 19)'],
               [1, 'rgb(139, 69, 19)'],
               [2, 'rgb(139, 69, 19)'],
               [3, 'rgb(139, 69, 19)'],
               [4, 'rgb(139, 69, 19)'],
               [5, 'rgb(139, 69, 19)'],
               ]

ct = 0.05 # copper_thickness
pl = float(5) # patch length
pw = float(5) # patch width
fl = float(2) # feeder length
fw = float(2) # feeder width

# height
ch = float(1)

tl = pl + fl # total length

tfp = (pw / 2) + (fw / 2)  # top feeder point
bfp = (pw / 2) - (fw / 2)  # bottom feeder point

data = [
    go.Mesh3d(
        x = [0, 0, tl, tl, 0, 0, tl, tl],
        y = [0, pw, pw, 0, 0, pw, pw, 0],
        z = [0, 0, 0, 0 ] + ([ct]*4),
        colorbar = go.ColorBar(
            title='ground'
        ),
        facecolor = coppercolor,
        intensity = intensity,
        i = i,
        j = j,
        k = k,
        name = 'ground',
        showscale = True
    ),
    go.Mesh3d(
        x = [0, 0, pl, pl, 0, 0, pl, pl],
        y = [0, pl, pl, 0, 0, pl, pl, 0],
        z = ([ch]*4) + ([ch + ct]*4),
        colorbar = go.ColorBar(
            title='patch_top'
        ),
        facecolor = coppercolor,
        intensity = intensity,
        i = i,
        j = j,
        k = k,
        name = 'patch_top',
        showscale = True
    ),
    go.Mesh3d(
        x = [pl, pl, tl, tl, pl, pl, tl, tl],
        y = [tfp, bfp]*4,
        z = ([ch]*4) + ([ch+ct]*4),
        colorbar = go.ColorBar(
            title='feeder_top'
        ),
        facecolor= coppercolor,
        i = i,
        j = j,
        k = k,
        name = 'feeder_top',
        showscale = True
    ),
    go.Mesh3d(
        x = [0, 0, tl, tl, 0, 0, tl, tl],
        y = [0, pw, pw, 0, 0, pw, pw, 0],
        z = ([0 + ct] * 4) + ([ch] * 4),
        colorbar = go.ColorBar(
            title='cavity'
        ),
        facecolor = cavitycolor,
        intensity = intensity,
        i = i,
        j = j,
        k = k,
        name='cavity',
        showscale=True
    ),
    go.Mesh3d()
]

layout = go.Layout(
    xaxis=go.XAxis(
        title='x'
    ),
    yaxis=go.YAxis(
        title='y'
    )

)

fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig)