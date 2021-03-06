import base64
import io
import os
import argparse

import flask
from flask_login.utils import login_required # import dash_html_components as html

import dash
from dash import dcc # import dash_core_components as dcc
from dash import html
from dash import dash_table
from dash.dash_table.Format import Group # ?
from dash.dependencies import Input, Output, State
import dash_dangerously_set_inner_html as dhtml

import plotly.express as px
import plotly.graph_objs as go

import numpy as np
from sklearn.decomposition import PCA

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor

descdict = dict(Descriptors._descList)

app = dash.Dash(__name__)

css_directory = os.getcwd()+'/static/'
stylesheets = ['bWLwgP.css']
static_css_route = '/static/'

def smi2svg(smi):
    mol = Chem.MolFromSmiles(smi)
    rdDepictor.Compute2DCoords(mol)
    mc = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mc)
    drawer = Draw.MolDraw2DSVG(300,300)
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    return svg

def parse_contents(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    with open('./static/temp.sdf', 'w') as f:
        f.write(decoded.decode('utf-8'))
    mols = [mol for mol in Chem.SDMolSupplier('./static/temp.sdf')]
    return mols

@app.server.route('{}<stylesheet>'.format(static_css_route))
def serve_stylesheet(stylesheet):
    if stylesheet not in stylesheets:
        raise Exception()
    return flask.send_from_directory(css_directory, stylesheet)

for stylesheet in stylesheets:
    app.css.append_css({"external_url": "/static/{}".format(stylesheet)})


app.layout = html.Div(children=[
    html.H1(children='Hello Chemoinfo'),
    dcc.Upload(
        id='upload-data',
        children=html.Div(
            ['Drag and Drop SDF or ',
            html.A('Select Files')
            ]),
            style={
                    'width': '80%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
            },
            multiple=True
        ),
    
    html.Div(children='''
    Dash : sample plot
    '''),

    html.Div([dcc.Dropdown(id='x-column',
                           value='MaxEStateIndex',
                           options=[{'label': key, 'value': key} for key in descdict.keys()],
                           style={'width':'48%', 'display':'inline-block'}),
              dcc.Dropdown(id='y-column',
                           value='MaxEStateIndex',
                           options=[{'label': key, 'value': key} for key in descdict.keys()],
                           style={'width':'48%', 'display':'inline-block'}),
                           ]),
    html.Div([
        html.Div([html.Div(id="molimg")], className="four columns"),
        html.Div([dcc.Graph(id='example-graph')], className="eight columns")
        ], className="row"),
    #html.Div([dcc.Graph(id='chemical-space')])
    ])

@app.callback(
    Output('example-graph', 'figure'),
     
    [Input('upload-data', 'contents'),
     Input('x-column', 'value'),
     Input('y-column', 'value')]
)
def update_graph(contents, x_column_name, y_column_name):
    mols = parse_contents(contents[0])
    for i, mol in enumerate(mols):
        AllChem.Compute2DCoords(mol)
    x = [descdict[x_column_name](mol) for mol in mols]
    y = [descdict[y_column_name](mol) for mol in mols]

    return {'data':[go.Scatter(
        x=x,
        y=y,
        #text=['mol_{}'.format(i) for i in range(len(mols))],
        text=[Chem.MolToSmiles(mol) for mol in mols],
        mode='markers',
        marker={
            'size':15,
            'opacity':0.5
        }
    )],
    'layout':go.Layout(
        xaxis={'title':x_column_name},
        yaxis={'title':y_column_name}
    )}


@app.callback(
    Output('molimg', 'children'),
    [Input('example-graph', 'hoverData'),
    ]
)
def update_img(hoverData1):
    try:
        svg = smi2svg(hoverData1['points'][0]['text'])
    except:
        svg = 'Select molecule'
    return dhtml.DangerouslySetInnerHTML(svg)

if __name__=='__main__':
    app.run_server(debug=True)