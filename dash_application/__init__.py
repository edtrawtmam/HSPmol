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
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdDepictor
descdict = dict(Descriptors._descList)

import pandas as pd
from db.df import *

css_directory = os.getcwd()+'/static/'
stylesheets = ['bWLwgP.css']
static_css_route = '/static/'

external_stylesheets = [os.getcwd()+'/static/bWLwgP.css']

dfu = dataQuali()

#pd.set_option('precision', 2)
#pd.options.display.float_format = '${:.2f}'.format

# Criado o dataframe com os dados da quali, define-se variaveis de configuração de ambiente e aleatórias de exibição
# Passa-se a manipular os dataframes e criar as views
# Antes, define-se algumas variáveis auxiliares:
# Na Versão da função de otimização contando Outliers
# E na Versão otmização em lote
# Rotina de otimização multimétodos - tentativa de melhora nas funções
# Nomes dos semáforos/listas de resultados experimentais de solubilidade:
semaforos = ['semaforo_NEV', 'semaforo_4BZC', 'semaforo_BZC', 'semaforo_TEO', 'semaforo_AS', 'semaforo_ADPC', 'semaforo_AAS'] #'semaforo_CAF', 'semaforo_ASC',  Não possuem bons solventes
# Métodos de otimização disponíveis:
metodos =['nelder-mead', 'COBYLA', 'SLSQP', 'BFGS', 'Mathieu', 'Abbott']
metodos_o =['nelder-mead', 'COBYLA', 'SLSQP', 'BFGS']

# Referência Bibliográfica
# 'Abbott'
# Determinação por contribuição de grupos funcionais
# 'Mathieu'
# Variáveis de apoio
m = metodos[0]
aps=''
apelido = ''
repeticao=1 # Quantidade de vezes que recalcula para todas as 15 substâncias.
p=0 # Passo de repetição
# Base de dados de solventes utilizada na otimização
solution=[] # Vai armazenar os resultados das otimizações. [metodo][dd][dp][dh][Ro]

# Define algumas váriáveis contendo constantes:
R = 8.314462618 # J mol-1 K-1 (exatamente) - https://pt.wikipedia.org/wiki/Constante_universal_dos_gases_perfeitos
T_celsius = 20.0 # Temperatura em Célsius
T = T_celsius + 273.15 # Temperatura em Kelvin

# Inicialmente um sub-dataframe de relatório das moléculas alvo 
#dfr = dfu[(dfu['classe'].isin(['cof', 'IFA'])) & dfu['Metodo'].isin(metodos) & (dfu.dd.isna() == False)].copy()
#dfr['molecula']=dfr['Nome']
# Insere as moléculas no dataframe do relatório
#Chem.PandasTools.AddMoleculeColumnToFrame(dfr, smilesCol='Smile', molCol='molecula', includeFingerprints=True)
#Chem.PandasTools.ChangeMoleculeRendering(dfr, renderer="SVG")
#Chem.PandasTools.RenderImagesInAllDataFrames(images=True)
# Projeta em html o dataframe relatório
#tabela_relatorio = dfr.set_index(['molecula', 'Nome', 'Metodo'])[['dd', 'dp','dh', 'Ro','dt', 'Vm','Eval','Iter']].style.format(format_semaforo).to_html()
#print(tabela_relatorio)


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

# Cria o aplicativo dash vinculado ao aplicativo flask
def creat_dash_application(flask_app):
    import dash_bootstrap_components as dbc

    dash_app = dash.Dash(
        server = flask_app,
        name = "Tabela",
        url_base_pathname = "/dash/",
        external_stylesheets=[dbc.themes.BOOTSTRAP]
    )

    #@dash_app.server.route('{}<stylesheet>'.format(static_css_route))
    #def serve_stylesheet(stylesheet):
    #    if stylesheet not in stylesheets:
    #        raise Exception()
    #    return flask.send_from_directory(css_directory, stylesheet)
    #for stylesheet in stylesheets:
    #    dash_app.css.append_css({"external_url": "/static/{}".format(stylesheet)})


    dash_app.layout = html.Div([

        html.H1(children='HSPmol nascente'),
        dcc.Upload(
            id='upload-data',
            children=html.Div(
                ['Arraste e solte o SDF ou ',
                html.A('Selecione os Arquivos')
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
        
        dcc.Graph(id='my_graph'), 

        dash_table.DataTable(
                id='proj-table',
                columns=[
                    {"name": i, "id": i, "deletable": False, "selectable": True, "hideable": True}
                    if i == "CID"  or i == "CAS" or i == "classe" or i == "dd"  or i == "dv"  or i == "dp"  or i == "dh" or i == "dt" or i == "Nome" or i == "Metodo" or i == "Smile"
                    else {"name": i, "id": i, "deletable": True, "selectable": True, "hideable": True}
                    for i in dfu.columns
                ],
                data=dfu.to_dict('records'),  # the contents of the table
                editable=True,              # allow editing of data inside all cells
                filter_action="native",     # allow filtering of data by user ('native') or not ('none')
                sort_action="native",       # enables data to be sorted per-column by user or not ('none')
                sort_mode="multi",         # sort across 'multi' or 'single' columns
                column_selectable="multi",  # allow users to select 'multi' or 'single' columns
                row_selectable="multi",     # allow users to select 'multi' or 'single' rows
                row_deletable=True,         # choose if user can delete a row (True) or not (False)
                selected_columns=[],        # ids of columns that user selects
                selected_rows=[],           # indices of rows that user selects
                hidden_columns = ['CID', 'CAS', 'Obs', 'Smile', 'ITER', 'EVAL',
                                'Ra_BFGS', 'RED_BFGS', 'X12_BFGS', 'out_ruimDentro_BFGS', 'out_bomFora_BFGS', 'Iter_BFGS', 'Eval_BFGS',
                                'Ra_COBYLA', 'RED_COBYLA', 'X12_COBYLA', 'out_ruimDentro_COBYLA', 'out_bomFora_COBYLA', 'Iter_COBYLA', 'Eval_COBYLA',
                                'Ra_SLSQP', 'RED_SLSQP', 'X12_SLSQP', 'out_ruimDentro_SLSQP', 'out_bomFora_SLSQP', 'Iter_SLSQP', 'Eval_SLSQP',
                                'Ra_nelder-mead', 'RED_nelder-mead', 'X12_nelder-mead', 'out_ruimDentro_nelder-mead', 'out_bomFora_nelder-mead', 'Iter_nelder-mead', 'Eval_nelder-mead',
                                'Ra_Mathieu', 'RED_Mathieu', 'X12_Mathieu', 'out_ruimDentro_Mathieu', 'out_bomFora_Mathieu', 
                                'Ra_Abbott', 'RED_Abbott', 'X12_Abbott', 'out_ruimDentro_Abbott', 'out_bomFora_Abbott',
                                'Semaforo_TEO', 'Eval', 'Iter', 
                                'semaforo_NEV', 'semaforo_4BZC', 'semaforo_BZC', 'semaforo_TEO', 'semaforo_AS', 'semaforo_ADPC', 'semaforo_AAS'],

                page_action="native",       # all data is passed to the table up-front or not ('none')
                page_current=0,             # page number that user is on
                #page_size=15,                # number of rows visible per page
                fixed_rows={'headers': True},
                style_as_list_view=True,
                style_table={'height': 500},  # defaults to 500

                style_cell={                # ensure adequate header width when text is shorter than cell's text
                    'minWidth': 95, 'maxWidth': 95, 'width': 95
                },
                style_cell_conditional=[    # align text columns to left. By default they are aligned to right
                    {
                        'if': {'column_id': c},
                        'textAlign': 'center'
                    } for c in ['Apelido', 'classe']
                ],
                style_data={                # overflow cells' content into multiple lines
                    'whiteSpace': 'normal',
                    'height': 'auto'
                }
            ),

        html.Button('Nova linha', id='editing-rows-button', n_clicks=0),
        html.Button('Exportar .csv', id='save_to_csv', n_clicks=0),

        # Create notification when saving to excel
        html.Div(id='placeholder', children=[]),
        dcc.Store(id="store", data=0),
        dcc.Interval(id='interval', interval=1000),

        dcc.Input(
                id='adding-rows-name',
                placeholder='Nome da nova coluna',
                value='',
                style={'padding': 10}
                ),
        html.Button('Adicionar nova coluna', id='adding-columns-button', n_clicks=0)], style={'height': 50})

    @dash_app.callback(
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


    @dash_app.callback(
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

    @dash_app.callback(
            Output('proj-table', 'columns'),
            [Input('adding-columns-button', 'n_clicks')],
            [State('adding-rows-name', 'value'),
            State('proj-table', 'columns')],
        )
    def add_columns(n_clicks, value, existing_columns):
            #print(existing_columns)
            if n_clicks > 0:
                existing_columns.append({
                    'name': value, 'id': value,
                    'renamable': True, 'deletable': True
                })
            #print(existing_columns)
            return existing_columns

    @dash_app.callback(
            Output('proj-table', 'data'),
            [Input('editing-rows-button', 'n_clicks')],
            [State('proj-table', 'data'),
            State('proj-table', 'columns')],
        )
    def add_row(n_clicks, rows, columns):
            # print(rows)
            if n_clicks > 0:
                rows.append({c['id']: '' for c in columns})
            # print(rows)
            return rows

    @dash_app.callback(
            Output('my_graph', 'figure'),
            [Input('proj-table', 'data')])
    def display_graph(data):
            df_fig = pd.DataFrame(data)
            # print(data)
            #print(df_fig)
            fig = px.scatter(df_fig, x='dv', y='dh', hover_data=df_fig)
            return fig

    @dash_app.callback(
            [Output('placeholder', 'children'),
            Output("store", "data")],
            [Input('save_to_csv', 'n_clicks'),
            Input("interval", "n_intervals")],
            [State('proj-table', 'data'),
            State('store', 'data')]
        )
    def df_to_csv(n_clicks, n_intervals, dataset, s):
        output = html.Plaintext("Os dados serão salvos na sua pasta.",
                                    style={'color': 'green', 'font-weight': 'bold', 'font-size': 'large'})
        no_output = html.Plaintext("", style={'margin': "0px"})

        input_triggered = dash.callback_context.triggered[0]["prop_id"].split(".")[0]

        if input_triggered == "save_to_csv":
                s = 6
                df = pd.DataFrame(dataset)
                df.to_csv("ProjetoAnaliseHSP.csv")
                return output, s
        elif input_triggered == 'interval' and s > 0:
                s = s-1
                if s > 0:
                    return output, s
                else:
                    return no_output, s
        elif s == 0:
                return no_output, s

    # Autenticação:
    for view_function in dash_app.server.view_functions:
        if view_function.startswith(dash_app.config.url_base_pathname):
                dash_app.server.view_functions[view_function] = login_required(dash_app.server.view_functions[view_function])
    
        return dash_app
