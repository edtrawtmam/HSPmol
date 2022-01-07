import dash
from dash import dcc # import dash_core_components as dcc
from dash import html
from dash import dash_table
from dash.dash_table.Format import Group # ?
from dash.dependencies import Input, Output, State
from flask_login.utils import login_required # import dash_html_components as html
import plotly.express as px
import pandas as pd
from db.df import *

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import PandasTools


dfu = dataQuali()

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
df_semaforo = pd.DataFrame()
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
dfr = dfu[(dfu['classe'].isin(['cof', 'IFA'])) & dfu['Metodo'].isin(metodos) & (dfu.dd.isna() == False)].copy()
dfr['molecula']=dfr['Nome']
# Insere as moléculas no dataframe do relatório
Chem.PandasTools.AddMoleculeColumnToFrame(dfr, smilesCol='Smile', molCol='molecula', includeFingerprints=True)
Chem.PandasTools.ChangeMoleculeRendering(dfr, renderer="SVG")
Chem.PandasTools.RenderImagesInAllDataFrames(images=True)
# Projeta em html o dataframe relatório
tabela_relatorio = dfr.set_index(['molecula', 'Nome', 'Metodo'])[['dd', 'dp','dh', 'Ro','dt', 'Vm','Eval','Iter']].style.format(format_semaforo).to_html()
#print(tabela_relatorio)

# Cria o aplicativo dash vinculado ao aplicativo flask
def creat_dash_application(flask_app):
    import dash_bootstrap_components as dbc

    dash_app = dash.Dash(
        server = flask_app,
        name = "Dashboard",
        url_base_pathname = "/dash/",
        external_stylesheets=[dbc.themes.BOOTSTRAP]
    )

    table_app = dash.Dash(
        server = flask_app,
        name = "Tabela",
        url_base_pathname = "/table/",
        external_stylesheets=[dbc.themes.BOOTSTRAP]
    )

    dash_app.layout = html.Div(children=[
        html.H1(children='Olá! Analise!'),

        html.Div(children='''
            HSPmol com Dash.
        '''),

    ])


    # Apresenta a tabela completa
    # Sorting operators (https://dash.plotly.com/datatable/filtering)
    # print(dfr.columns)
    table_app.layout = html.Div([
        html.Div([
            dcc.Input(
                id='adding-rows-name',
                placeholder='Nome da nova coluna',
                value='',
                style={'padding': 10}
                ),
            html.Button('Add Column', id='adding-columns-button', n_clicks=0)
            ], style={'height': 50}),

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
            hidden_columns = [],
            page_action="native",       # all data is passed to the table up-front or not ('none')
            page_current=0,             # page number that user is on
            page_size=15,                # number of rows visible per page
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
    html.Button('Exportar para o Excel', id='save_to_csv', n_clicks=0),

    # Create notification when saving to excel
    html.Div(id='placeholder', children=[]),
    dcc.Store(id="store", data=0),
    dcc.Interval(id='interval', interval=1000),

    dcc.Graph(id='my_graph')

    ])


    def update_styles(selected_columns):
        return [{
            'if': {'column_id': i},
            'background_color': '#D2F3FF'
        } for i in selected_columns]


    @table_app.callback(
        Output('proj-table', 'columns'),
        [Input('adding-columns-button', 'n_clicks')],
        [State('adding-rows-name', 'value'),
        State('proj-table', 'columns')],
    )
    def add_columns(n_clicks, value, existing_columns):
        print(existing_columns)
        if n_clicks > 0:
            existing_columns.append({
                'name': value, 'id': value,
                'renamable': True, 'deletable': True
            })
        print(existing_columns)
        return existing_columns


    @table_app.callback(
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


    @table_app.callback(
        Output('my_graph', 'figure'),
        [Input('proj-table', 'data')])

    def display_graph(data):
        df_fig = pd.DataFrame(data)
        # print(data)
        print(df_fig)
        fig = px.bar(df_fig, x='dv', y='dh')
        return fig

    @table_app.callback(
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


