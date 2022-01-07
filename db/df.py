import pandas as pd


def dataQuali():
    # Função que transpõe os dados da quali para o app
    dbQuali_source = 'https://github.com/edtrawtmam/HSPmol/raw/main/db/BaseDadosQuali.xlsx?raw=true'
    dfu = pd.read_excel(dbQuali_source)
    dfu.drop('Unnamed: 0', axis=1, inplace=True)
    return dfu





# Formatações gerais para os dataframes
format_semaforo = {
   'semaforo_NEV':'{0:,.0f}',
   'semaforo_4BZC':'{0:,.0f}',
   'semaforo_BZC':'{0:,.0f}',
   'semaforo_TEO':'{0:,.0f}',
   'semaforo_AS':'{0:,.0f}',
   'semaforo_CAF':'{0:,.0f}',
   'semaforo_ASC':'{0:,.0f}',
   'semaforo_ADPC':'{0:,.0f}',
   'semaforo_AAS':'{0:,.0f}',
   'dd':'{0:,.2f}',
   'dp':'{0:,.2f}',
   'dh':'{0:,.2f}',
   'dv':'{0:,.2f}',
   'dt':'{0:,.2f}',
   'Ra':'{0:,.2f}',
   'Ro':'{0:,.2f}',
   'RED':'{0:,.2f}',
   'mM':'{0:,.2f}',
   'Vm':'{0:,.2f}',
   'Rm':'{0:,.2f}',
#   'column2':'format2'
   }

def highlight_semaforo(val):
    #print(type(val))
    if val == 0:
      color = 'red' 
    elif val == 1:
      color = 'green'
    elif val == 2:
      color = 'yellow'
    elif (val != 0) & (val != 1) & (val != 2) & (type(val) != 'str'):
      color = 'grey'
    elif (type(val) == 'str'):
     color = 'black'
    else:
     color = 'black'
    background_color = '#000066'
    sai = 'color: '+str(color)+', background-color:'+str(background_color) 
    return 'color: %s' %color

# Configura a exportação de gráficos
config = {
  'toImageButtonOptions': {
    'format': 'svg', # one of png, svg, jpeg, webp
    'filename': 'custom_image',
    #'height': 900,
    #'width': 800,
    'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
  }
}

# Configura parâmetros comuns a todos os gráficos
# Define as cores das classes
classe_color = {"IFA": "blue",
                "solv": "orange",
                "cof": "purple",
                "pol": "Brown",
                }

metodo_color = {"Abbott": "green",
                "Mathieu": "brown",
                'nelder-mead':'darkviolet', 
                'COBYLA':'blue', 
                'SLSQP':'black', 
                'BFGS':'violet'
                }

# Define as cores do experimento:
solub_color = {'s': 'green', 
              'n': 'red',
              'Bib': 'blue',
              'Alvo': 'brown'
              }

semaforo_color = {'Não dissolve':'red',
                  'Dissolve':'green',
                  'Dissolve após calor':'yellow',
              }

outlier_color = {0:'pink',
                 1:'purple',
              }


def line_style_metodo(x):
  ''' Função para estilizar os círculos de Hansen em função do método.
      Entrada, string método
      Saída, lista[estilo da linha, cor da linha]
      '''
  x = x
  line_dash = 'solid'
  line_color = 'aliceblue'
  if (x == 'Hansen'):
    line_dash = 'solid'
    line_color = 'green'
  if (x == 'Abbott'):
    line_dash = 'solid'
    line_color = 'green'
  elif (x == 'Mathieu'):
    line_dash = 'solid'#"dot"
    line_color = "brown"
  elif (x == 'nelder-mead'):
    line_dash = 'solid'#"dash"
    line_color = "darkviolet"
  elif (x == 'COBYLA'):
    line_dash = 'solid'#"longdash"
    line_color = "blue"
  elif (x == 'SLSQP'):
    line_dash = 'solid'#"dashdot"
    line_color = "black"
  elif (x == 'BFGS'):
    line_dash = 'solid'#"longdashdot"
    line_color = "violet"
  return [line_dash, line_color]

# Cria uma biblioteca de simbolos para uso nos gráficos
from plotly.validators.scatter.marker import SymbolValidator

raw_symbols = SymbolValidator().values
namestems = []
namevariants = []
symbols = []
for i in range(0,len(raw_symbols),3):
    name = raw_symbols[i+2]
    symbols.append(raw_symbols[i])
    namestems.append(name.replace("-open", "").replace("-dot", ""))
    namevariants.append(name[len(namestems[-1]):])

# Aqui acabam as configurações aleatórias de exibição

