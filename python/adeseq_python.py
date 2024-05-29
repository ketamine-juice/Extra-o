#%% md
# # Bibliotecas
#%%
"""

------------------------------------------------------------------------
#                                                                      #
#    Secção de verificação, importação e instalação de bibliotecas.    #
#                                                                      #
------------------------------------------------------------------------

"""
import os
import pip

def import_or_install(package):
    """
    Tenta importar um pacote. Se não for possível importá-lo, tenta instalar
    o pacote usando o pip e, em seguida, importá-lo novamente.

    Parâmetros:
        package (str): O nome do pacote a ser importado ou instalado.

    Exemplo:
        >>> import_or_install('numpy')
        Tenta importar o pacote 'numpy'. Se não estiver disponível, tenta
        instalar usando o pip e, em seguida, importá-lo novamente.
    """
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package, '--quiet'])
        


# Ir populando com packages à medida que vão sendo acrescentadas ao projeto.

pkgs = [
    'pandas',         # Biblioteca para manipulação e análise de dados em Python
    'seaborn',        # Biblioteca para criação de gráficos estatísticos mais atraentes e informativos
    'matplotlib',     # Biblioteca para criação de gráficos estáticos e interativos em Python
    'scipy',          # Biblioteca para funções e algoritmos matemáticos e estatísticos em Python
    'numpy',          # Biblioteca para manipulação de arrays e matrizes multidimensionais em Python
    'scikit-learn',   # Biblioteca para aprendizado de máquina e mineração de dados em Python
    'pydeseq2',       # Interface Python para o pacote DESeq2, usado para análise de expressão diferencial
    'gseapy',         # Biblioteca para análise de enriquecimento de genes
    'sanbomics',      # Biblioteca para análise bioinformática em Python
    'fastcluster'     # Biblioteca para agrupamento rápido de grandes conjuntos de dados
]

for pkg in pkgs:
    import_or_install(pkg)
    
import os
import pandas
import sklearn
import numpy
import seaborn
import matplotlib
import scipy
import gseapy
import pydeseq2
import sanbomics
#%% md
# # Upload dos Dados
#%% md
# Os dados para esta parte de python foram obtidos por conversão das variáveis em R. Usamos os packages TCGA BioLinks para transferir os dados necessários e depois estes foram exportados.
# 
# A lógica em R:
# 
# ```r
# write_csv(seqdata, 'seqdata.csv') # Para excluir os rownames e manter os nomes do gene
# write_csv(meta, 'meta.csv', rownames = FALSE) # Para manter o barcode
# ```
#%%
print('A criar pasta para output')

if not os.path.exists('./output'):
    os.mkdir('./output')    
#%% md
# # Pré-processamento de dados
#%% md
# Abaixo é realizado o pré-processamento dos dados importados. É utilizada uma pipeline comum com a biblioteca `pandas`, sendo realizadas várias operações de modificação de índices, nomes de colunas e filtragem de dados. O resultado final é um *dataframe* com 248 colunas, correspondente às amostras com subtipo de adenocarcinoma do pulmão e 60660 linhas correspondentes a todos os genes analisados. Os passos seguintes consistirão na filtragem dos dados com base em CPM (*counts per million*).
#%%
floc = './dataset/seqdata.csv'
print(f'A importar o dataframe na localização {floc}\n')
seqdata = pandas.read_csv('./dataset/seqdata.csv')
print('A renomear colunas e a definir genes como índice.\n')
seqdata = seqdata.rename(columns={'Unnamed: 0':'Genes'})
seqdata = seqdata.set_index('Genes')
print('Concluído!',
      'Dataframe seqdata com:',
      f'    - {seqdata.shape[1]} colunas',
      f'    - {len(seqdata)} linhas',
      sep= '\n'
      )
print('A mostrar head do dataframe\n',seqdata.head())
#%% md
# # Metadados
#%% md
# ## Carregamento dos metadados 1.0
#%%
# Carregamento de metadados
floc = './dataset/meta.csv'
print(f'A importar o dataframe na localização {floc}\n')
meta = pandas.read_csv('./dataset/meta.csv')
print('A renomear colunas e a definir barcode como índice.\n')
meta = meta.set_index('barcode')
print('Concluído!',
      'Dataframe seqdata com:',
      f'    - {meta.shape[1]} colunas',
      f'    - {len(meta)} linhas',
      sep= '\n'
      )
print('A mostrar head do dataframe\n',meta.head())
#%% md
# ### Verificação de valores omissos
#%%
# Verifica a existência de dados omissos
print('A verificar a existência de dados omissos:', meta.isna().any(), sep='\n')
#%%
# Remover NAs do expr. subtype
print('A remover dados omissos do dataframe de metadados.')
meta_clean = meta.dropna(how='any')

if len(meta_clean) == 248:
    print(f'data set de subtipo de expressão contém 248 valores como esperado, com tamanho {len(meta_clean)}')
    
print('Concluído!',
      'A mostrar head do dataframe tratado:',
      meta_clean.head(),
      sep='\n')
#%% md
# ### Filtragem das colunas do dataframe seqdata de acordo com os metadados existentes
#%%
# Filtrar seqdata com base no filtro
print('A reorganizar e filtrar dataframe seqdata com base nos índices dos metadados.')
seqdata_clean = seqdata[meta_clean.index] # Não é preciso validar pois as colunas já estão pelo nome
print('Concluído!')
#%% md
# ## Carregamento dos metadados 2.0
#%% md
# Foi feita um segundo objeto com metadados relacionados com os hábitos de fumo e com o estado vital.
#%% md
# ### Hábitos de fumo
#%%
# Carregamento de hábitos de fumo
floc = './dataset/fumo.csv'
print(f'A importar o dataframe na localização {floc}\n')
smoke_status = pandas.read_csv('./dataset/fumo.csv')
print('A renomear colunas.\n')
smoke_status = smoke_status.rename(columns={
    'data.paper_Smoking.Status' : 'Smoke Habits',
    'data.years_smoked' : 'Years Smoked'})
print('A corrigir índice de acordo com o índice do dataframe de metadados.')
smoke_status.drop(columns=['Unnamed: 0', 'Years Smoked'], inplace=True)
smoke_status['barcode'] = meta.index

# Filtramos por Smoke Habits
print('A limpar dados:')
print('Removendo dados omissos.')
smoke_status.isna().any()
smoke_status_clean = smoke_status[(smoke_status['Smoke Habits'].notna()) & (smoke_status['Smoke Habits'] != '[Not Available]')]
print('A gerar dataset limpo')
smoke_status_clean = smoke_status_clean[['barcode', 'Smoke Habits']]
print('Concluído!',
      'Dataframe seqdata com:',
      f'    - {smoke_status_clean.shape[1]} colunas',
      f'    - {len(smoke_status_clean)} linhas',
      sep= '\n'
      )
print('A mostrar head do dataframe\n',smoke_status_clean.head())


#%% md
# ### Estado vital
#%%
vital = pandas.read_csv('./dataset/vital.csv')

vital.rename(columns={'data.barcode' : 'barcode', 'data.vital_status' : 'Vital Status'}, inplace=True)

vital_vs_smoke = pandas.merge(
    smoke_status_clean,
    vital,
    on = 'barcode',
    how = 'inner'
)

vital_vs_smoke = vital_vs_smoke[['barcode', 'Smoke Habits', 'Vital Status']]
vital_vs_smoke['Smoke Habits'] = vital_vs_smoke['Smoke Habits'].replace(
    {'Current reformed smoker for > 15 years' : 'EX > 15 y',
    'Current reformed smoker for < or = 15 years' : 'EX <= 15 y',
    'Current smoker': 'Smokes',
    'Lifelong Non-smoker': 'Never smoked'
})
vital_vs_smoke
#%% md
# ## Carregamento dos metadados 3.0
#%%
stage = pandas.read_csv('./dataset/stage.csv')
stage.drop(columns=['Unnamed: 0'], inplace=True)

stage = stage.rename(columns={'data$paper_Tumor.stage' : 'Tumor Stage'})

stage['barcode'] = meta.index

print(f'encontrados {stage.iloc[:, 0].isna().sum()} valores omissos',
           'foi criada uma versão filtrada')

stage_clean = stage.dropna(how='any')
stage_clean[stage_clean['barcode'].isin(meta.index)]
stage_clean
#%% md
# ## Conjuntos
#%% md
# No final, foram criados três conjuntos de metadados:
# 
# * Conjunto 1 (**meta**): género e tipo de expressão
# 
# * Conjunto 2 (**vital_vs_smoke**): hábitos de fumo e estado vital
# 
# * Conjunto 3 (**stage_clean**): estágio de tumor
#%% md
# # Análise inicial dos dados
#%% md
# À semelhança do trabalho realizado em R, são realizadas algumas análise estatísticas básicas com o intuito do correlacionar dados, nomeadamente:
# * Distribuição de estado vital;
# * Distribuição de hábitos de fumo;
# * Distribuição de sexo biológico;
# * Correlação Hábitos de fumo / Estado vital;
# 
# Os resultados obtidos e os gráficos demonstrados são iguais aos demonstrados ao longo do relatório desenvolvido utilizando R, sendo a principal diferença o tipo de bibliotecas utilizadas.
#%% md
# ## Estado vital
#%%
print('A preparar tabela de contingência para estado vital.')
vital_status = vital['Vital Status'].value_counts()

# Plotting the bar plot
print('A preparar o gráfico.')
seaborn.barplot(x=vital_status.index, y=vital_status.values, palette="viridis")

matplotlib.pyplot.title("Gráfico de Barras do Status Vital")
matplotlib.pyplot.xlabel("Status Vital")
matplotlib.pyplot.ylabel("Contagem")


matplotlib.pyplot.savefig('./output/vital_status.pdf')
# Showing the plot
print('Concluído, a exibir gráfico abaixo.', 'Plot gravado em ./output\n',sep='\n')
matplotlib.pyplot.show()



#%% md
# ## Hábitos de fumo
#%%
# Adjusting smoke habit data
print('A preparar tabela de contingência para hábitos de fumo.')
smoke_habit = smoke_status_clean['Smoke Habits'].value_counts()
smoke_habit = smoke_habit.rename({
    'Current reformed smoker for > 15 years': 'EX > 15 y',
    'Current reformed smoker for < or = 15 years': 'EX <= 15 y',
    'Current smoker': 'Smokes',
    'Lifelong Non-smoker': 'Never smoked',
})

print('A preparar o gráfico.')
# Plotting smoke habit data
seaborn.barplot(x=smoke_habit.index, y=smoke_habit.values, palette="mako")

matplotlib.pyplot.title("Gráfico de Barras do Hábito de Fumar")
matplotlib.pyplot.xlabel("Hábito de Fumar")
matplotlib.pyplot.ylabel("Contagem")

matplotlib.pyplot.savefig('./output/smoke_habit.pdf')
# Showing the plot
print('Concluído, a exibir gráfico abaixo.', 'Plot gravado em ./output/\n',sep='\n')
matplotlib.pyplot.show()
#%% md
# ## Hábitos de fumo vs Estado vital
#%%
print('Correlacionar hábitos de fumo com estado vital através de teste qui-quidrado.',
      'A preparar tabela de contingência para alfa = 0.05',
      sep='\n')
contingency_table = pandas.crosstab(vital_vs_smoke['Vital Status'], vital_vs_smoke['Smoke Habits'])
alpha = 0.05
print('A realizar teste qui-quadrado.')
chi2, p, dof, expected = scipy.stats.chi2_contingency(contingency_table)
print('Concluído, resultados abaixo:',
      f'Estatística de teste:{chi2}',
      f'Graus de liberdade:{dof}',
      f'Valor de prova:{p}',
      f'Para alfa = {0.05}:',
      sep = '\n')

# para 95% IC
if p < 0.05:
    print('É plausível considerar correlação')
    
else:
    print('Não é plausível considerar correlação')

print('A criar um gráfico para os resultados.')
seaborn.countplot(x='Smoke Habits', hue = 'Vital Status', data = vital_vs_smoke)
matplotlib.pyplot.title('Distribuição de Estado Vital vs Hábitos de Fumo')
matplotlib.pyplot.xlabel('Hábitos de Fumo')
matplotlib.pyplot.ylabel('Contagem')
matplotlib.pyplot.legend(title = 'Estado Vital')

matplotlib.pyplot.savefig('./output/vital_vs_smoke.pdf')
print('Concluído, a exibir gráfico abaixo.', 'Plot gravado em ./output/\n',sep='\n')
matplotlib.pyplot.show()


    
#%% md
# # Análise de Expressão Diferencial utilizando PyDESEQ2
# 
# Para a análise de expressão genética diferencial em Python utilizamos a biblioteca `PyDESEQ2`, uma adaptação da biblioteca R DESEQ2. 
#%%
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

SAVE = False  # whether to save the outputs of this notebook

if SAVE:
    # Replace this with the path to directory where you would like results to be
    # saved
    OUTPUT_PATH = "../output_files/synthetic_example"
    os.makedirs(OUTPUT_PATH, exist_ok=True)  # Create path if it doesn't exist
#%%
# Pipeline
# https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_pandas_io_example.html
# https://pydeseq2.readthedocs.io/en/latest/auto_examples/plot_step_by_step.html
# Credits to mousepixels: https://github.com/mousepixels/sanbomics_scripts/blob/main/PyDeseq2_DE_tutorial.ipynb
# And PyDeSEQ2 documentation

import pickle
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

#%%
print('A preparar dados para DESEQ:,'
      ' - Transport tabela de counts para cada gene',
      sep='\n')
counts = seqdata_clean
counts_df = counts.T
print('A criar tabela de metadados e verificar se todos os dados correspondem corretamente.')
metadata = meta_clean[['expr']]
# Verificar se dados correspondem
if (counts_df.index == meta_clean.index).all():
    print('Dados correspondem corretamente, a filtrar dataset retendo genes com counts > 10')
    genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
    print('A converter em dataframe')
    counts_df = counts_df[genes_to_keep]
#%%
import pydeseq2

inference = pydeseq2.default_inference.DefaultInference(n_cpus=8)
print('A gerar objeto DeseqDataSet')
dds = pydeseq2.dds.DeseqDataSet(
    counts = counts_df,
    metadata = metadata,
    design_factors='expr',
    refit_cooks=True,
    inference = inference
)

if type(dds) == pydeseq2.dds.DeseqDataSet:
    print('Dataset criado, a correr operações DESEQ:')
    dds.deseq2()

#%%
# Guardar para acelerar o processo mais á frente com o load (estilo RDA/RDS)
# open("dds.pkl", "wb") as f:
#    pickle.dump(dds, f)
    
# 
#with open('dds.pkl', 'rb') as f:
#    pickle.load(f)
#%% md
# ## Testes de Wald
# 
# Na célula abaixo é feito o tratamento estatístico dos resultados, sendo realizados testes de Wald para comparar a expressão genética de diferentes genes analisados com o objetivo de verificar e comparar se existem genes diferencialmente expressos entre as diferentes amostras. Neste caso concreto são comparadas as condições Proximal Proliferative (Prof. Prolif) vs Terminal Respiratory Unit (TRU) e Proximal Inflamatory (Prox. Inflam) vs Terminal Respiratory Unit (TRU). 
# 
# O resultado da análise consiste numa tabela com:
#  * os níveis médios da expressão genética;
#   * a taxa de sobre/sobexpressão sob a forma $log_{2}(\frac{Expressão na condição A}{Expressão na condição B}$;
#   * resultados do teste estatístico, nomeadamente os valores de prova;
# 
#%%
# Crie um objeto DeseqStats usando o conjunto de dados DESeq fornecido (dds) e os parâmetros especificados
# alpha: nível de significância para o teste (padrão 0.05)
# cooks_filter: se aplicar ou não a filtragem de distância de Cook para remover possíveis outliers
# independent_filter: se aplicar ou não a filtragem independente para aumentar o poder, removendo genes com contagem baixa
# inference: especifica o tipo de inferência estatística a ser usada (por exemplo, 'paramétrica', 'não paramétrica')
print('A executar tratamento estatístico dos dados.')
stat_res = pydeseq2.ds.DeseqStats(dds, alpha=0.05, cooks_filter=True, independent_filter=True, inference=inference)

# Execute o teste de Wald para análise de expressão diferencial no conjunto de dados DESeq
# O teste de Wald é um teste paramétrico usado para determinar se há uma diferença significativa na expressão
# Ele calcula valores de p e log fold changes para cada gene no conjunto de dados
print('A executar testes de Wald.')
stat_res.run_wald_test()

# Recupere os resultados do teste de Wald como um DataFrame do pandas
# Este DataFrame contém os resultados da análise de expressão diferencial, incluindo:
# - log fold changes
# - valores de p
# - valores de p ajustados (corrigidos para testes múltiplos)
# - estatísticas adicionais e metadados para cada gene
print('Concluído, a exibir resultados:')
stat_res.summary()


#%% md
# Em seguida filtramos os dados, definindo como filtro:
#  * valor de prova seja tal que seja plausível considerar que são diferencialmente expressos (para um intervalo de confiança de 95%);
#   * valor de $| log_{2}FC |$ superior a  2, indicativo de um nível de sobre/sobexpressão
#   * CPM superior a 20
#%%
#MA Plot
print('A criar MA Plot .')
matplotlib.pyplot.figure()

stat_res.plot_MA()
matplotlib.pyplot.axhline(y=2, color='green', linestyle='--')
matplotlib.pyplot.axhline(y=-2, color='green', linestyle='--')
print('MA Plot gravado em ./output/')
matplotlib.pyplot.savefig('./output/MA_Plot.pdf')
#%% md
# Ao analisar o gráfico acima vemos que existem genes diferencialmente expressos, sinalizados a vermelho (valores cujo valor de prova é inferior a 0.05). No gráfico optamos por implementar um thresold para valores de log2FC de 2 e -2 respetivamente como forma de destacar os genes com alteração mais significativa de expressão. É possível visualizar que existem alguns genes com valores log2FC discrepantes, que poderão ser outliers. 
#%%
# Extracting the results DataFrame from the statistical results object
print('A isolar dataframe de resultados')
rdf = stat_res.results_df

# Print the initial number of genes present in the results DataFrame
# This provides a count of the total number of rows (genes) before filtering
print('Genes inicialmente presentes -->', len(rdf))  # Expected to be > 52310

# Apply filtering criteria to select significant genes:
# - Adjusted p-value (padj) < 0.05
# - Absolute value of log2 fold change (|log2FoldChange|) > 2
# - Base mean expression (baseMean) > 20
significant = rdf[(rdf.padj < 0.05) & (abs(rdf.log2FoldChange) > 2) & (rdf.baseMean > 20)]

# Print the number of genes that meet the differential expression criteria
# This gives an idea of the subset size that meets the specified thresholds
print(f'Genes mais diferencialmente expressos (log2FoldChange > 2, CPM > 20, padj < 0.05 -->', len(significant))  

# Sort the significant genes by the adjusted p-value (padj) in ascending order
# This ensures that the most significant genes are at the top of the DataFrame
significant = significant.sort_values(by='padj')

# Display the top 100 genes based on the sorted DataFrame
# This returns the first 100 rows (genes) after sorting
significant_top_100 = significant.head(100)

# Print the indexes of the top 100 differentially expressed genes
# This provides a list of the most significant genes based on the filtering and sorting criteria
#print('Lista dos 100 genes mais diferencialmente expressos:', *significant_top_100.index, sep='\n')
#%%
# Isolar os 20 genes mais significativos
# Compute log-transformed counts and store them in a new layer 'logcounts'
# np.log1p(x) computes log(1 + x) which helps in managing zeros in the data
dds.layers['logcounts'] = numpy.log1p(dds.layers['normed_counts'])

# Display the 'logcounts' layer to verify the transformation
dds.layers['logcounts']

# Subset the DESeqDataSet to include only the first 50 rows and the significant indices
dds_sigs = dds[:, significant.index[:20,]]

# Create a DataFrame for graphing with the log-transformed counts
# The DataFrame's index is the variable names (var_names)
# The DataFrame's columns are the observation names (obs_names)
grapher = pandas.DataFrame(
    dds_sigs.layers['logcounts'].T,
    index=dds_sigs.var_names,
    columns=dds_sigs.obs
)

# Plot a clustered heatmap of the data using seaborn's clustermap function
# z_score=0 standardizes the rows (genes) by subtracting the mean and dividing by the standard deviation
# cmap='RdYlBu_r' sets the color map to a diverging palette from red to yellow to blue, reversed
clustermap = seaborn.clustermap(grapher, z_score=0, cmap='RdYlBu_r')

clustermap.savefig('./output/heatmap.pdf')
#%% md
# 
#%%
# Isolar os 20 genes mais significativos
# Compute log-transformed counts and store them in a new layer 'logcounts'
# np.log1p(x) computes log(1 + x) which helps in managing zeros in the data
dds.layers['logcounts'] = numpy.log1p(dds.layers['normed_counts'])

# Display the 'logcounts' layer to verify the transformation
dds.layers['logcounts']

# Subset the DESeqDataSet to include only the first 50 rows and the significant indices
dds_sigs = dds[:30, significant.index[:20,]]

# Create a DataFrame for graphing with the log-transformed counts
# The DataFrame's index is the variable names (var_names)
# The DataFrame's columns are the observation names (obs_names)
grapher = pandas.DataFrame(
    dds_sigs.layers['logcounts'].T,
    index=dds_sigs.var_names,
    columns=dds_sigs.obs
)

# Plot a clustered heatmap of the data using seaborn's clustermap function
# z_score=0 standardizes the rows (genes) by subtracting the mean and dividing by the standard deviation
# cmap='RdYlBu_r' sets the color map to a diverging palette from red to yellow to blue, reversed
clustermap = seaborn.clustermap(grapher, z_score=0, cmap='RdYlBu_r')

clustermap.savefig('./output/heatmap.pdf')
#%% md
# 
#%% md
# Ao analisar o heatmap acima, conseguimos identificar a existência de possíveis padrões para os genes apresentados, havendo uma clara subexpressão nos genes FLJ22447, UBE2C, ANLN, LINC00973, NFE4, SBSN, SPOCD1, SPINK4, HMGCS2, REG4, TFF2, DPCR1, CASR, RPL26P30, PCSK2, CYP4B1, MACROD2, PLA2G4F, PGC, CACNA2D2. 
# 
# Os genes FLJ22447, UBE2C, ANLN e LINC00973 são genes que têm sido associados aos fatores de prognóstico de doentes com LUAD. Estes genes estão envolvidos na regulação do ciclo celular, na divisão celular e na proliferação celular, e pensa-se que a sua sobreexpressão contribui para a natureza agressiva dos tumores LUAD. Os genes NFE4, SBSN, SPOCD1, SPINK4 e HMGCS2 são genes cujos papeis na doença são menos bem compreendidos, podendo estar envolvidos no processo metastático e supressão tumoral. Estes genes podem estar envolvidos em processos como a diferenciação celular, o metabolismo e a resposta ao stress. Por outro lado, REG4, TFF2, DPCR1, CASR, RPL26P30, PCSK2, CYP4B1, MACROD2, PLA2G4F e PGC são genes associados a uma função pulmonar normal e pensa-se que são genes supressores de tumores. A sua desregulação pode contribuir para o desenvolvimento e progressão do LUAD. O CACNA2D2 é um gene que foi identificado como um potencial marcador de prognóstico no LUAD. A sua expressão está associada a uma melhor sobrevivência global e a uma sobrevivência sem recidivas em doentes com LUAD. Globalmente, estes genes representam um conjunto diversificado de intervenientes moleculares no LUAD, estando alguns associados a uma doença mais agressiva e outros a um melhor prognóstico. 
# 
# Embora o quadro acima não permita fazer propriamente fazer uma distinção entre as diferentes condições (pois é possível ver que há alguns perfis semelhantes para condições diferentes) é possível perceber os perfis genéticos de cada amostra. Em todo o caso há um possível padrão e diferenças para as condições proximal inflamatory / proximal proliferative a TRU. É possível ainda ver que há efetivamente alteração nos perfis de expressão de genes reguladores críticos, desde genes supressores tumorais a genes relacionados com a agressividade e gravidade do tumor. Como próximos passos seria importante fazer uma análise mais detalhada para cada amostra e cada gene, constituindo estes potenciais alvos terapêuticos.
# 
# 
# Fontes:
# 
# - [LINC03033 long intergenic non-protein coding RNA 3033 Homo sapiens (human) - Gene - NCBI](https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=400221)
# - [Expression of UBE2C in lung adenocarcinoma based on database analysis and its clinical significance](https://pubmed.ncbi.nlm.nih.gov/33051417/)
# - [Comprehensive molecular profiling of lung adenocarcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4231481/)
# - [Prognostic significance of ANLN in lung adenocarcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6036425/)
# - [LINC00973 Induces Proliferation Arrest of Drug-Treated Cancer Cells by Preventing p21 Degradation - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7664178/)
# - [Identification of candidate genes associated with triple negative breast cancer - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5620011/)
# - [SBSN drives bladder cancer metastasis via EGFR/SRC/STAT3 signalling](https://pubmed.ncbi.nlm.nih.gov/35484216/)
# - [SPOCD1 promotes the proliferation and metastasis of glioma cells by up-regulating PTX3 - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5934553/)
# - [SPINKs in Tumors: Potential Therapeutic Targets - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8873584/)
# - [3158 - Gene ResultHMGCS2 3-hydroxy-3-methylglutaryl-CoA synthase 2](https://www.ncbi.nlm.nih.gov/gene/3158)
# - [REG4 is an indicator for KRAS mutant lung adenocarcinoma with TTF-1 low expression](https://pubmed.ncbi.nlm.nih.gov/31428934/)
# - [Diagnostic utility of trefoil factor families for the early detection of lung cancer and their correlation with tissue expression - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9996639/)
# - [Differences in the early stage gene expression profiles of lung adenocarcinoma and lung squamous cell carcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6865721/)
# - [CYP4B1 is a prognostic biomarker and potential therapeutic target in lung adenocarcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7886130/)
# - [MACROD2 Haploinsufficiency Impairs Catalytic Activity of PARP1 and Promotes Chromosome Instability and Growth of Intestinal Tumors](https://pubmed.ncbi.nlm.nih.gov/29880585/)
# - [Integrated analysis of gene expression and DNA methylation datasets identified key genes and a 6-gene prognostic signature for primary lung adenocarcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8596225/)
# - [Expression and Prognostic Significance of m6A-Related Genes in Lung Adenocarcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7049251/)
# - [Prognostic value of immune-related genes in the tumor microenvironment of lung adenocarcinoma and lung squamous cell carcinoma - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7138544/)
# 
# ## Análise de enriquecimento
# 
# Abaixo foi desenvolvida a análise de enriquecimento utilizando a biblioteca GSEApy. 
# 
# A análise de enriquecimento de conjuntos de genes (GSEA, do inglês Gene Set Enrichment Analysis) é uma abordagem estatística utilizada em genómica funcional para identificar padrões biológicos e processos celulares associados a conjuntos de genes específicos. A GSEA avalia se um conjunto predefinido de genes, agrupados com base em funções biológicas, vias metabólicas ou outras características, está significativamente enriquecido em uma amostra experimental em comparação com amostras de controle. Neste caso são comparadas as condições PP/TRU e PI/TRU. 
#%%
print('A mostrar a lista de bibliotecas de genes.')
gseapy.get_library_name()
#%%
lib = 'MSigDB_Oncogenic_Signatures'
print('A executar análise de enriquecimento e atribuição de rankings com base nos valores.',
      f'Irá ser utilizada a biblioteca {lib} por conter os oncogenes mais comuns.', sep='\n')
ranking = rdf[['stat']].dropna().sort_values('stat', ascending=False)

pre_res = gseapy.prerank(rnk = ranking,
                     gene_sets= ['MSigDB_Oncogenic_Signatures'],
                     seed = 6, permutation_num=100
                     )
print(
    'Concluído!',
    'A mostrar resultados de enriquecimento: ',
    '   - 3 genes com maior enrichment score (positivamente associados\n',
    pre_res.res2d[['Term','NES']].head(3),
    '   \n- 3 genes com menor enrichment score (negativamente associados)\n',
    pre_res.res2d[['Term','NES']].tail(3),
    sep='\n'
)


#%% md
# Abaixo são construídos os diferentes GSEA plots para os top 3 / bottom 3 genes.
# 
# Os GSEA plots, ou gráficos de análise de enriquecimento de conjuntos de genes, são uma representação visual essencial na interpretação dos resultados da análise de enriquecimento GSEA. Esses gráficos exibem a distribuição acumulada de genes dentro de conjuntos funcionais pré-definidos ao longo de uma lista ordenada de genes, geralmente classificados por sua relevância em uma comparação de condições experimentais. Iremos analisar e identificar picos ou vales significativos na curva de enriquecimento, indicando conjuntos de genes que estão sobre-representados ou sub-representados para as condições analisadas. 
#%%
terms = pre_res.res2d.Term
top_bottom_3 = [0,1,2,len(terms) - 1, len(terms) - 2, len(terms) - 3]
print(f'A mostrar plot para os termos:',*top_bottom_3, '(Top 3 / Bottom 3)')
i = 0
import os
for term in top_bottom_3:
    if not os.path.exists('./output/gsea_plots'):
        os.makedirs('./output/gsea_plots')
    gseapy.gseaplot(rank_metric=pre_res.ranking, term=terms[term], **pre_res.results[terms[term]])
    print(f'A gravar plot para {terms[term]}')
    gseapy.gseaplot(rank_metric=pre_res.ranking, term=terms[term], ofname=f'./output/gsea_plots/{terms[term]}.pdf', **pre_res.results[terms[term]])
    i += 1
#%% md
# Ao analisar os resultados obtidos, podemos concluir que os conjuntos de genes no "top 3" - CSR LATE UP, VEGFA UP e RB P107 UP - são upregulated estando estes conjuntos de genes relacionados com a modulação da resposta imune por parte do tumor, afetando a sua proliferação e resistência à eliminação e, consequente, o prognóstico. Os conjuntos de genes no "bottom 3" - GLI1 UP UP, DCA UP DN e RAPA EARLY  UP DN - apresentam uma variabilidade nos seus níveis de expressão, sendo upregulated em alguns casos e downregulated noutros - estes resultados podem dever-se a diferentes padrões e perfis de genéticos das diferentes condições, tornando estes genes significativos para o estudo e comparação dos diferente subtipos de tumor. Estes resultados fazem sentido uma vez que estes conjuntos de genes relacionam-se com fatores fenotípicos como o tipo de tecido tumoral, inibição enzimática (como o caso da inibição da pyruvate dehydrogenase kinase com DCA) e reguladores de crescimento tumoral como o caso de RAPA EARLY. 
# 
# Fontes:
# - [How to interpret GSEA results and plot](https://www.youtube.com/watch?v=Yi4d7JIlAsM)
# - [Identification of a novel intermittent hypoxia-related prognostic lncRNA signature and the ceRNA of lncRNA GSEC/miR-873-3p/EGLN3 regulatory axis in lung adenocarcinoma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10573295/)
# - [Increased VEGF-A in solid type of lung adenocarcinoma reduces the patients’ survival](https://www.nature.com/articles/s41598-020-79907-6)
# - [GLI1-Amplifications Expands the Spectrum of Soft Tissue Neoplasms Defined by GLI1 Gene Fusions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6821565/)
# - [Human Gene Set: DCA_UP.V1_DN](https://www.gsea-msigdb.org/gsea/msigdb/cards/DCA_UP.V1_DN)
# - [Human Gene Set: RAPA_EARLY_UP.V1_UP](https://www.gsea-msigdb.org/gsea/msigdb/cards/RAPA_EARLY_UP.V1_UP)
# 
# 
# 
#%% md
# Abaixo é gerado o volcano plot para os resultados.
#%%
from sanbomics.plots import volcano 
print('A criar volcano plot')
rdf_copy = rdf.copy()
rdf_copy['Symbol'] = rdf_copy.index
matplotlib.pyplot.figure()
volcano(rdf_copy, symbol='Symbol')
print('Concluído, a exibir plot. Plot gravado em ./output/')
matplotlib.pyplot.savefig('./output/volcano_plot.pdf')
matplotlib.pyplot.show()
#%% md
# Acima vemos o volcano plot para a visualização da análise de enriquecimento de forma mais fácil. É possível identificar no 2º quadrante do referencial os genes downregulated e no 1º quadrante os genes upregulated, comparando as condições PP/TRU e PI/TRU. É possível ainda identificar alguns genes de particular relevância, sinalizados no gráfico, que representam genes particularmente significativos no que respeita aos seus padrões de expressão conforme temos vindo a falar sobre eles. 
# 
# https://www.youtube.com/watch?v=7aWAdw2jhj0
#%% md
# # Clustering
# 
# Abaixo realizamos o clustering com o objetivo de organizar genes com base nas suas caraterísticas partilhadas. Iremos utilizar o clustering com distância de Manhatan dada a sua ampla utilização. 
# 
#%%

print('A preparar dados para clustering: escalar, estimar distância, clusterizar.')

# Assume 'rdf' is a DataFrame with a 'log2FoldChange' column
res_ord = rdf.sort_values(by='log2FoldChange').log2FoldChange.to_frame(name='log2FoldChange').head(20)

# Scaling the selected data
scaled = sklearn.preprocessing.scale(res_ord)

# Performing hierarchical/agglomerative clustering
Z = scipy.cluster.hierarchy.linkage(scaled, metric='cityblock')

# Create a figure
matplotlib.pyplot.figure(figsize=(10, 8))

# Plot the dendrogram
dendro = scipy.cluster.hierarchy.dendrogram(
    Z,
    labels=list(res_ord.index),
    leaf_rotation=90,
    leaf_font_size=10,
    color_threshold=0.7 * max(Z[:, 2]),  # Set a color threshold for better visual separation
    above_threshold_color='gray'  # Color for clusters above the threshold
)

# Add title and labels
matplotlib.pyplot.title("Dendrograma de Clusterização Hierárquica", fontsize=16)
matplotlib.pyplot.xlabel("Genes", fontsize=14)
matplotlib.pyplot.ylabel("Distância", fontsize=14)

# Customize the font and layout
matplotlib.pyplot.xticks(fontsize=10)
matplotlib.pyplot.yticks(fontsize=12)
matplotlib.pyplot.grid(True)

# Save the dendrogram plot to a PDF file
matplotlib.pyplot.savefig('./output/dendrograma.pdf', bbox_inches='tight')

# Show the plot
print('A exibir dendrograma. Dendrograma gravado em ./output/')
matplotlib.pyplot.show()

#%% md
# Acima é representado o dendrograma gerado para os resultados do clustering, refletindo a estrutura hierárquica dos clusters. É possível identificar dois grandes grupos. O grupo à esquerda parece consistir de genes tipicamente relacionados com fatores associados ao tecido tumoral, como codificação de proteínas transmembranas, patogenicidade, etc. O grupo à direita aparenta estar relacionado com fatores de agressividade do tumor como tem vindo a ser discutido ao longo desta análise. Estes genes poderão constituir biomarcadores ou alvos terapêuticos interessantes, estando já identificados na literatura, corroborando os resultados desta análise. 
# 
# - [Identification of potential diagnostic and prognostic biomarkers for LUAD based on TCGA and GEO databases](https://portlandpress.com/bioscirep/article/41/6/BSR20204370/228708/Identification-of-potential-diagnostic-and)
#%% md
# # Conclusões
# 
# Os resultados obtidos ao longo desta análise são diferentes dos obtidos na análise feita em R. Isto pode dever-se, em parte, ao facto de ter sido utilizada uma abordagem com DESEQ2, sendo feito um tratamento estatístico diferente dos dados, assim como a comparação entre categorias de dados ser diferente - sendo comparado PP vs TRU e PI vs TRU em vez de comparar categoria a categoria. Isto implica que as restantes análises sejam feitas partindo de dados pré-tratados de forma diferente. Em todo o caso os resultados são congruentes, sendo possível identificar diferenças entre as diferentes condições conforme seria de esperar e permitindo ainda identificar mais alvos e biomarcadores potenciais, pelo que acrescenta riqueza e variabilidade aos resultados deste trabalho e abre caminho para novas investigações de um conjunto alargado de genes e conjuntos de genes. Importa também referir que para a análise de enriquecimento foi utilizada também uma biblioteca diferente, sendo adotado uma biblioteca com gene sets focados em oncogenes. Concluímos assim que existem efetivamente diferenças entre os diferentes sub-tipos de adenocarcinoma do pulmão, sendo estes resultados corroborados pela literatura, como demonstramos ao longo do documento.
# 