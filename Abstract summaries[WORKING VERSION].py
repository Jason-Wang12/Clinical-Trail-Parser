from Bio import Entrez
from Bio import Medline
from tqdm import tqdm
import pandas as pd
pd.set_option('display.max_colwidth', -1)
import numpy as np
#Modules for bigram analysis
import nltk
from nltk import bigrams, trigrams
from nltk.corpus import stopwords
from nltk.tokenize import sent_tokenize, word_tokenize 
import collections
import re
#Graphing modules
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns

# Change this email to your email address
Entrez.email = "jwang@fortressbiotech.com"
'''
disease_df = pd.read_excel('file:///C:/Users/jwang/Desktop/Finance/Surrogate Endpoints FDA 12_16-19.xlsx', sheet_name = ' Adult SE table')#Only 
disease_df= disease_df[disease_df['Type of approval appropriate for'].str.contains('Accelerated')] #Filters for only diseases that qualify accelerated approval pathways
disease_df[' Patient Population '] = disease_df[' Patient Population '].apply(lambda x: re.findall(r'(?<=with )(?s)(.*$)', str(x))) #Remove the word all the words before the disease name
disease_df = disease_df[disease_df.astype(str)[' Patient Population '] != '[]'] #Removes blank cells with only brackets
disease_df.iloc[1]
disease_df.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/test4.xlsx')
disease_list = list(disease_df[' Patient Population '])
disease_list
'''
#############################################################################
disease_list = [input('Please enter a search term: ')]

#search and return total number of publications 
def search(x):
    Entrez.email=Entrez.email
    for x in disease_list:
        keyword = x
        handle = Entrez.esearch(db ='pubmed',
                                retmax=10,
                                retmode ='text',
                                term = keyword)
        results= Entrez.read(handle)
        print('Total number of publications that contain the term {}: {}'.format(keyword, results['Count']))    
    return results

def batch(x):

    keyword = disease_list
    handle = Entrez.esearch(db = 'pubmed',
                            retmax = results['Count'],
                            term = keyword,
                            mindate = 2010,
                            maxdate = 2020
                            )
    result = Entrez.read(handle)
    ids = result['IdList']
    batch_size = 100
    batches = [ids[x: x + 100] for x in range(0, len(ids), batch_size)]
    return batches

def fetch(x):
    record_list = []
    for batch in tqdm(batches):
        handle = Entrez.efetch(db="pubmed", 
                               id=batch, 
                               rettype="medline", 
                               retmode="text")
        records = Medline.parse(handle)
        record_list.extend(list(records))
        print('Complete.')
    return record_list
    
if __name__ == '__main__':
    results = search(disease_list)
    batches= batch(results)
    record_list = fetch(batches)

##################################################################################################################################################

abstract_table = pd.DataFrame(record_list, columns = ['AB']).dropna()


def processing(abstracts):
    #Lower case all words
    abstract_table['AB'] = abstract_table['AB'].str.lower()
    #remove numbers
    abstract_table['AB'] = abstract_table['AB'].replace(r'[0-9]+', '', regex = True)
    #remove special characters
    abstract_table['AB'] = abstract_table['AB'].replace(r'[!"#$%&()*+,-./:;<=>?@[\]^_`{|}~]', '', regex=True)    
    return abstract_table

if __name__ == '__main__':
    table_abstracts = [processing(record_list)]
    #print(table_abstracts[:2])

def stop_words(abstracts):
    #Remove Stop words
    stop_words_list = stopwords.words('english')
    abstract_table['AB'] = abstract_table['AB'].str.lower()
    abstract_table['AB'] = abstract_table['AB'].apply(lambda x: ' '.join([word for word in x.split() 
    if word not in (stop_words_list)]))
    return abstract_table
  
def collection_words(abstracts):
    #Remove collection words
    collection_words = ['natural', 'killer', 'nk', 'cells', 'mononuclear', 'blood', 'adaptive', 'innate', 'immune', 'acute', 'myeloid', 'leukemia', 'aml']
    abstract_table['AB'] = abstract_table['AB'].str.lower()
    abstract_table['AB'] = abstract_table['AB'].apply(lambda x: ' '.join([word for word in x.split() if word not in (collection_words)]))
    return abstract_table

def clean_abstracts(abstracts):
    table_tweet = stop_words(abstract_table)
    table_tweet = collection_words(abstract_table)
    return table_tweet

if __name__ == '__main__':
    table_abstracts = clean_abstracts(table_abstracts)

##################################################################################################################################################
#Combine the abstracts together

abx_nsw= [abx.split() for abx in table_abstracts['AB']]
term_bigrams = [list(bigrams(abx)) for abx in abx_nsw]
bigrams_flat = [word for words in term_bigrams for word in words]
bigrams_count = collections.Counter(bigrams_flat)
bigrams_df= pd.DataFrame(bigrams_count.most_common(50), columns = ['bigrams', 'count'])

#Now graph the bigrams out
d = bigrams_df.set_index('bigrams').T.to_dict('records')
G = nx.Graph()
for k,v in d[0].items():
    G.add_edge(k[0], k[1], weight = (v *10))
fig, ax = plt.subplots(figsize=(20, 14))
pos = nx.spring_layout(G, k=0.75,iterations=20)

# Plot networks
nx.draw_networkx(G, pos,
                 font_size=16,
                 width=3,
                 edge_color='lightgray',
                 node_color='lightblue',
                 with_labels = False,
                 ax=ax)

                  
# Create offset labels
for key, value in pos.items():
    x, y = value[0]+.135, value[1]+.045
    ax.text(x, y,
            s=key,
            bbox=dict(facecolor='red', alpha=0.25),
            horizontalalignment='center', fontsize=12)
    
plt.show()

sns.set_style("white")
plt.figure(figsize=(20,12), dpi=150)
plt.subplot(2, 2, 1)
sns.set(font_scale = 0.5)
sns.barplot(x="count", y="bigrams", data=bigrams_df, palette="RdBu_r")
plt.title('Bigrams')

#Combine the abstracts together

abx_nsw= [abx.split() for abx in table_abstracts['AB']]
term_trigrams = [list(trigrams(abx)) for abx in abx_nsw]
trigrams_flat = [word for words in term_trigrams for word in words]
trigrams_count = collections.Counter(trigrams_flat)
trigrams_df= pd.DataFrame(trigrams_count.most_common(30), columns = ['trigrams', 'count'])

#Now graph the trigrams out
d = trigrams_df.set_index('trigrams').T.to_dict('records')
G = nx.Graph()
for k,v in d[0].items():
    G.add_edge(k[0], k[1], weight = (v *10))
fig, ax = plt.subplots(figsize=(20, 14))
pos = nx.spring_layout(G, k=0.75,iterations=20)

# Plot networks
nx.draw_networkx(G, pos,
                 font_size=16,
                 width=3,
                 edge_color='lightgray',
                 node_color='lightblue',
                 with_labels = False,
                 ax=ax)

                  
# Create offset labels
for key, value in pos.items():
    x, y = value[0]+.135, value[1]+.045
    ax.text(x, y,
            s=key,
            bbox=dict(facecolor='red', alpha=0.25),
            horizontalalignment='center', fontsize=12)
    
plt.show()

sns.set_style("white")
plt.figure(figsize=(12, 8), dpi=150)
plt.subplot(2, 2, 1)
sns.set(font_scale = 0.5)
sns.barplot(x="count", y="trigrams", data=trigrams_df, palette="RdBu_r")
plt.title('trigrams')

'''
filter for only key terms you care about
'''
title_table= pd.DataFrame(record_list, columns = ['TI', 'AB', 'AD']).dropna() #Pull out title, abstract, contact info
clean_table = title_table
clean_table['AB'] = clean_table['AB'].str.lower()
clean_table['TI'] = clean_table['TI'].str.lower()
search_terms = ['phase', 'trial', 'patients']
#Since it's a dataframe, we need to convert a list of strings into a series that can be applied to the dataframe
clean_table= clean_table[clean_table['TI'].str.contains('|'.join(search_terms))] #Make sure it's a clincial study
#clean_table = clean_table[~clean_table.TI.str.contains('Newly, newly')]  #Turn this on if you want to filter out studies looking at front line therapy patients
obv_study_terms = ['demographic','association', 'associated', 'retrospective', 'prognostic', 'prediction', 'detection']  
clean_table = clean_table[~clean_table['AB'].str.contains('|'.join(obv_study_terms))]
clean_table['AB'] = clean_table['AB'].apply(lambda x: re.findall(r'results:.*|findings:.*', str(x))) # remove all the text ahead of results
clean_table = clean_table[clean_table.astype(str)['AB'] != '[]'] #Remove all the blank cels that result 
#clean_table = clean_table[clean_table['AD'].str.contains('@')] #Keep only results that have an email address
clean_table['AD'] = clean_table['AD'].apply(lambda x: re.findall(r'[\w\.-]+@[\w\.-]+', str(x))) #Remove all institution/affiliation info
#Save as a CSV so that if you need to pull up the data again, you don't have to pull from pubmed during working hrs

clean_table.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/{}.xlsx'.format(disease_list))
################################################################
'''
Filter out only publications which have results
'''

#Extract results number after results text
#turn abbreviation into full word for next line
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np


#clean_table = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/nsclc treatment.xlsx')
df = clean_table

#convert column to a string
df['AB'] = df['AB'].astype(str)
df['TI'] = df['TI'].astype(str)
df['drug names'] = df['AB'].str.strip('[]') #remove the brackets, so it's recognized as a string
df['drug names'] = df['drug names'].str.replace(r'NCT\d+$', '')
df['drug names'] = df['drug names'].str.extract(r'([A-Za-z]{2,6}\d{2,8}\b|\w+mab\b | \w+nib\b |\w+nib\b|\w+sertib\b|\w+stat\b|\w+staurin\b|\w+bine\b|\w+cef\b|\w+asone\w+bital\b|\w+cycline\b|\w+zole\b|\w+platin\b)')
#df['drug names'] = df['drug names'].str.extract(r'([^nct][a-z]{1,4}\d{2,7})|\w+stat\b|')
df.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/test.xlsx')

'''
Parse out late stage drug resuls
'''

df['mOS'] = df['AB'].str.replace('os', 'overall survival')
#df['mOS'] = df['AB'].str.extract('(?<=overall)*[ a-zA-Z]* (\d+.\d+ month)')
#df['mOS'] = df['AB'].str.extract('(?<=overall survival).*? (\d+.\d+ month)')
#df['mOS'] = df['AB'].str.extract('(?<=overall survival).*? (\d+.\d+ month)')
df['mOS'] = df['mOS'].str.extract('(?<=overall)*[^/]* (\d+.\d+ month)')
df['mOS']= df['mOS'].str.strip('month')


df['PFS'] = df['AB'].str.replace('pfs', 'progression free survival')
#df['PFS'] = df['AB'].str.extract('(?<=progression free survival)*[^/]* (\d+.\d+ month)')
df['PFS'] = df['PFS'].str.extract('(?<=progression free).*? (\d+.\d+ month)')
df['PFS']= df['PFS'].str.strip('month')

'''
Parse out early stage drug results
'''
df['ORR'] = df['AB'].str.replace('orr | overall response', 'objective response')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.+?%)')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.|$)%')
df['ORR'] = df['ORR'].str.extract('(?<=objective response).*?(\d.|$%)')
df['ORR']= df['ORR'].str.strip('%')

df['PR'] = df['AB'].str.replace('partial response', 'pr')
df['PR']= df[~df['PR'].str.contains('pruritus', na=False)]
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.+?%)')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.|$)%')
df['PR'] = df['PR'].str.extract('(?<=pr).*?(\d.|$%)')
df['PR']= df['PR'].str.strip('%')

df['CR'] = df['AB'].str.replace('complete response', 'cr')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.+?%)')
#df['ORR'] = df['AB'].str.extract('(?<=objective response).+?(\d.|$)%')
df['CR'] = df['CR'].str.extract('(?<=cr).*?(\d.|$%)')
df['CR']= df['CR'].str.strip('%')

#Results table for survival
survival = pd.DataFrame()
survival = survival.append(df['drug names'])
survival =survival.append(df['mOS'])
survival = survival.T.dropna()
survival = survival[~survival.mOS.str.contains('-')]
survival = survival.T
print(survival)

progression = pd.DataFrame()
progression = progression.append(df['drug names'])
progression =progression.append(df['PFS'])
progression = progression.T.dropna()
progression = progression[~progression.PFS.str.contains('-')]
progression = progression.T
print(progression)

#Results table for tumor shrinkage
shrink = pd.DataFrame()
shrink = shrink.append(df['drug names'])
shrink = shrink.append(df['ORR'])
shrink = shrink.T.dropna()
shrink = shrink[shrink['ORR'].apply(lambda x: str(x).isdigit())]
shrink = shrink.T
print(shrink)

#Results for partial and complete response
Response = pd.DataFrame()
Response = Response.append(df['drug names'])
Response = Response.append(df['PR'])
Response = Response.append(df['CR'])
Response = Response.fillna(0)
Response = Response.T
Response = Response[Response['PR'].apply(lambda x: str(x).isdigit())]
Response = Response[Response['CR'].apply(lambda x: str(x).isdigit())]
print(Response)


survival.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/survival.xlsx')
progression.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/progression.xlsx')
shrink.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/shrink.xlsx')  
Response.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/response.xlsx')

late_stage = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/survival.xlsx', sheet_name = 'Sheet1', header = 1, index_col =0)
progression = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/progression.xlsx', sheet_name = 'Sheet1', header = 1, index_col = 0)
early_stage = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/shrink.xlsx', sheet_name = 'Sheet1', header = 1, index_col =0)
response = pd.read_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/response.xlsx', sheet_name = 'Sheet1', header = 0, index_col = 0)

'''
Overall Survival
'''
OS = late_stage.loc['mOS']
#OS = OS.dropna()
OS = OS.sort_values(ascending = True)
OS[100:].plot(kind = 'barh', figsize = (8,6), fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'mOS', legend = False)
#OS.iloc[np.r_[1:10, 15:20].plot(kind = 'barh', figsize = (8,6), fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'mOS', legend = False)
plt.show()
'''
PFS
'''
PFS = progression.loc['PFS']
#PFS = PFS.dropna()
PFS = PFS.sort_values(ascending = True)
#PFS.iloc[np.r_[49:59,85:89]].plot(kind = 'barh', figsize = (15, 10), fontsize = 12, cmap = plt.get_cmap('plasma'), title = 'PFS', legend = False)
PFS.plot(kind = 'barh', figsize = (15, 10), fontsize = 12, cmap = plt.get_cmap('plasma'), title = 'PFS', legend = False)
plt.show()
'''
Objective Response Rate
'''
ORR = early_stage.loc['ORR']
ORR = ORR.dropna()
ORR = ORR.sort_values(ascending = False)
#ORR.iloc[np.r_[1:5, 10:15]].plot(kind = 'bar', figsize = (8, 6),fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'Overall Response Rate', legend = False)
ORR[:15].plot(kind = 'bar', figsize = (8, 6),fontsize = 12, cmap = plt.get_cmap('viridis'), title = 'Overall Response Rate', legend = False)
plt.show()

#######################
'''
CR + PR
'''
CR = response.loc['CR'].fillna(0)
PR = response.loc['PR'].fillna(0)
ORR = pd.concat([CR, PR], axis = 1)
ORR['sum'] = ORR['CR']+ORR['PR']
ORR['sum'] = ORR[ORR['sum'] !='-']
ORR = ORR.sort_values(by = 'sum', ascending = False)
plt.rc('xtick', labelsize=20) 
ORR['sum'].plot(kind = 'bar', stacked = True, figsize = (10,8))
plt.show()
