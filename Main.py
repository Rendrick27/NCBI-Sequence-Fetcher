#!/usr/bin/env python3
# Importar o módulo Bio.Entrez
from Bio import Entrez

# Definir o email e a ferramenta para a API Entrez
Entrez.email = "rendrickcarreira.social@gmail.com"
Entrez.tool = "meu_script"

# Receber os argumentos da linha de comando
import sys
database = sys.argv[1] # escolher a base de dados
term = sys.argv[2] # passar o termo de pesquisa

# Fazer uma consulta usando a função esearch com o recurso history
handle = Entrez.esearch(db=database, term=term, usehistory="y")
record = Entrez.read(handle)
handle.close()

# Obter os identificadores das sequências encontradas
webenv = record["WebEnv"]
query_key = record["QueryKey"]
count = int(record["Count"])

# Fazer um pedido usando a função efetch para obter as sequências em formato FASTA
batch_size = 10 # definir o tamanho do lote
for start in range(0, count, batch_size):
    end = min(count, start+batch_size)
    print("A obter as sequências %i até %i" % (start+1, end))
    fetch_handle = Entrez.efetch(db=database,
                                 rettype="fasta",
                                 retmode="text",
                                 retstart=start,
                                 retmax=batch_size,
                                 webenv=webenv,
                                 query_key=query_key)
    data = fetch_handle.read()
    fetch_handle.close()
    # Escrever os dados para STDOUT
    sys.stdout.write(data)

##Link: https://stuntspt.gitlab.io/asb2023/classes/class_05/index.html#/
##Criado por
#David
#Helder
#Joana
#Rendrick