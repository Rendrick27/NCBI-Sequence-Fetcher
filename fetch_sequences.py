# Linux shebang line
#!/usr/bin/env python3

# Importa a biblioteca Bio.Entrez para acessar bancos de dados do NCBI
from Bio import Entrez

# Importa o módulo sys para lidar com argumentos de linha de comando
import sys

# Define a função fetch_sequences que busca sequências em um banco de dados
def fetch_sequences(database, term):

    # Definir um email para a API do Entrez
    Entrez.email = "fakeemail@exemplo.com"

    # Definir o nome da ferramenta que está sendo usada para acessar a API do Entrez
    Entrez.tool = "fetch_sequences.py"

    # Usa o método esearch para pesquisar o banco de dados com o termo de busca e salva os resultados em um histórico
    search_results = Entrez.esearch(db=database, term=term, usehistory="y")
    search_results = Entrez.read(search_results)

    # Salva o WebEnv e QueryKey para uso em chamadas API posteriores
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    # Usa o método efetch para buscar as sequências em lotes e escrevê-las no console
    batch_size = 10 # Número de registros para buscar de uma vez
    for start in range(0, int(search_results["Count"]), batch_size):
        end = min(int(search_results["Count"]), start+batch_size) # Calcula o final do lote atual
        fetch_results = Entrez.efetch(db=database, rettype="fasta", retmode="text",
                                      retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        sys.stdout.write(fetch_results.read()) # Escreve os resultados no console
        fetch_results.close() # Fecha a conexão para liberar recursos

# Verifica se o script está sendo executado como o programa principal
if __name__ == '__main__':

    # Verifica se o número de argumentos está correto
    if len(sys.argv) != 3:
        # Se o número de argumentos estiver incorreto, imprime informações de uso e sai com um código de erro
        print("Uso: python fetch_sequences.py [nome do banco de dados] [termo de busca]")
        sys.exit(1)

    # Obtém o nome do banco de dados e o termo de busca a partir dos argumentos da linha de comando
    database = sys.argv[1]
    term = sys.argv[2]

    # Chama a função fetch_sequences com os argumentos fornecidos pelo usuário
    fetch_sequences(database, term)