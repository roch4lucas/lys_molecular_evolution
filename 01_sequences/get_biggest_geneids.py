#!/usr/bin/env python3
"""
Description: Python script to get the longest protein sequence for a list of 
NCBI GeneIDs, excluding pseudogenes, using NCBI Entrez. Script verbose is in Portuguese

- Input: TSV file with 2 columns: GeneIDs / Organisms
- Output: FASTA file

Don't forget to add a valid NCBI Entrez email (Entrez.mail)

Author: Lucas Rocha
"""

import os
import time
import concurrent.futures
from Bio import Entrez, SeqIO
import pandas as pd
import threading

# Email config
Entrez.email = "your@email.com"  # add your own

print_lock = threading.Lock()

def ler_arquivo_geneids(arquivo):
    """Lê o arquivo TSV contendo GeneIDs e organismos."""
    try:
        df = pd.read_csv(arquivo, sep='\t', header=None)
        # Assumindo que o arquivo tem duas colunas: GeneID e organismo
        if len(df.columns) >= 2:
            df.columns = ['gene_id', 'organismo'] + [f'coluna_{i}' for i in range(2, len(df.columns))]
        else:
            df.columns = ['gene_id']
            df['organismo'] = None
        return df
    except Exception as e:
        with print_lock:
            print(f"Error reading GeneID file: {e}")
        return None

def eh_pseudogene(record):
    """
    Verifica se um registro de proteína está associado a um pseudogene.
    
    Args:
        record: Um objeto SeqRecord do Biopython com informações da proteína
    
    Returns:
        bool: True se for um pseudogene, False caso contrário
    """
    # Palavras-chave que podem indicar um pseudogene
    keywords = ['pseudogene', 'pseudo', 'nonfunctional', 'non-functional']
    
    # Verificar no título/descrição
    if any(keyword in record.description.lower() for keyword in keywords):
        return True
        
    # Verificar nas anotações
    for feature in record.features:
        # Verificar qualificadores como note, product, gene_desc
        for qualifier in ['note', 'product', 'gene_desc']:
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    if any(keyword in value.lower() for keyword in keywords):
                        return True
                        
        # Verificar no qualificador pseudogene específico
        if 'pseudogene' in feature.qualifiers:
            return True
            
    # Verificar em outras anotações
    for annotation_key, annotation_value in record.annotations.items():
        if isinstance(annotation_value, str):
            if any(keyword in annotation_value.lower() for keyword in keywords):
                return True
        elif isinstance(annotation_value, list):
            for item in annotation_value:
                if isinstance(item, str) and any(keyword in item.lower() for keyword in keywords):
                    return True
    
    return False

def processar_gene(gene_info):
    """
    Processa um único gene para encontrar a maior proteína não-pseudogene.
    
    Args:
        gene_info: Tupla contendo (gene_id, organismo)
    
    Returns:
        tuple: (gene_id, proteína, protein_id) ou (gene_id, None, None) se não encontrado
    """
    gene_id, organismo = gene_info
    
    with print_lock:
        print(f"Processando gene {gene_id}...", end=" ", flush=True)
    
    try:
        # Construir a consulta
        query = f"{gene_id}[Gene ID]"
        if organismo:
            query += f" AND {organismo}[Organism]"
        query += " NOT pseudogene[All Fields]"  # Tentar excluir pseudogenes na consulta
        
        # Primeiro buscar os IDs de proteínas associados ao gene
        search_handle = Entrez.esearch(db="protein", term=query, retmax=100)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        protein_ids = search_results.get("IdList", [])
        
        if not protein_ids:
            with print_lock:
                print(f"Nenhuma proteína encontrada")
            return gene_id, None, None
        
        # Buscar informações para todas as proteínas encontradas
        fetch_handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="gb", retmode="text")
        records = list(SeqIO.parse(fetch_handle, "genbank"))
        fetch_handle.close()
        
        # Filtrar pseudogenes
        filtered_records = []
        for record in records:
            if not eh_pseudogene(record):
                filtered_records.append(record)
            else:
                with print_lock:
                    print(f"  Pseudogene ignorado: {record.id}")
        
        # Encontrar a sequência mais longa entre os registros filtrados
        if filtered_records:
            maior_proteina = max(filtered_records, key=lambda rec: len(rec.seq))
            # Adicionar o GeneID à descrição da sequência
            maior_proteina.description = f"GeneID:{gene_id} | {maior_proteina.description}"
            with print_lock:
                print(f"OK - Proteína: {maior_proteina.id} ({len(maior_proteina.seq)} aa)")
            return gene_id, maior_proteina, maior_proteina.id
        else:
            with print_lock:
                print(f"Nenhuma proteína válida (todos pseudogenes)")
            return gene_id, None, None
            
    except Exception as e:
        with print_lock:
            print(f"Erro: {e}")
        return gene_id, None, None

def main():
    arquivo_entrada = "geneids.tsv"
    arquivo_saida = "proteinas_maiores.fasta"
    
    # Número de workers para processamento paralelo
    # Ajuste este valor com base na sua conexão de internet e CPU
    # Um valor entre 4-10 normalmente é bom para equilibrar desempenho e 
    # respeitar os limites de requisições do NCBI
    num_workers = 5
    
    start_time = time.time()
    
    # Ler o arquivo de entrada
    df = ler_arquivo_geneids(arquivo_entrada)
    if df is None:
        return
    
    print(f"Processando {len(df)} genes em paralelo com {num_workers} workers...")
    
    # Criar lista de tuplas (gene_id, organismo) para processamento
    genes_para_processar = [(str(row['gene_id']), row.get('organismo')) for _, row in df.iterrows()]
    
    # Lista para armazenar os resultados
    resultados = {}
    
    # Executar o processamento paralelo
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submeter todas as tarefas
        future_to_gene = {
            executor.submit(processar_gene, gene_info): gene_info[0] 
            for gene_info in genes_para_processar
        }
        
        # Processar os resultados conforme eles são concluídos
        completed = 0
        for future in concurrent.futures.as_completed(future_to_gene):
            gene_id = future_to_gene[future]
            try:
                gene_id, proteina, protein_id = future.result()
                resultados[gene_id] = proteina
                
                # Atualizar contador de progresso
                completed += 1
                if completed % 5 == 0 or completed == len(genes_para_processar):
                    print(f"Progresso: {completed}/{len(genes_para_processar)} genes processados")
                
            except Exception as e:
                print(f"Erro ao processar gene {gene_id}: {e}")
    
    # Filtrar apenas proteínas encontradas
    proteinas_encontradas = [proteina for proteina in resultados.values() if proteina is not None]
    
    # Salvar os resultados no arquivo FASTA
    if proteinas_encontradas:
        with open(arquivo_saida, "w") as output_handle:
            SeqIO.write(proteinas_encontradas, output_handle, "fasta")
        
        end_time = time.time()
        duracao = end_time - start_time
        
        print(f"\nProcessamento concluído em {duracao:.2f} segundos.")
        print(f"{len(proteinas_encontradas)} proteínas foram salvas em '{arquivo_saida}'.")
        print(f"Genes não encontrados ou que possuíam apenas pseudogenes: {len(df) - len(proteinas_encontradas)}")
    else:
        print("\nNenhuma proteína foi encontrada.")

if __name__ == "__main__":
    main()
