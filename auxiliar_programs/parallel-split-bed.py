import pandas as pd
import pyranges as pr
import multiprocessing as mp
import sys

# Dividir o arquivo CSV em n_chunks partes
def divide_csv_file(file_path, n_chunks):
    df = pd.read_csv(file_path)
    lengthDF = len(df)
    chunk_size = lengthDF // n_chunks
    dividedDF = [df.iloc[i*chunk_size:(i+1)*chunk_size] for i in range(n_chunks)]
    
    # Incluir quaisquer linhas restantes no último chunk
    if lengthDF % n_chunks != 0:
        dividedDF[-1] = pd.concat([dividedDF[-1], df.iloc[n_chunks*chunk_size:]])
    
    return dividedDF

# Função para somar os valores da coluna 'Start' em um DataFrame
def sum_df_start(args):

    df, WINDOW_TSS_DOWNSTREAM,WINDOW_TSS_UPSTREAM, MAX_SIZE_FRAGMENT = args

    for ind in df.index:
        gene_name = df['Name'][ind]
        bed_file_full = df['ssd-file'][ind]
        bed_file = bed_file_full
        bed_file = bed_file_full.rsplit('/', 1)[-1]
        antibody = df['antibody'][ind]
        tss_chromosome = df['Chromosome'][ind]
        tss = int(df['TSS'][ind])
        individual_id = df['sample_id'][ind]
        strand = df['Strand'][ind]

        WINDOW_TSS_DOWNSTREAM = int(WINDOW_TSS_DOWNSTREAM)
        WINDOW_TSS_UPSTREAM = int(WINDOW_TSS_UPSTREAM)
        MAX_SIZE_FRAGMENT = int(MAX_SIZE_FRAGMENT)

        if strand == '-':
            window_tss_start = tss - WINDOW_TSS_DOWNSTREAM
            window_tss_end = tss + WINDOW_TSS_UPSTREAM
        else:
            window_tss_start = tss - WINDOW_TSS_UPSTREAM
            window_tss_end = tss + WINDOW_TSS_DOWNSTREAM

        bed_file_whole = pr.read_bed(bed_file)
        df_bed_file_filtered = bed_file_whole[tss_chromosome,window_tss_start:window_tss_end].as_df()

        if len(df_bed_file_filtered)>0:
            df_bed_file_filtered['fragment_size'] = df_bed_file_filtered['End'] - df_bed_file_filtered['Start']
            df_bed_file_filtered = df_bed_file_filtered[df_bed_file_filtered['fragment_size']<MAX_SIZE_FRAGMENT]
            df_bed_file_filtered = df_bed_file_filtered.drop(columns=['fragment_size'])

            qtd_frags = len(df_bed_file_filtered)

            if qtd_frags > 0:

                prBED = pr.PyRanges(df_bed_file_filtered)
                bedFileOut = str(individual_id) + "-" + antibody + '-' + gene_name + '.bed'
                prBED.to_bed(bedFileOut)

                df.loc[ind,'SSD-TSS-ANTIBODY-BEDFILE'] = bedFileOut
                df.loc[ind,'window_tss_start'] = window_tss_start
                df.loc[ind,'window_tss_end'] = window_tss_end

    return df

# Função para realizar o processamento em paralelo
def process_in_parallel(file_path, n_chunks,WINDOW_TSS_DOWNSTREAM,WINDOW_TSS_UPSTREAM,MAX_SIZE_FRAGMENT):
    # Dividir o CSV em partes
    dividedDF = divide_csv_file(file_path, n_chunks)
    
    # Criar um pool de processos
    pool = mp.Pool(processes=n_chunks)

    args = [(chunk,WINDOW_TSS_DOWNSTREAM,WINDOW_TSS_UPSTREAM, MAX_SIZE_FRAGMENT) for chunk in dividedDF]
    
    # Mapear a função para os chunks
    dfs_results = pool.map(sum_df_start, args)
    
    # Fechar o pool e esperar pelos processos
    pool.close()
    pool.join()
    
    # Somar os resultados de todos os chunks
    prefix = '4-Samples'
    strFile = file_path.replace('3-Samples',prefix)

    combined_df = pd.concat(dfs_results, ignore_index=True)
    strFile1 = strFile.rsplit('/', 1)[-1]
    combined_df.to_csv(strFile1)


if __name__ == "__main__":
    # Nome do arquivo CSV e número de chunks

    if len(sys.argv) != 6:
        print("Usage: python program_name.py <filename>")
        sys.exit(1)

    file_path = sys.argv[1]
    WINDOW_TSS_DOWNSTREAM = sys.argv[2]
    WINDOW_TSS_UPSTREAM = sys.argv[3]
    MAX_SIZE_FRAGMENT = sys.argv[4]
    n_cpus = int(sys.argv[5])
    
    # Chamar a função de processamento paralelo 
    process_in_parallel(file_path, n_cpus,WINDOW_TSS_DOWNSTREAM,WINDOW_TSS_UPSTREAM,MAX_SIZE_FRAGMENT)
