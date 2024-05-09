import sys
import pandas as pd
import pyranges as pr

def split_tss_samples(csvFile,WINDOW_TSS_DOWNSTREAM,WINDOW_TSS_UPSTREAM,MAX_SIZE_FRAGMENT):
   
    df = pd.read_csv(csvFile,delimiter= ',',na_filter=False)

    for ind in df.index:
        gene_name = df['Name'][ind]
        bed_file_full = df['ssd-file'][ind]
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

    csvFileOut = "4-Splited-Samples-Selected-TSS-" + str(individual_id) + '.csv'
    df.to_csv(csvFileOut)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python program_name.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    WINDOW_TSS_DOWNSTREAM = sys.argv[2]
    WINDOW_TSS_UPSTREAM = sys.argv[3]
    MAX_SIZE_FRAGMENT = sys.argv[4]

    split_tss_samples(filename,WINDOW_TSS_DOWNSTREAM,WINDOW_TSS_UPSTREAM,MAX_SIZE_FRAGMENT)
