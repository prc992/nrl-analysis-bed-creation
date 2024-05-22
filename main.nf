 #! /usr/bin/env nextflow
nextflow.enable.dsl=2

process split_tss_samples{
  tag "Sample - $sampleId" 
  //maxForks 3

  //Docker Image
  container = 'prc992/pyranges:v1.2'
  label 'default_mem'

  publishDir "$bedFolderOut", mode : 'copy', pattern : '*.bed'
  publishDir "$csvFolderOut", mode : 'copy', pattern : '*.csv'


  input:
  tuple val(sampleId), path(csvFile),path(bedFileIn),val(bedFolderOut),val(csvFolderOut),val(WINDOW_TSS_DOWNSTREAM),val(WINDOW_TSS_UPSTREAM)

  output:
  path("*.bed")
  path("*.csv")

  script:
  """
  python $params.python_prog $csvFile $WINDOW_TSS_DOWNSTREAM $WINDOW_TSS_UPSTREAM $params.MAX_SIZE_FRAGMENT $params.cpu
  """
}

workflow {

    chSampleInfo = Channel.fromPath(params.samples) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId,row.csvFile,row.bedFileIn,row.bedFolderOut,row.csvFolderOut\
                      ,row.WINDOW_TSS_DOWNSTREAM,row.WINDOW_TSS_UPSTREAM) }

    split_tss_samples(chSampleInfo)
}

