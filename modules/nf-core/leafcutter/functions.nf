
process BAM_JUNC{
  label 'process_medium'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/huvec/splicing/leafcutter_25_09_2023.img"
  } else {
    container "wtsihgi/nf_cellbender_container:3cc9983"
  }


  input:
    tuple(val(name),path(bam),path(bai))

  output:
    path("*.junc", emit: junction_file)

   script:
    """
    regtools junctions extract -s XS -a 8 -m 50 -M 500000 ${bam} -o ${bam}.junc 
    
    """

}


process LEAFCUTTER_CLUSTERING{

  label 'process_medium'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/huvec/splicing/leafcutter_25_09_2023.img"
  } else {
    container "wtsihgi/nf_cellbender_container:3cc9983"
  }

  input:
    path(junction_files)

  output:
    path('combined_perind.counts.gz.qqnorm_phen_aggregated', emit: r1)
    path('combined_perind.counts.gz.qqnorm_chrALL_aggregated', emit: r2)
    path('clustered_junction_file_refined', emit: r3)
    path('clustered_junction_file_perind.counts.gz.PCs', emit: r4)
    path('clustered_junction_file_perind.counts.gz.ave', emit: r5)
    path('clustered_junction_file_perind.counts.gz', emit: r6)

  script:
    """
    for junc_file in ${junction_files}; do
        echo \$junc_file >> juncfiles.txt
    done
    
    leafcutter_cluster_regtools -j juncfiles.txt -m 50 -o clustered_junction_file -l 500000    
    prepare_phenotype_table clustered_junction_file_perind.counts.gz -p 10


    agg_chr_file=combined_perind.counts.gz.qqnorm_chrALL_aggregated
    head -1 clustered_junction_file_perind.counts.gz.qqnorm_1 > \$agg_chr_file
    for f in *_perind.counts.gz.qqnorm*; do
        echo \$f
      if [ \$f != \$agg_chr_file ]; then tail -n+2 \$f >> \$agg_chr_file; fi;
    done;


    agg_chr_file=combined_perind.counts.gz.qqnorm_phen_aggregated
    head -1 clustered_junction_file_perind.counts.gz.phen_1 > \$agg_chr_file
    for f in *_perind.counts.gz.phen*; do
        echo \$f
      if [ \$f != \$agg_chr_file ]; then tail -n+2 \$f >> \$agg_chr_file; fi;
    done;

    """
}
