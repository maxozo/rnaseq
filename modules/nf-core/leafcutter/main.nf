
include {
    BAM_JUNC;
    LEAFCUTTER_CLUSTERING;
} from "./functions.nf"

workflow LEAFCUTTER{
    take:
        name_bam_bai
    main:
        BAM_JUNC(name_bam_bai)
        LEAFCUTTER_CLUSTERING(BAM_JUNC.out.junction_file.collect())
}
