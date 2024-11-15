#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"

multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    case 'latest':
    gatk_image="gatk4500.sif";
    default:
    gatk_image="gatk4400.sif";
    break;
}


switch (params.server) {

    case 'lnx02':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"


        //modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"

        modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
    case 'rgi01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"

        //modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules/";
    break;
}


switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "/data/shared/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "/data/shared/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        break;


    case 'hg38':
        assembly="hg38"
        smncaller_assembly="38"
        grch_assembly="grch38"
        pcgr_assembly="grch38"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta = "/data/shared/genomes/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/GRCh38.primary.dict"
        genome_version="V1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version="V2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="V3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared//genomes/hg38/inhouse_DBs/hg38v3_primary/"
        }

        // Gene and transcript annotation files:
        refFlat="/data/shared/genomes/hg38/gene.annotations/refFlat.txt"

        gencode_gtf = "/data/shared/genomes/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "/data/shared/genomes/hg38/gene.annotations/gencode.v36.annotation.gff3"
     
        //Program  files:
        msisensor_list="/data/shared/genomes/hg38/program_DBs/msisensor/hg38v3_msisensor_scan.txt"
        
        accucopy_config="/data/shared/genomes/hg38/accucopy/accucopy.docker.nextflow.conf"
        cnvradar_anno="/data/shared/genomes/hg38/program_DBs/cnvradar/All_20180418.vcf.gz"
        cnvradar_anno_idx="/data/shared/genomes/hg38/program_DBs/cnvradar/All_20180418.vcf.gz.tbi"
        cnvradar_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed" 

        cnvradar_roisum_dir="/data/shared/genomes/hg38/program_DBs/cnvradar/inhouse_roi_summaries/"
       


        //Structural variants
        delly_exclude="/data/shared/genomes/hg38/program_DBs/delly/human.hg38.excl.tsv"
        
        smoove_exclude="/data/shared/genomes/hg38/interval.files/smoove/smoove.hg38.excluderegions.bed"
        smoove_gff="/data/shared/genomes/hg38/gene.annotations/GRCh38_latest_genomic.gff.gz"



        //Repeat Expansions:
        expansionhunter_catalog="/data/shared/genomes/hg38/program_DBs/expansionHunter/expansionHunter_hg38_stripy.variant_catalog.json"
        hipSTR_bed="/data/shared/genomes/hg38/interval.files/STRs/GRCh38.hipstr_reference.bed"

        // Somatic calling files:
        gatk_wgs_pon="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"


        pcgr_data_dir="/data/shared/genomes/hg38/program_DBs/PCGR/"
        pcgr_VEP="/data/shared/genomes/hg38/program_DBs/PCGRv2/VEP_112_GRCh38_merged/"
        pcgr_data_dir2="/data/shared/genomes/hg38/program_DBs/PCGRv2/20240621/"
        pcgr_data_dir3="/data/shared/genomes/hg38/program_DBs/PCGRv2/20240927/"
        hmftools_data_dir_v534="/data/shared/genomes/hg38/program_DBs/hmftools/v5_34/ref/38"
        hmftools_data_dir_v60="/data/shared/genomes/hg38/program_DBs/hmftools/v6_0/ref/38"

        // Program indexes:
        
        sequenza_cg50_wig="/data/shared/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"


        // Regions & variants:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        inhouse127_geneIntervals="/data/shared/genomes/hg38/interval.files/geneIntervals/241022_inhouse127genes.3col.SM.bed"

        //ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.bed"

        callable_regions="/data/shared/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions="/data/shared/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        dbsnp="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KGindels="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
        KGindels_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

        KGmills="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        KGmills_idx="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        KG_p1_High_snps="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

        hapmap="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
        omni="/data/shared/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
        //WES_ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"

        break;
}
/*
switch (params.panel) {
    case "WES_2":
        ROI="${WES_ROI}";
        panelID="WES"
    break;

    case "WES":
        ROI="${WES_ROI}";
        panelID="WES_subpanel"
    break;

    default: 
        ROI="${WES_ROI}";
        panelID="WES"
    break;
}
*/

dataStorage="/lnx01_data3/storage/";
variantStorage="${dataStorage}/variantStorage/${params.genome}/"
outputDir="${params.outdir}/"



Channel
    .fromPath(params.intervals_list)
    .map { it -> tuple(it.baseName,it)}
    .set { HTC_interval_list }


// SYMLINK PROCESSES

process inputFiles_symlinks_fq{
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/fq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    input:
    tuple val(caseID), val(sampleID), path(r1),path(r2), val(type)// from read_input2
    
    output:
    tuple path(r1),path(r2)
    script:
    """
    """
}



process inputFiles_symlinks_cram{
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/cram_TN_symlinks/", mode: 'link', pattern: '*.{ba,cr}*'
    publishDir "${caseID}/${outputDir}/variantcalls/Alignment_symlinks/", mode: 'link', pattern: "*.{ba,cr}*"

    input:

    tuple val(caseID), val(sampleID), path(cram), path(crai),val(type)    
    output:
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram")
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai")
    script:
    """
    mv ${cram} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram
    mv ${crai} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai
    """
}


///////////////////////// PREPROCESS MODULES ////////////////////////////////////

process tb_fastq_to_ubam {
    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 20
    maxForks 10

    input:
    tuple val(caseID),val(sampleID), path(r1), path(r2), val(type)

    output:
    tuple val(caseID),val(sampleID), path("${sampleID}.unmapped.from.fq.bam"),val(type)
    //tuple path(r1),path(r2)
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${r1} \
    -F2 ${r2} \
    -SM ${sampleID} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    -O ${sampleID}.unmapped.from.fq.bam
    """
}

process tb_markAdapters {

    input:
    tuple val(caseID),val(sampleID), path(uBAM),val(type)
    
    output:
    tuple val(caseID),val(sampleID), path("${sampleID}.ubamXT.bam"), path("${sampleID}.markAdapterMetrics.txt"), val(type)
    
    script:

    """
    ${gatk_exec} MarkIlluminaAdapters \
    -I ${uBAM} \
    -O ${sampleID}.ubamXT.bam \
    --TMP_DIR ${tmpDIR} \
    -M ${sampleID}.markAdapterMetrics.txt
    """


}

process tb_align {
    tag "$sampleID"

    maxForks 5
    errorStrategy 'ignore'
    cpus 20

    input:
    tuple val(caseID),val(sampleID), path(uBAM), path(metrics),val(type)

    output:
    tuple val(caseID),val(sampleID), path("${sampleID}.${type}.${params.genome}.${genome_version}.QNsort.BWA.clean.bam"),val(type)
    
    script:
    """
    ${gatk_exec} SamToFastq \
    -I ${uBAM} \
    -INTER \
    -CLIP_ATTR XT \
    -CLIP_ACT 2 \
    -NON_PF \
    -F /dev/stdout \
    |  singularity run -B ${s_bind} ${simgpath}/bwa0717.sif bwa mem \
    -t ${task.cpus} \
    -p \
    ${genome_fasta} \
    /dev/stdin \
    | ${gatk_exec} MergeBamAlignment \
    -R ${genome_fasta} \
    -UNMAPPED ${uBAM} \
    -ALIGNED /dev/stdin \
    -MAX_GAPS -1 \
    -ORIENTATIONS FR \
    -SO queryname \
    -O ${sampleID}.${type}.${params.genome}.${genome_version}.QNsort.BWA.clean.bam
    """
}

/*
process tb_markDup_v2_bam_cram {
    errorStrategy 'ignore'
    maxForks 6
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    
    input:
    tuple val(caseID),val(sampleID), path(bam), val(type)
    
    output:
    tuple val(caseID), val(sampleID), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam"), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD*bai"), val(type), emit: markDup_bam

    tuple val(caseID),val(sampleID),  path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram"), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD*crai"),val(type), emit: markDup_cram


    script:
    """
    samtools view -h ${bam} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam /dev/stdin
    sambamba index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam

    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam

    samtools index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram

    """
}
*/

process markDup_cram {
    errorStrategy 'ignore'
    maxForks 6
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/CRAM/", mode: 'copy', pattern: "*.BWA.MD.cr*"
    publishDir "${caseID}/${outputDir}/variantcalls/Alignment_symlinks/", mode: 'link', pattern: "*.BWA.MD.cr*"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    input:
    tuple val(caseID),val(sampleID), path(bam), val(type)
    
    output:
    tuple val(caseID),val(sampleID),  path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram"), path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD*crai"), val(type), emit: markDup_output
    
    script:
    """
    samtools view -h ${bam} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=/data/TMP/TMP.${user}/ -o /dev/stdout /dev/stdin \
    |  samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram -

    samtools index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.cram
    """
}

/*
process tb_cram_bam {

   // publishDir "${caseID}/${params.outdir}/cram_TN_symlinks/", mode: 'link', pattern: '*.symedit*'
    input:
    tuple val(caseID), val(sampleID), path(cram), path(crai),val(type)

    output:
    tuple val(caseID), val(sampleID), path("*.bam"), path("*.bai"), val(type), emit:bam
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram")
    path("${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai")
    script:
    """

    samtools view \
    -b \
    -T ${genome_fasta} \
    -o ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam ${cram}

    samtools index ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.bam
    mv ${cram} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram
    mv ${crai} ${sampleID}.${type}.${params.genome}.${genome_version}.BWA.MD.symedit.cram.crai
    """
}
*/

///////////////////////////////// QC MODULES ////////////////////////////
// TB QC input channel structure:
// tuple val(caseID), val(sampleID),  path(aln), path(aln_index),val(type)

process tb_samtools {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/QC/", mode: 'copy'

    input:
    tuple val(caseID), val(sampleID),  path(bam), path(bai),val(type) 
    output:
    path("${sampleID}.samtools.sample.stats.txt")

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/samtools.sif samtools stats \
    ${bam} > ${sampleID}.samtools.sample.stats.txt
    """
}
/*
process tb_qualimap {
    errorStrategy 'ignore'
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/QC/bamQC", mode: 'copy'

    cpus 10
    input:
    tuple val(caseID), val(sampleID),  path(bam), path(bai),val(type)
    //path(targetBED) from ch_qualimap_ROI
    output:
    path ("${sampleID}/")

    script:
    use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''

    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif \
    qualimap --java-mem-size=5G bamqc \
    -nt ${task.cpus} \
    -outdir ${sampleID} \
    -bam ${bam} $use_bed -sd -sdmode 0
    """
}

process tb_fastqc_bam {
    errorStrategy 'ignore'
    tag "$sampleID"
    cpus 1
    publishDir "${caseID}/${outputDir}/QC/", mode: 'copy',saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }
    input:
    tuple val(caseID), val(sampleID),  path(bam), path(bai),val(type)
    
    output:
    path "*_fastqc.{zip,html}"
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/fastqc.sif --quiet --threads ${task.cpus} ${bam}
    """
}
*/


process multiQC {
    errorStrategy 'ignore'
    publishDir "${launchDir}", mode: 'copy'
    //publishDir "${launchDir}/*/${params.outdir}", mode: 'copy'
    input:
    path("*_fastqc.*") 
    path("${sn}.${sampleID_type}.samtools.sample.stats.txt")
    path("bamQC/*")
    output:
    path ("*.multiQC.report.html")
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -n ${date}.TN_WES.multiQC.report.html \
    -f -q  ${launchDir}/*/${outputDir}/QC/
    """

}



/////////////////////////////// VARIANT CALLING MODULES ///////////////////////

process tb_haplotypecaller {
    errorStrategy 'ignore'
    cpus 4
    if (params.server=="lnx01") {
    maxForks 2
    }
    tag "$sampleID"
    publishDir "${caseID}/${outputDir}/variantcalls/", mode: 'copy', pattern: "*.haplotypecaller.*"
    publishDir "${caseID}/${outputDir}/variantcalls/gvcf/", mode: 'copy', pattern: "*.g.*"
    publishDir "${variantStorage}/gVCF/tumorBoard/", mode: 'copy', pattern:'*.g.vc*' //

    input:
    tuple val(caseID), val(sampleID), path(cram), path(crai),val(type) 
    
    output:
    tuple val(caseID), val(sampleID),  path("${caseID}.${sampleID}.${type}.haplotypecaller.vcf.gz"), path("${caseID}.${sampleID}.${type}.haplotypecaller.vcf.gz.tbi"),emit: sample_gvcf
    path("${caseID}.${sampleID}.${type}.g.*")
    path("${crai}")
    path("${cram}")
    
    script:
    def datatype=params.wgs ? "": "-L ${ROI}"
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" HaplotypeCaller \
    -I ${cram} \
    -R ${genome_fasta} \
    -ERC GVCF \
    -L ${ROI} \
    --smith-waterman FASTEST_AVAILABLE \
    --native-pair-hmm-threads 30 \
    -pairHMM FASTEST_AVAILABLE \
    --dont-use-soft-clipped-bases \
    -O ${caseID}.${sampleID}.${type}.g.vcf.gz 
    
    ${gatk_exec} GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${caseID}.${sampleID}.${type}.g.vcf.gz \
    -O ${caseID}.${sampleID}.${type}.haplotypecaller.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
    """
}
// NB: Disabled full wGS germline var calling for now (too slow - requires implementation of splitIntervals: $datatype \)
process mutect2 {
    tag "$caseID"
    if (params.server=="lnx01") {
    maxForks 2
    }
    publishDir "${caseID}/${outputDir}/variantcalls/", mode: 'copy'
    //publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*.for.VarSeq.*"
    publishDir "${caseID}/${outputDir}/QC/mutect2_filtering/", mode: 'copy', pattern: "*.{table,stats,tsv}"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    input:
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
 
    output:

    tuple val(caseID), path("${caseID}.mutect2.for.VarSeq.vcf.gz"), path("${caseID}.mutect2.for.VarSeq.vcf.gz.tbi"), emit: mutect2_ALL

    tuple val(caseID), val(sampleID_tumor), path("${caseID}.mutect2.PASSonly.vcf.gz"), path("${caseID}.mutect2.PASSonly.vcf.gz.tbi"),emit: mutect2_PASS 
  
    tuple val(caseID), path("${caseID}.mutect2.PASSonly.TUMORonly.vcf.gz"), path("${caseID}.mutect2.PASSonly.TUMORonly.vcf.gz.tbi"),emit: mutect2_tumorPASS 

    tuple val(caseID), path("${caseID}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz"), path("${caseID}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz.tbi"), emit: mutect2_PASS_TMB_filtered 

    tuple val(caseID), path("${caseID}.mutect2.PASSonly.snpeff.vcf"), emit: mutect2_snpEFF 

    // for HRD:
    tuple val(caseID), path("${caseID}.mutect2.PASSonly.chr1_22_XY.vcf.gz"), path("${caseID}.mutect2.PASSonly.chr1_22_XY.vcf.gz.tbi"), emit: mutect2_PASS_reduced

    script:
    //def datatype=params.wgs ? "-L ${ROI}" : "-L ${ROI}"
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" Mutect2 \
    -R ${genome_fasta} \
    -I ${bamT} \
    -I ${bamN} \
    -normal ${sampleID_normal} \
    --germline-resource ${mutect_gnomad} \
    --panel-of-normals ${gatk_wgs_pon} \
    -L ${ROI} \
    --dont-use-soft-clipped-bases \
    --native-pair-hmm-threads 30 \
    -pairHMM FASTEST_AVAILABLE \
    --smith-waterman FASTEST_AVAILABLE \
    -O ${caseID}.mutect2.raw.vcf.gz \
    --f1r2-tar-gz ${caseID}.within.f1r2.tar.gz
    
    ${gatk_exec} LearnReadOrientationModel \
    -I ${caseID}.within.f1r2.tar.gz \
    -O ${caseID}.within.ROmodel.tar.gz
    
    ${gatk_exec} GetPileupSummaries -I ${bamT} \
    -R ${genome_fasta} \
    -V ${gatk_contamination_ref} \
    -L ${gatk_contamination_ref} \
    -O ${caseID}.within.getpileupsummaries.table
    
    ${gatk_exec} CalculateContamination \
    -I ${caseID}.within.getpileupsummaries.table \
    -tumor-segmentation ${caseID}.segments.table \
    -O ${caseID}.contamination.table
    
    ${gatk_exec} FilterMutectCalls \
    -V ${caseID}.mutect2.raw.vcf.gz \
    -R ${genome_fasta} \
    --tumor-segmentation ${caseID}.segments.table \
    --contamination-table ${caseID}.contamination.table \
    --min-allele-fraction 0.001 \
    -O ${caseID}.mutect2.for.VarSeq.vcf.gz
    
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${caseID}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered \
    -O ${caseID}.mutect2.PASSonly.vcf.gz

    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${caseID}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered \
    -L /data/shared/genomes/hg38/chrom1_22_XY.1col.list \
    -O ${caseID}.mutect2.PASSonly.chr1_22_XY.vcf.gz
    
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${caseID}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered -xl-sn ${sampleID_normal} --exclude-non-variants \
    -O ${caseID}.mutect2.PASSonly.TUMORonly.vcf.gz

    java -jar /data/shared/programmer/snpEff5.2/snpEff.jar GRCh38.99 ${caseID}.mutect2.PASSonly.vcf.gz > ${caseID}.mutect2.PASSonly.snpeff.vcf

    cat ${caseID}.mutect2.PASSonly.snpeff.vcf | java -jar /data/shared/programmer/snpEff5.2/SnpSift.jar filter \
    "(ANN[0].EFFECT has 'missense_variant'| ANN[0].EFFECT has 'frameshift_variant'| ANN[0].EFFECT has 'stop_gained'| ANN[0].EFFECT has 'conservative_inframe_deletion'|  ANN[0].EFFECT has 'disruptive_inframe_deletion'|ANN[0].EFFECT has 'disruptive_inframe_insertion'|ANN[0].EFFECT has 'conservative_inframe_insertion') & (GEN[${sampleID_tumor}].AF >=0.01 & GEN[${sampleID_tumor}].DP>25)" > ${caseID}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf

    bgzip ${caseID}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf
    bcftools index -t ${caseID}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz

    """
}



process strelka2 {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/variantcalls/strelka2", mode: 'copy'
    cpus 10

    if (params.server=="lnx01"){
            conda '/data/shared/programmer/miniconda3/envs/py310'
    }
    
    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

    output:
    path("*.strelka2.*")
    tuple val(caseID), path("${caseID}.strelka2.merged.vaf.vcf.gz"), emit: strelkarenameVCF

    script:
    def datatype=params.wgs ? "": "--exome"
    """
    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif /tools/strelka2/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${bamN} \
    --tumorBam ${bamT} \
    --referenceFasta  ${genome_fasta} \
    $datatype \
    --runDir strelka

    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif python2 strelka/runWorkflow.py \
    -j ${task.cpus} \
    -m local

    python /data/shared/programmer/VCFpytools/add_vaf_strelka2.py \
    --input strelka/results/variants/somatic.indels.vcf.gz \
    --output ${caseID}.strelka2.indels.vaf.vcf \
    --variant indel

    python /data/shared/programmer/VCFpytools/add_vaf_strelka2.py \
    --input strelka/results/variants/somatic.snvs.vcf.gz \
    --output ${caseID}.strelka2.snvs.vaf.vcf \
    --variant snv

    ${gatk_exec} MergeVcfs \
    -I ${caseID}.strelka2.snvs.vaf.vcf \
    -I ${caseID}.strelka2.indels.vaf.vcf \
    -O ${caseID}.strelka2.merged.vaf.vcf.gz
    """
}

process strelka2_edits {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/variantcalls/strelka2", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN), val(typeN), val(sampleID_tumor),path(bamT), path(baiT), val(typeT), path(strelkavcf)

    output:
    
    tuple val(caseID), path("${caseID}.strelka2.for.VarSeq.gz"),path("${caseID}.strelka2.for.VarSeq.gz.tbi"),emit: strelka2_ALL    
    
    tuple val(caseID), path("${caseID}.strelka2.PASSonly.vcf.gz"),path("${caseID}.strelka2.PASSonly.vcf.gz.tbi"), emit: strelka2_PASS 
    
    tuple val(caseID), path("${caseID}.strelka2.PASSonly.TUMORonly.vcf.gz"),path("${caseID}.strelka2.PASSonly.TUMORonly.vcf.gz.tbi"), emit: strelka2_TUMOR_PASS

    tuple val(caseID), path("${caseID}.strelka2.PASSonly.snpeff.vcf"), emit: strelka2_PASS_snpeff
    
    tuple val(caseID), path("${caseID}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz"),path("${caseID}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz.tbi"), emit: strelka2_PASS_TMB_filtered

    path("*.strelka2.*")
    
    shell:
    '''
    printf "TUMOR !{sampleID_tumor}_TUMOR" >> !{caseID}.strelka_rename.txt
    printf "\nNORMAL !{sampleID_normal}_NORMAL" >> !{caseID}.strelka_rename.txt

    bcftools reheader \
    --samples !{caseID}.strelka_rename.txt \
    -o !{caseID}.strelka2.for.VarSeq.gz !{strelkavcf}

    bcftools index -t !{caseID}.strelka2.for.VarSeq.gz
    
    !{gatk_exec} SelectVariants -R !{genome_fasta} \
    -V !{caseID}.strelka2.for.VarSeq.gz \
    --exclude-filtered \
    -O !{caseID}.strelka2.PASSonly.vcf.gz

    !{gatk_exec} SelectVariants -R !{genome_fasta} \
    -V !{caseID}.strelka2.for.VarSeq.gz \
    --exclude-filtered -xl-sn !{sampleID_normal}_NORMAL --exclude-non-variants \
    -O !{caseID}.strelka2.PASSonly.TUMORonly.vcf.gz

    java -jar /data/shared/programmer/snpEff5.2/snpEff.jar GRCh38.99 !{caseID}.strelka2.PASSonly.vcf.gz > !{caseID}.strelka2.PASSonly.snpeff.vcf

    cat !{caseID}.strelka2.PASSonly.snpeff.vcf | java -jar /data/shared/programmer/snpEff5.2/SnpSift.jar filter \
    "(ANN[0].EFFECT has 'missense_variant'| ANN[0].EFFECT has 'frameshift_variant'| ANN[0].EFFECT has 'stop_gained'| ANN[0].EFFECT has 'conservative_inframe_deletion'|  ANN[0].EFFECT has 'disruptive_inframe_deletion'|ANN[0].EFFECT has 'disruptive_inframe_insertion'|ANN[0].EFFECT has 'conservative_inframe_insertion') & (GEN[!{sampleID_tumor}_TUMOR].VAF >=0.05 & GEN[!{sampleID_tumor}_TUMOR].DP>25 & GEN[!{sampleID_normal}_NORMAL].VAF<0.001)" > !{caseID}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf
    
    bgzip !{caseID}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf
    bcftools index -t !{caseID}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz
    
    '''
}


process msisensor {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/MSIsensor/", mode: 'copy'
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*_msi"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/msisensorPro120'

    if (params.server=="lnx01") {
    maxForks 2
    }
    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT), val(typeT)

    output:
    path("*_msi*")
 
    script:
    def datatype=params.wgs ? "": "-e ${ROI}"
    """
    msisensor-pro msi \
    -d ${msisensor_list} \
    -n ${bamN} -t ${bamT} \
    $datatype \
    -g ${genome_fasta} \
    -o ${caseID}_msi
    """
}

process sequenza_conda {
    errorStrategy 'ignore'
    tag "$caseID"
    if (params.server=="lnx01") {
    maxForks 2
    }

    if (params.server=="rgi01") {
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenza30'
    }

    input:
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    output:
    tuple val(caseID), path("${caseID}.seqz.final.gz") 
    
    script:
    """
    sequenza-utils bam2seqz \
    -n ${bamN} -t ${bamT} \
    --fasta ${genome_fasta} \
    -gc ${sequenza_cg50_wig} \
    -o ${caseID}.seqz.phase1.gz
    sequenza-utils seqz_binning --seqz ${caseID}.seqz.phase1.gz \
    -w 50 -o ${caseID}.seqz.final.gz 
    """
}

process sequenza_R_output_conda {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/", mode: 'copy'
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*_{segments,alternative_fit,genome_view}.{txt,pdf}"
    publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "sequenza_conda/*_{alternative_fit,genome_view}.pdf"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenzaEnv'
    input:
    tuple val(caseID),  path(seqz)

    output:
    path("sequenza_conda/*")

    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    t1 = sequenza.extract("${seqz}",verbose=F)
    cp = sequenza.fit(t1)
    sequenza.results(sequenza.extract = t1, cp.table = cp, sample.id = "${caseID}", out.dir = "sequenza_conda" )
    """
}

process sequenza_R_output_conda_editPARAMS {
    errorStrategy 'ignore'
    tag "$caseID"
    publishDir "${caseID}/${outputDir}/", mode: 'copy'
    //publishDir "${caseID}/${outputDir}/tumorBoard_files/", mode: 'copy', pattern: "*_{segments,alternative_fit,genome_view}.{txt,pdf}"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenzaEnv'
    input:
    tuple val(caseID),  path(seqz)

    output:
    path("sequenza_conda_editPARAMS/*")

    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    t1 = sequenza.extract("${seqz}",verbose=F)
    cp = sequenza.fit(t1, segment.filter=1e6)
    sequenza.results(sequenza.extract = t1, cp.table = cp, sample.id = "${caseID}", out.dir = "sequenza_conda_editPARAMS", CNt.max=1000)
    """
}
/*
process pcgr_v141 {
    tag "$caseID"
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR141/", mode: 'copy', pattern: "*.pcgr_acmg.*"
    //publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr_acmg.*")
    
    script:
    def datatype=params.wgs ? "--assay WGS": "--assay WES"
    //tumorsite=${pcgr_tumor} ? "--tumor_site ${pcgr_tumor}" : ""
    """
    singularity run -B ${s_bind} ${simgpath}/pcgr141.sif pcgr \
    --input_vcf ${vcf} \
    --pcgr_dir ${pcgr_data_dir} --output_dir . \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_TMB_NonSyn \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb --estimate_msi_status \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    --tumor_site ${pcgr_tumor} \
    --estimate_signatures \
    --tmb_algorithm nonsyn \
    --include_trials
    """
}


process pcgr_v203_mutect2 {
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR203/mutect2/", mode: 'copy', pattern: "*.pcgr.*"
    //publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
   
  
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr203'


    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr.*.{xlsx,tsv,html}")
    
    script:
    def datatype=params.wgs ? "--assay WGS": "--assay WES"
    """
    pcgr \
    --input_vcf ${vcf} \
    --refdata_dir  ${pcgr_data_dir2} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_mutect2 \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display missense_only \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    --tumor_site ${pcgr_tumor} \
    --estimate_signatures
    """
}
*/

process pcgr_v212_strelka2 {
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR212/strelka2/", mode: 'copy', pattern: "*.pcgr.*"
    //publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
   
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr212'  
   
    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr.*.{xlsx,tsv,html}")
    
    script:
    
    def datatype=params.wgs ? "--assay WGS": "--assay WES"
    """
    pcgr \
    --input_vcf ${vcf} \
    --refdata_dir  ${pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_strelka2 \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display missense_only \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    --tumor_site ${pcgr_tumor} \
    --pcgrr_conda /lnx01_data3/shared/programmer/miniconda3/envs/pcgrr212 \
    --estimate_signatures
    """
}

process pcgr_v212_strelka2_manualFilter {
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR212/strelka2_manual", mode: 'copy', pattern: "*.pcgr.*"
    //publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
   
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr212'  
   
    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr.*.{xlsx,tsv,html}")
    
    script:
    
    def datatype=params.wgs ? "--assay WGS": "--assay WES"
    """
    pcgr \
    --input_vcf ${vcf} \
    --refdata_dir  ${pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${caseID}_strelka2_manual \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display missense_only \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    --tumor_site ${pcgr_tumor} \
    --pcgrr_conda /lnx01_data3/shared/programmer/miniconda3/envs/pcgrr212 \
    --estimate_signatures
    """
}

process pcgr_v212_mutect2 {
    errorStrategy 'ignore'
    publishDir "${caseID}/${outputDir}/PCGR212/mutect2/", mode: 'copy', pattern: "*.pcgr.*"
    publishDir "${caseID}/${outputDir}/tumorBoard_files/",mode: 'copy', pattern:"*.html"
    //publishDir "${caseID}/${params.outdir}/tumorBoard_files/", mode: 'copy', pattern: "*.flexdb.html"
   
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr212'  
   
    input:
    tuple val(caseID),  path(vcf), path(idx), val(pcgr_tumor)

    output:
    path("*.pcgr.*")
    
    script:
    def tumorsite=params.pcgrtumor ? "--tumor_site ${params.pcgrtumor}" : ""
    def rnaexp=params.rnaExp ? "--input_rna_expression ${rna}" : ""
    """

    pcgr \
    --input_vcf ${vcf} \
    --refdata_dir  ${pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${grch_assembly} \
    --sample_id ${caseID}_pcgr212 \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display coding_non_silent \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    --assay WES \
    --pcgrr_conda /lnx01_data3/shared/programmer/miniconda3/envs/pcgrr212 \
    --estimate_signatures \
    --tumor_site ${pcgr_tumor}
    """
//bcftools index -f -t ${vcf}
}


/////// WGS ONLY PROCESSES ////////////////

process cnvkit_somatic {
    errorStrategy 'ignore'
    tag "$caseID"

    cpus 12
    if (params.server=="lnx01") {
    maxForks 2
    }

    publishDir "${caseID}/${outputDir}/NEWTOOLS/cnvkit_somatic/", mode: 'copy'

    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

    output:
    path("${caseID}.cnvkit/*")
    //path("*.targetcoverage.cnn"), emit: cnvkit_cnn_out
    tuple val(caseID), path("${caseID}.cnvkit/*.call.cns"), emit: CNVcalls
    tuple val(caseID), path("${caseID}.cnvkit/*.cnr"), emit: CNVcnr
    //path("${caseID}.cnvkit/*.cnn")
    
    // touch ${index}
    script:
    """
    mv ${baiN} intermediate_crai
    cp intermediate_crai ${baiN}
    rm intermediate_crai

    mv ${baiT} intermediate_crai2
    cp intermediate_crai2 ${baiT}
    rm intermediate_crai2

    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py batch \
    ${bamT} \
    -n ${bamN} \
    -m wgs \
    -p ${task.cpus} \
    -f ${genome_fasta} \
    --annotate ${refFlat} \
    --scatter --diagram \
    -d ${caseID}.cnvkit/
    cp ${caseID}.cnvkit/*.targetcoverage.cnn .
    """
}

process cnvkitExportFiles {
    errorStrategy 'ignore'
    tag "$caseID"
   // publishDir "${inhouse_SV}/CNVkit/raw_calls/", mode: 'copy', pattern: '*.cnvkit.vcf'
    publishDir "${caseID}/${outputDir}/NEWTOOLS/cnvkit_somatic/", mode: 'copy'

    input:
    tuple val(caseID), path(cnvkit_calls)// from cnvkit_calls_out
    tuple val(caseID), path(cnvkit_cnr)// from cnvkit_cnr_out

    output:
    path("*.vcf")
    path("*.seg")
    tuple val("${caseID}"), path("${caseID}.cnvkit.somatic.vcf"), emit: cnvkitForSVDB

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py export vcf \
    ${cnvkit_calls} \
    -i ${caseID} \
    -o ${caseID}.cnvkit.somatic.vcf

    singularity run -B ${s_bind} ${simgpath}/cnvkit.sif cnvkit.py export seg \
    ${cnvkit_cnr} \
    -o ${caseID}.cnvkit.cnr.somatic.seg
    """
}

process manta_somatic {
    errorStrategy 'ignore'
    tag "$caseID"

    publishDir "${caseID}/${outputDir}/NEWTOOLS/manta_somatic/", mode: 'copy'
    //publishDir "${outputDir}/structuralVariants/manta/", mode: 'copy', pattern: "*.{AFanno,filtered}.*"

    cpus 12

    if (params.server=="lnx01") {
    maxForks 1
    }
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 
    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

    output:
    tuple val(caseID), path("${caseID}.manta.*.{vcf,vcf.gz,gz.tbi}")
    tuple val(caseID), path("${caseID}.manta.somaticSV.vcf.gz"), emit: mantaSV_all
    tuple val(caseID), path("${caseID}.manta.somaticSV.PASSonly.vcf.gz"), emit: mantaSV_pass
    //    tuple val(caseID), path("${caseID}.manta.somaticSV.PASSonly.Inhouse127.vcf.gz"), emit: mantaSV_pass_inhouse
     //  tuple val(caseID), path("${caseID}.manta.somaticSV.bcftools.Inhouse127.vcf.gz"), emit: mantaSV_inhouse

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif configManta.py \
    --normalBam ${bamN} \
    --tumorBam ${bamT} \
    --referenceFasta ${genome_fasta} \
    --callRegions ${manta_callable_regions} \
    --runDir manta

    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif ./manta/runWorkflow.py -j ${task.cpus}

    mv manta/results/variants/candidateSmallIndels.vcf.gz \
    ${caseID}.manta.candidateSmallIndels.vcf.gz
    
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    ${caseID}.manta.candidateSmallIndels.vcf.gz.tbi
    
    mv manta/results/variants/candidateSV.vcf.gz \
    ${caseID}.manta.candidateSV.vcf.gz
    
    mv manta/results/variants/candidateSV.vcf.gz.tbi \
    ${caseID}.manta.candidateSV.vcf.gz.tbi

    mv manta/results/variants/diploidSV.vcf.gz \
    ${caseID}.manta.diploidSV.vcf.gz
    
    mv manta/results/variants/diploidSV.vcf.gz.tbi \
    ${caseID}.manta.diploidSV.vcf.gz.tbi

    gzip -dc ${caseID}.manta.diploidSV.vcf.gz > ${caseID}.manta.diploidSV.vcf

    mv manta/results/variants/somaticSV.vcf.gz \
    ${caseID}.manta.somaticSV.vcf.gz
    
    mv manta/results/variants/somaticSV.vcf.gz.tbi \
    ${caseID}.manta.somaticSV.vcf.gz.tbi


    bcftools view \
    -i 'FILTER="PASS" | FILTER="."' \
    ${caseID}.manta.somaticSV.vcf.gz > ${caseID}.manta.somaticSV.PASSonly.vcf

    bgzip ${caseID}.manta.somaticSV.PASSonly.vcf
    bcftools index -t ${caseID}.manta.somaticSV.PASSonly.vcf.gz


    """
}
/*
    bcftools filter \
    -R ${inhouse127_geneIntervals} \
    -o ${caseID}.manta.somaticSV.bcftools.Inhouse127.vcf ${caseID}.manta.somaticSV.vcf.gz

    bgzip ${caseID}.manta.somaticSV.bcftools.Inhouse127.vcf
    bcftools index -t ${caseID}.manta.somaticSV.bcftools.Inhouse127.vcf.gz
*/
//bcftools filter -R {inhouse127_geneIntervals} -o ${caseID}.manta.somaticSV.PASSonly.bcftools.Inhouse127.vcf.gz ${caseID}.manta.somaticSV.vcf.gz

process hrd_scores_fullSV {
    errorStrategy 'ignore'
    tag "$caseID"
  
    publishDir "${caseID}/${outputDir}/NEWTOOLS/HRD_FULLSV/", mode: 'copy'
    publishDir "${caseID}/${outputDir}/tumorBoard_files/",mode: 'copy', pattern: "*.txt"
    cpus 4
    maxForks 3
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sigrap011/'

    input:
    tuple val(caseID), path(snv), path(cnv), path(sv)

    output:
    tuple val(caseID), path("${caseID}.HRD_SCORES.txt")
    path("*.txt")
    script:
    """
    #!/usr/bin/env Rscript
    library(sigrap)
    library(BSgenome.Hsapiens.UCSC.hg38)
 
    chord=sigrap::chord_run(vcf.snv="${snv}",vcf.sv="${sv}", sv.caller="manta", sample.name="${caseID}", ref.genome="hg38")
    hrdetect=sigrap::hrdetect_run(snvindel_vcf="${snv}",sv_vcf="${sv}", nm="${caseID}", cnv_tsv="${cnv}", genome="hg38")
   
    c1=as.data.frame(chord[2])
    c2=c1[,1:6]
    names(c2)=c("sample","CHORD_p_hrd","CHORD_hr_status","CHORD_hrdtype","CHORD_pBRCA1", "CHORD_pBRCA2")

    d1=data.frame(hrdetect[1], hrdetect[2])
    names(d1)=c("sample","HRDETECT_p_hrd")

    m1=merge(c2,d1,by="sample")
    if (m1[["HRDETECT_p_hrd"]]>0.7) m1[["HRDETECT_verdict"]]="HR_DEFICIENT" else m1[["HRDETECT_verdict"]]="HR_proficient"
    write.table(m1,file="${caseID}.HRD_SCORES.txt",sep="\t",quote=F,row.names=F)

    """
}

process hrd_scores_PASS {
    errorStrategy 'ignore'
    tag "$caseID"
  
    publishDir "${caseID}/${outputDir}/NEWTOOLS/HRD_PASSvariants/", mode: 'copy'

    cpus 4
    maxForks 3
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sigrap011/'

    input:
    tuple val(caseID), path(snv), path(cnv), path(sv)

    output:
    tuple val(caseID), path("${caseID}.HRD_SCORES.txt")
    path("*.txt")
    script:
    """
    #!/usr/bin/env Rscript
    library(sigrap)
    library(BSgenome.Hsapiens.UCSC.hg38)
 
    chord=sigrap::chord_run(vcf.snv="${snv}",vcf.sv="${sv}", sv.caller="manta", sample.name="${caseID}", ref.genome="hg38")
    hrdetect=sigrap::hrdetect_run(snvindel_vcf="${snv}",sv_vcf="${sv}", nm="${caseID}", cnv_tsv="${cnv}", genome="hg38")
   
    c1=as.data.frame(chord[2])
    c2=c1[,1:6]
    names(c2)=c("sample","CHORD_p_hrd","CHORD_hr_status","CHORD_hrdtype","CHORD_pBRCA1", "CHORD_pBRCA2")

    d1=data.frame(hrdetect[1], hrdetect[2])
    names(d1)=c("sample","HRDETECT_p_hrd")

    m1=merge(c2,d1,by="sample")
    if (m1[["HRDETECT_p_hrd"]]>0.7) m1[["HRDETECT_verdict"]]="HR_DEFICIENT" else m1[["HRDETECT_verdict"]]="HR_proficient"
    write.table(m1,file="${caseID}.HRD_SCORES.txt",sep="\t",quote=F,row.names=F)

    """
}


process cobalt {
    publishDir "${caseID}/${outputDir}/NEWTOOLS/cobalt_amber_sage_purple/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$caseID"

    if (params.server=="lnx01") {
    maxForks 1
    }
    cpus 16
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/hmftools/'

    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    
    output: 
    tuple val(caseID), path("${caseID}_cobaltOut/")
    
    script:
    """
    cobalt "-Xmx16G" \
    -reference ${sampleID_normal} \
    -reference_bam ${bamN} \
    -tumor ${sampleID_tumor} \
    -tumor_bam ${bamT} \
    -ref_genome ${genome_fasta} \
    -output_dir ${caseID}_cobaltOut \
    -threads ${task.cpus} \
    -gc_profile ${hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp
    """
}

process amber {
    publishDir "${caseID}/${outputDir}/NEWTOOLS/cobalt_amber_sage_purple/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$caseID"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/hmftools/'

    cpus 12
    
    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)

    output:
    tuple val(caseID), val(sampleID_normal), val(sampleID_tumor), path("${caseID}_amberOut/")
    
    script:
    """
    amber "-Xmx16G" \
    -reference ${sampleID_normal} \
    -reference_bam ${bamN} \
    -tumor ${sampleID_tumor} \
    -tumor_bam ${bamT} \
    -ref_genome ${genome_fasta} \
    -output_dir ${caseID}_amberOut \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -loci ${hmftools_data_dir_v534}/copy_number/GermlineHetPon.38.vcf.gz
    """
}

process sage {
    publishDir "${caseID}/${outputDir}/NEWTOOLS/cobalt_amber_sage_purple/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$caseID"

    if (params.server=="lnx01") {
    maxForks 2
    }
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 
    cpus 16
    input: 
    tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    
    output: 
    tuple val(caseID), path("${caseID}_sage.somatic.SNV_INDELS.PASS.vcf.gz"),emit:sage_pass
    tuple val(caseID), path("${caseID}_sage.somatic.SNV_INDELS.vcf.gz"),emit:sage_all
    script:
    """
    java -jar /data/shared/programmer/hmftools/sage_v3.4.4.jar \
    -reference ${sampleID_normal} \
    -reference_bam ${bamN} \
    -tumor ${sampleID_tumor} \
    -tumor_bam ${bamT} \
    -ref_genome ${genome_fasta} \
    -ref_genome_version 38 \
    -threads ${task.cpus} \
    -ensembl_data_dir ${hmftools_data_dir_v534}/common/ensembl_data/ \
    -high_confidence_bed ${hmftools_data_dir_v534}/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz \
    -hotspots ${hmftools_data_dir_v534}/variants/KnownHotspots.somatic.38.vcf.gz \
    -panel_bed ${hmftools_data_dir_v534}/variants/ActionableCodingPanel.38.bed.gz \
    -output_vcf ${caseID}_sage.somatic.SNV_INDELS.vcf.gz

    bcftools view \
    -i 'FILTER="PASS" | FILTER="."' \
    ${caseID}_sage.somatic.SNV_INDELS.vcf.gz > ${caseID}_sage.somatic.SNV_INDELS.PASS.vcf

    bgzip ${caseID}_sage.somatic.SNV_INDELS.PASS.vcf
    bcftools index -t ${caseID}_sage.somatic.SNV_INDELS.PASS.vcf.gz

    """
}
process purple_full {
    publishDir "${caseID}/${outputDir}/NEWTOOLS/cobalt_amber_sage_purple/PURPLE_FULL", mode: 'copy'
    publishDir "${caseID}/${outputDir}/tumorBoard_files/",mode: 'copy', pattern: "*.{png,purity.tsv}"
    errorStrategy 'ignore'
    tag "$caseID"
    cpus 12

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/'

    input:
    tuple val(caseID), val(sampleID_normal), val(sampleID_tumor), path(amber),path(cobalt),path(manta_sv), path(sage)

    output:
    tuple val(caseID), path("${caseID}_purple/")
    tuple val(caseID), path("${caseID}.purple.cnv.somatic.tsv"), emit: purple_full_for_hrd
    tuple path("${caseID}.purple.qc"), path("${caseID}.purple.purity.tsv"),path("${caseID}.purple.circos.png")
    script:
    """
    purple "-Xmx16G" \
    -reference ${sampleID_normal} \
    -tumor ${sampleID_tumor} \
    -ref_genome ${genome_fasta} \
    -output_dir ${caseID}_purple \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -somatic_sv_vcf ${manta_sv} \
    -somatic_vcf ${sage} \
    -gc_profile ${hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp \
    -ensembl_data_dir ${hmftools_data_dir_v534}/common/ensembl_data/ \
    -amber ${amber} \
    -cobalt ${cobalt} \
    -circos /lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/bin/circos

    cp ${caseID}_purple/${sampleID_tumor}*.somatic.tsv ${caseID}.purple.cnv.somatic.tsv
    cp ${caseID}_purple/${sampleID_tumor}*.qc ${caseID}.purple.qc
    cp ${caseID}_purple/${sampleID_tumor}*.purity.tsv ${caseID}.purple.purity.tsv
    cp ${caseID}_purple/plot/${sampleID_tumor}.circos.png ${caseID}.purple.circos.png
    """
}

process purple_pass {
    publishDir "${caseID}/${outputDir}/NEWTOOLS/cobalt_amber_sage_purple/PURPLE_PASS/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$caseID"
    cpus 12

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/'

    input:
    tuple val(caseID), val(sampleID_normal), val(sampleID_tumor), path(amber),path(cobalt),path(manta_sv), path(sage)

    output:
    tuple val(caseID), path("${caseID}_purple/")
    tuple val(caseID), path("${caseID}.purple.cnv.somatic.tsv"), emit: purple_pass_for_hrd
    tuple path("${caseID}.purple.qc"), path("${caseID}.purple.purity.tsv"),path("${caseID}.purple.PASS.circos.png")
    script:
    """
    purple "-Xmx8G" \
    -reference ${sampleID_normal} \
    -tumor ${sampleID_tumor} \
    -ref_genome ${genome_fasta} \
    -output_dir ${caseID}_purple \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -somatic_sv_vcf ${manta_sv} \
    -somatic_vcf ${sage} \
    -gc_profile ${hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp \
    -ensembl_data_dir ${hmftools_data_dir_v534}/common/ensembl_data/ \
    -amber ${amber} \
    -cobalt ${cobalt} \
    -circos /lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/bin/circos

    cp ${caseID}_purple/${sampleID_tumor}*.somatic.tsv ${caseID}.purple.cnv.somatic.tsv
    cp ${caseID}_purple/${sampleID_tumor}*.qc ${caseID}.purple.qc
    cp ${caseID}_purple/${sampleID_tumor}*.purity.tsv ${caseID}.purple.purity.tsv
    cp ${caseID}_purple/plot/${sampleID_tumor}.circos.png ${caseID}.purple.PASS.circos.png
    """
//    purple "-Xmx16G" \
}


/////////// SUBWORKFLOWS

workflow SUB_DNA_PREPROCESS {

    take:
    case_fastq_input_ch     // caseid, sampleID, ,R1, R2, type
   
    main:
    inputFiles_symlinks_fq(case_fastq_input_ch)
    tb_fastq_to_ubam(case_fastq_input_ch)
    tb_markAdapters(tb_fastq_to_ubam.out)
    tb_align(tb_markAdapters.out)
    markDup_cram(tb_align.out)

    emit:
    finalAln=markDup_cram.out.markDup_output //caseID, sampleID, CRAM, CRAI,type

}

workflow SUB_DNA_QC {
    take:
    cram_per_sample_ch
 
    main:
    tb_samtools(cram_per_sample_ch)
   // tb_qualimap(cram_per_sample_ch)
    //tb_fastqc_bam(cram_per_sample_ch)
 //   multiQC(tb_samtools.out.collect())

}

workflow SUB_PAIRED_TN {
    take:
    tumorNormal_cram_ch
    caseID_pcgrID
    main:
    if (!params.hrdOnly) {
    mutect2(tumorNormal_cram_ch)
    strelka2(tumorNormal_cram_ch)
    strelka2_edits(tumorNormal_cram_ch.join(strelka2.out.strelkarenameVCF))
    //tumorNormal_bam_ch.join(strelka2.out.strelkarenameVCF).view()

    msisensor(tumorNormal_cram_ch)

    sequenza_conda(tumorNormal_cram_ch)
    sequenza_R_output_conda(sequenza_conda.out)
    sequenza_R_output_conda_editPARAMS(sequenza_conda.out)
   // pcgr_v141(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
  //  pcgr_v203_mutect2(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
    pcgr_v212_strelka2(strelka2_edits.out.strelka2_PASS.join(caseID_pcgrID))
    pcgr_v212_strelka2_manualFilter(strelka2_edits.out.strelka2_PASS_TMB_filtered.join(caseID_pcgrID))
    pcgr_v212_mutect2(mutect2.out.mutect2_tumorPASS.join(caseID_pcgrID))
    }

    if (params.wgs || params.hrdOnly) {
       // cnvkit_somatic(tumorNormal_cram_ch)
        //cnvkitExportFiles(cnvkit_somatic.out.CNVcalls, cnvkit_somatic.out.CNVcnr)
        manta_somatic(tumorNormal_cram_ch)
        amber(tumorNormal_cram_ch)
        cobalt(tumorNormal_cram_ch)
        sage(tumorNormal_cram_ch)

        amber.out.join(cobalt.out).join(manta_somatic.out.mantaSV_all).join(sage.out.sage_all)
        | set {purple_full_input}

        amber.out.join(cobalt.out).join(manta_somatic.out.mantaSV_pass).join(sage.out.sage_pass)
        | set {purple_pass_input}

        purple_full(purple_full_input)
        purple_pass(purple_pass_input)
       // purple_inhouse(purple_inhouse_input)
        
        sage.out.sage_all.join(purple_full.out.purple_full_for_hrd).join(manta_somatic.out.mantaSV_all)
        | set {hrd_full_input}

        sage.out.sage_pass.join(purple_pass.out.purple_pass_for_hrd).join(manta_somatic.out.mantaSV_pass)
        | set {hrd_PASS_input}

        hrd_scores_fullSV(hrd_full_input)
        hrd_scores_PASS(hrd_PASS_input)

    }
   // emit:    
    //mutect2_out=mutect2.out.mutect2_ALL
}
