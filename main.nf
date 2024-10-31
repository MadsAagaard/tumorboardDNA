#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


params.rundir                           ="${launchDir.baseName}" 
params.gatkTEMP                         ="${launchDir.baseName}/gatkTEMP"
params.server                           ="lnx01"
params.genome                           ="hg38" 
params.outdir                           ="TN_WES_results"
params.panel                            ="WES_2"    // set ROI to full WES
params.fastq                            =null
params.cram                             =null
params.fastqInput                       =null
params.help                             =false
params.pcgr_tumor                       =null
params.qualimap                         =null
params.samplesheet                      =null
params.hg38v1                           =null
params.hg38v2                           =null
params.skipQC                           =null
params.archiveStorage                   =null
params.keepwork                         =null
params.nomail                           =null
params.wgs                              =null
params.gatk                             ="new"
params.fastqTMB                         =null
params.cramTMB                          =null
params.NGC                              =null
params.assaytype                        =null
params.hrdOnly                          =null
//outdir_full_path= "${launchDir}/${params.outdir}/"

runtype = "paired_TN"


switch (params.server) {
    case 'lnx02':
     //   modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        dataArchive="/lnx01_data2/shared/dataArchive";        
    break;
    case 'lnx01':
     //   modules_dir="/home/mmaj/scripts_lnx01/nextflow_lnx01/dsl2/modules";
        dataArchive="/lnx01_data2/shared/dataArchive";        
    break;
}
/*
if (params.wgs) {
    datapattern="WG4_NGC,WG4,WG3"
}

if (!params.wgs) {
    datapattern="EV8,EV8_BEH"
}
*/

switch (params.assaytype) {

    case 'NGC':
    datapattern="WG4_NGC"
    break;

    case 'BEH':
    datapattern="EV8_BEH"
    break;
    default:
    datapattern=""
    break;
}

def helpMessage() {
    log.info"""

    Generel info:
    Requires a samplesheet containing 5 columns in specific order (tab separated), without headerline:
    1) caseID, 2) NPN normal WES, 3) NPN tumor WES, 4) NPN tumor RNA, 5) PCGR tumor value

    Example samplesheet:

    johnDoe 112217976652	111184925465	111184925473    23

    The script will automatically look for fastq or cram files in subfolders at /lnx01_data2/shared/dataArchive/. This location contains read-only access to the data archive. Theres no need to copy or move any input data.

    The user can point to a specific folder containing input data using the --fastq or --cram option. 

    This is only needed if input data exists outside the data archive (e.g. if data are in personal folders or stored at other KG Vejle servers).

    Usage:

    Main options:
      --help                print this help message
      
      --genome              hg19 or hg38
                                Default: hg38

      --outdir              Select which folder to write output to.
                                Default: TN_WES_results

      --samplesheet         path to case samplesheet. Can contain multiple patients/cases (each case in a separate line). 


      --server              Select which server the analysis is performed on (kga01 or lnx01)
                                Default: lnx01

      --fastq               Path to dir with fastq files
                                Default: data storage dirs at lnx01 server or kga01 server

      --fastqInput          Use Fastq as input, automatically search for relevant fastq files at KG Vejle data archive


      --skipQC              Do not run QC module


      --keepwork            keep the workfolder generated by the nextflow script.
                                Default: not set - removes the Work folder generated by nextflow

      --nomail              Does not send a mail-message when completing a script
                                Default: not set - sends mail message if the user is mmaj or raspau and only if the script has been running longer than 20 minutes.


    """.stripIndent()
}
if (params.help) exit 0, helpMessage()


def errorMessage1() {

    log.info"""

    USER INPUT ERROR: If no samplesheet is selected, the user needs to point to a folder containing relevant fastq or CRAM files... 
    Run the script with the --help parameter to see available options
    
    """.stripIndent()
}

if (!params.samplesheet && !params.fastq && !params.cram) exit 0, errorMessage1()

def errorMessage2() {

    log.info"""

    USER INPUT ERROR: Choose either fastq or CRAM as input... Not both. 
    Run the script with the --help parameter to see available options
    
    """.stripIndent()
}

if (params.fastq && params.cram) exit 0, errorMessage2()


///////////////////////////// SAMPLESHEET CHANNELS /////////////////////////////

// Samplesheet cols (fixed order)
// 0: CaseID, 1: WES.blood, 2: WES.tumor, 3: RNA tumor, 4: pcgr_tumor_code

////////////////////////////////////////////////////////////////////////////////

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[1], row[0])}
    .set { normalID_caseID }
//above: Normal sampleID (NPN), caseID

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[2], row[0])}
    .set { tumorID_caseID }

//above: tumor sampleID (NPN), caseID

channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[0], row[4])} 
    .set { caseID_pcgrID }



channel.fromPath(params.samplesheet)
    .splitCsv(sep:'\t')
    .map { row -> tuple(row[0], row[1]+"_${datapattern}")}
    .set { caseID_normalID } // use for Mutect2 --normal

///////////////// END: SAMPLESHEET CHANNELS ////////////////////////


////////////////// INPUT DATA (FASTQ) CHANNELS ///////////////////




    if (params.fastq) {
        params.reads = "${params.fastq}/**{.,-}{${datapattern}}{.,-}*R{1,2}*{fq,fastq}.gz"
    }

    if (!params.cram && !params.fastq && params.fastqInput) {
        params.reads="${dataArchive}/{lnx01,lnx02,tank_kga_external_archive}/**/*{.,-}{${datapattern}}{.,-}*R{1,2}*{fq,fastq}.gz"
    }

    if (!params.cram && (params.fastqInput ||params.fastq)) {
        channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .map { it -> [it[0], file(it[1][0]),file(it[1][1])] }
        .set { read_pairs_ch }
        // above sampleID, r1, r2
        normalID_caseID
        .join(read_pairs_ch)
        .map {tuple(it[1],it[0]+"_${datapattern}",it[2],it[3],"NORMAL")}
        .set { NN1 }
        //above: caseid, NPN_sampletype(NORMAL), normal R1, normal R2

        tumorID_caseID
        .join(read_pairs_ch)
        .map {tuple(it[1],it[0]+"_${datapattern}",it[2],it[3],"TUMOR")}
        .set { CF1 }
        //above: caseid, sampleID, ,R1, R2, type

        NN1.concat(CF1)
        .set { case_fastq_input_ch }
        //above: NN2 and CF2 in the same channel (same structure as NN2 and CF2)
    }

    ////////////////// INPUT DATA (CRAM) CHANNELS ///////////////////

    if (!params.cram && !params.fastqInput && !params.fastq) {
        cramfiles="${dataArchive}/{lnx01,lnx02,tank_kga_external_archive}/**/*{_,-}{${datapattern}}*.cram"
        craifiles="${dataArchive}/{lnx01,lnx02,tank_kga_external_archive}/**/*{_,-}{${datapattern}}*.crai"
    }

    if (params.cram ) {
        cramfiles="${params.cram}/*{_,-}{${datapattern}}*.cram"
        craifiles="${params.cram}/*{_,-}{${datapattern}}*.crai"
    }

    if (!params.fastqInput && !params.fastq) {
        Channel
        .fromPath(cramfiles)
        .map { tuple(it.baseName.tokenize('_').get(0),it) }
        .set { sampleID_cram }
        // above: sampleID, sampleCRAM
        Channel
        .fromPath(craifiles)
        .map { tuple(it.baseName.tokenize('_').get(0),it) }
        .set { sampleID_crai }
        // above: sampleID, sampleCRAI

        // Join with samplesheet:
        normalID_caseID // sampleID normal, caseID
        .join(sampleID_cram).join(sampleID_crai)
        .map {tuple(it[1],it[0]+"_${datapattern}", it[2],it[3],"NORMAL")}
        .set { cram_normal }
        //above structure: caseID, NPN_EV8, CRAM, CRAI, NORMAL
        

        tumorID_caseID
        .join(sampleID_cram).join(sampleID_crai)
        .map {tuple(it[1],it[0]+"_${datapattern}",it[2],it[3],"TUMOR")}
        .set { cram_tumor }
        //above structure: caseID, NPN_EV8, CRAM, CRAI, TUMOR
        
        cram_normal.concat(cram_tumor)
        .set { case_npn_cram_crai_ch }
        // caseID, NPN, CRAM, CRAI

        case_npn_cram_crai_ch
        .filter{it =~ /NORMAL/}
        .set { normal_cram_ch }

        case_npn_cram_crai_ch 
        .filter{it =~ /TUMOR/}
        .set { tumor_cram_ch }
        
        normal_cram_ch
        .join(tumor_cram_ch)
        .set { tumorNormal_cram_ch }

        normal_cram_ch
        .concat(tumor_cram_ch)
        .set { cram_per_sample_ch }

    }

log.info """\

========================================================
KGA Vejle paired tumor-normalWES nextflow pipeline v3
========================================================
results     : $params.outdir
user        : $user
rundir      : $params.rundir
runtype     : $runtype
runID       : $date.$user
"""

include { 

         inputFiles_symlinks_cram;
         tb_haplotypecaller;
         pcgr_v212_mutect2;
         SUB_DNA_PREPROCESS;
         SUB_DNA_QC;
         SUB_PAIRED_TN } from "./modules/tumorBoard.modules.v1.nf" 


workflow {
    if (params.fastqInput || params.fastq) {

        SUB_DNA_PREPROCESS(case_fastq_input_ch)

        SUB_DNA_PREPROCESS.out.finalAln
        .filter{it =~ /NORMAL/}  
        .set {normal_cram_ch }

        SUB_DNA_PREPROCESS.out.finalAln
        .filter{it =~ /TUMOR/}  
        .set {tumor_cram_ch }

        normal_cram_ch.join(tumor_cram_ch)
        .set { tumorNormal_cram_ch }

        normal_cram_ch.concat(tumor_cram_ch)
        .set {cram_per_sample_ch}
    }
    
    if (!params.fastqInput && !params.fastq) {
        inputFiles_symlinks_cram(cram_per_sample_ch)
    }

    if (!params.skipQC) {
        SUB_DNA_QC(cram_per_sample_ch)
    }
    tb_haplotypecaller(normal_cram_ch)
    SUB_PAIRED_TN(tumorNormal_cram_ch, caseID_pcgrID)
}


/*
    if (!params.fastqInput && !params.fastq) {
        inputFiles_symlinks_cram(cram_per_sample_ch)
        tb_haplotypecaller(case_npn_cram_crai_ch)  // caseid, npn, cram, crai, type

        

        tb_cram_bam.out.bam
        .filter{it =~ /NORMAL/}
        .set { bam_normals_ch }
        // above: caseid, npn, cram, crai, type
        
        tb_cram_bam.out.bam
        .filter{it =~ /TUMOR/}
        .set { bam_tumor_ch }

        bam_normals_ch
        .join(bam_tumor_ch)
        .set { tumorNormal_cram_ch }
      // above structure: tuple val(caseID), val(sampleID_normal), path(bamN), path(baiN),val(typeN), val(sampleID_tumor),path(bamT), path(baiT),val(typeT)
    
        
    }
*/

        



/*
workflow.onComplete {
    // Read samplesheet and determine format
    def samplesheetLines = new File(params.samplesheet).readLines()
    def numColumns = samplesheetLines[0].tokenize('\t').size()
    def germlineOnly = numColumns == 2 // Assume 'germline only' if there are only two columns
    
    // Extract the first six digits from the samplesheet name
    def samplesheetName = new File(params.samplesheet).getName()
    def samplesheetDate = samplesheetName.find(/\d{6}/)

    // Extract names from the first column of the samplesheet before the first uppercase letter
    def names = samplesheetLines.collect { line ->
        def match = (line.split('\t')[0] =~ /^[a-z]+/)[0]
        return match ? match : null
    }.findAll { it != null } // Filter out nulls, which represent no match found

    // Only send email if --nomail is not specified and duration is longer than 20 minutes
    if (!params.nomail && workflow.duration > 1200000 && workflow.success) {
        if (System.getenv("USER") in ["raspau", "mmaj"]) {
            
            def workDirMessage = params.keepwork ? "WorkDir             : ${workflow.workDir}" : "WorkDir             : Deleted"
            
            def body = """\
            Pipeline execution summary
            ---------------------------
            Pipeline completed  : Tumorboard DNA ${samplesheetDate}
            ${germlineOnly ? "Germline only" : "Paired tumor-normal"}
            Duration            : ${workflow.duration}
            Completed at        : ${workflow.complete}
            Success             : ${workflow.success}
            ${workDirMessage}
            OutputDir           : ${params.outdir ?: 'Not specified'}
            Exit status         : ${workflow.exitStatus}
            Names               : ${names.join(', ')}
            """.stripIndent()


            // Send email using the built-in sendMail function
            sendMail(to: 'Andreas.Braae.Holmgaard@rsyd.dk,Annabeth.Hogh.Petersen@rsyd.dk,Isabella.Almskou@rsyd.dk,Jesper.Graakjaer@rsyd.dk,Lene.Bjornkjaer@rsyd.dk,Martin.Sokol@rsyd.dk,Mads.Jorgensen@rsyd.dk,Rasmus.Hojrup.Pausgaard@rsyd.dk,Signe.Skou.Tofteng@rsyd.dk', subject: 'Tumorboard pipeline Update', body: body)
        }
    }

    // Handle deletion of WorkDir based on --keepwork parameter
    if (!params.keepwork && workflow.duration > 1200000 && workflow.success) {
        println("Deleting work directory: ${workflow.workDir}")
        "rm -rf ${workflow.workDir}".execute()
    }
}
*/
