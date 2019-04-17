# variant-discovery-pipeline
**CURRENTLY A WORK IN PROGRESS**

`variant-discovery-pipeline` is a pipeline for variant calling WES and WGS data. The Broad's GATK software is used
for much of the preprocessing and variant calling. Annovar annotates the identified variants. The `README.md` documents
all necessary steps to run variant calling on raw FASTQ data from quality control to variant annotation. In its current
state, `variant-discovery-pipeline` only works on the Sherlock HPC cluster at Stanford. 

1. [Quality Control](#quality-control)
2. [Preprocessing](#preprocessing)
3. [Variant Calling](#variant-calling)
4. [Annotation](#annotation)
5. [Analysis](#analysis)

## Quality Control
Before preprocessing can begin, the raw FASTQ files should be evaluated for quality control (QC).
Sequencing can introduce artefacts and other bias that can affect variant calling and downstream
analysis. It is vital we are aware of any potential defects early on. We use two tools, [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiqc](https://multiqc.info/)
to produce QC reports for FASTQ data. 

First, use `fastqc` to generate individual QC reports for each FASTQ file.
```
fastqc [input fastq]
```
evaluates a single FASTQ file. If all the FASTQ files are stored in a single directory type
```
fastqc *.fastq.gz
```
This will produce reports for each individual FASTQ file. Once these reports have been generated,
`multiqc` can summarize all the reports. Type
```
multiqc [data directory]
```
and replace `[data directory]` with the folder containing the `fastqc` reports. The file
`multiqc_report.html` will contain the summarized report info. This 
[link](https://multiqc.info/docs/#using-multiqc-reports) provides an overview of how to use
the report. 

Its up to the user how to proceed once a QC problem has been identified. Depending
on the quality of the sample, the user may decide to discard it all together. Although in the past
trimming programs were often recommending, GATK 4.0 guidelines [recommend **against** trimming as
it can hinder the BaseRecalibration preprocessing step](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_analysis_CollectBaseDistributionByCycle.php).
Instead, BaseRecalibration can usually account
for poor quality reads used by Mutect2. 

## Preprocessing
`variant-discovery-pipeline` follows a modified version of the Broad's 
[Best Practices Pipelines for Variant Discovery](https://software.broadinstitute.org/gatk/best-practices/workflow).
Not all scatter operations are performed to save time waiting on excess SLURM jobs. 

Preprocessing requires 2 items from the user. First, the user must provide an inputs file describing each sample, its FASTQ files,
Library ID, and the directory where the preprocessed files should be saved. This info should be saved in a `.csv` or `.xlsx` file. 
A reference file must also be specified that tells `variant-discovery-pipeline` where to look for items like the reference genome fasta,
`cromwell` jar, and email address to notify the user of pipeline failure. Reference genome files can be downloaded from the [Broad's
Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle) over FTP. Cromwell executables can be found on its 
[Github](https://github.com/broadinstitute/cromwell/releases).

```
Input File Columns (.csv or .xslx):
  Sample
  Type
  FASTQ1
  FASTQ2
  OutputDirectory
  Library
  SamplePrefix (Optional but still need to include column even if it'll be left blank)
  
Reference File Columns (Must be .csv formatted):
  Reference
  Path
```
Use absolute paths for the reference files. See `preprocessing_references_template.csv` for a template file with all the required
references you need to specify.
To call the preprocessing pipeline, type
```
python preprocess.py [inputs.csv] [references.csv] [log_directory]
```
`[log_directory]` is the directory to store all the cromwell logs and intermediate files. Cromwell produces a lot of intermediate
files and logs, so its best to specify this directory as somewhere with large storage volume (like `$SCRATCH`) and to delete often 
once you no longer need the logs or intermediates. 

## Variant Calling


## Annotation


## Analysis


Modified implementation of the Broad's [Best Practices Pipelines for Variant Discovery](https://software.broadinstitute.org/gatk/best-practices/workflow) for use on Stanford's Sherlock Compute Cluster. The project consists of two pipelines: PreProcessing and VariantCalling. PreProcessing performs data preprocessing on an individual sample. VariantCalling performs Mutect2 variant calling for a normal/tumor pair of samples. 

## Details
`/scripts` contains the WDL files that specify the workflows.  
`/inputs` contains the .json files that dictate which files are processed by the workflow.  
`preprocess_samples.py` preprocesses each each sample within a specified directory by running each sample through the pipeline. Given a directory with the -d parameter, the script creates a customized input .json file for each sample labeled `preprocessing_SAMPLENAME.json` that is sent to an sbatch job that controls the preprocessing workflow for that sample. All logs are stored in the cromwell-monitoring-logs directory in `$PI_SCRATCH`. 
## Workflows
`preprocessing.wdl` specifies the preprocessing tasks to create an analysis ready bam given paired end FASTQ files as input. 
## Inputs
`preprocessing.json` specifies the variables used for the `preprocessing.wdl` workflow.
## Example Commands
### Using `preprocess_samples.py`
```
python /home/groups/carilee/software/variant-discovery-pipeline/preprocess_samples.py -d /scratch/groups/carilee/cromwell-test/short-data/
```
This creates a customized input file for CTR119_short and then launches an sbatch job to control the cromwell process.  
### Running cromwell manually
With the local backend
```
java -jar /home/groups/carilee/software/cromwell-35.jar run /home/groups/carilee/software/variant-discovery-pipeline/scripts/preprocessing.wdl -i /home/groups/carilee/software/variant-discovery-pipeline/inputs/preprocessing.json
```
To use the SLURM backend, specify the config files as `your.conf`
```
java -Dconfig.file=/home/groups/carilee/software/variant-discovery-pipeline/your.conf -jar /home/groups/carilee/software/cromwell-35.jar run /home/groups/carilee/software/variant-discovery-pipeline/scripts/preprocessing.wdl -i /home/groups/carilee/software/variant-discovery-pipeline/inputs/preprocessing.json
```
