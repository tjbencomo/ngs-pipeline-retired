# variant-discovery-pipeline
CURRENTLY UNDER DEVELOPMENT AND NOT SUITED FOR PRODUCTION USE!
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
