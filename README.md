# SNPs_buddy
Basic scripts to help with VCF files 

## Extract genes of interest
To extract genes of interest from vcf files there is the snakemake file `Extract_SNPs.smk`

As the output of the pipeline, there are 99 files, and finding the right gene could be painful. 

So here is a small script to help you with it.
By defautl input dir is '/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED', but you can change it in by passing `--config input_dir=MYNICEPATH` to snakemake command


In constants.py you can find path to the directory with vcf  files 

### Instruction



0. install conda or miniconda [instructions here](https://docs.anaconda.com/free/miniconda/miniconda-install/)
1. clone code to the directory of choice `git clone https://github.com/holstegelab/SNPs_buddy.git` or find it on Spider `/project/holstegelab/Share/SNPs_buddy`
    If you want to update - run the command `git pull` from the directory with the code
3. install conda environment `conda env create -f envs/snp_buddies.yaml` and activate it
4. To run this code you have to move to a directory with vcf files. (/project/holstegelab/Share/NL_VUMC_batches1_9 for example)
5. start snakemake script with command `snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --rerun-incomplete --conda-frontend conda --use-conda -c 2 --config gene=MyFancyGene`. After -c you have to mention amount of cores that will be used. I recommend to run it on nodes and use all required amount of cores. It could be an interactive node or just use `srun` for example `srun -p normal -c 32 -t 12:0:0 snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --rerun-incomplete --conda-frontend conda --use-conda -c 32 --config gene=MyFancyGene`
6. You will get a directory with a genename and some vcf-files and stat files inside


## Check samples that have a SNPs

1.  To use this script type `python Extract_samples.py`
2. After that script will ask to print the path to the vcf file
        *NOTE: there is no "autotabulation" yet, so use the command `realpath MYNICEVCF.vcf`, copy it, and paste this path as input to Python script
3. The output will be a table in the shell. The script will ask for an output file. Type a path to the output file.

## Extract by sample_name
To extract SNPs by sample name you can use the script `Sample_by_samplename.smk`
usage:
`snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Sample_by_samplename.smk --rerun-incomplete --conda-frontend conda --use-conda -c 2 --config sample_names=Sample1,Sample2,SampleN dir_name=MyFancyResults`

To change core numbers change `-c 2` to actual core numbers

use `--config sample_names` to pass all sample names that you want to insclude in output. Use comma-separated list without space.
If there are more than 1 sample it's handy to pass a directory name for output with `--config dir_name=MyFancyResluts`
If dir_name is not there output dir name will be name of 1st samplename.


snakemake --jobs 20 --cluster "sbatch -n {resources.n} --mem {resources.mem_mb} -p {resources.partition} -t {resources.time_min}" --rerun-incomplete --conda-frontend conda --use-conda --snakefile ~/Share/SNPs_buddy/Sample_by_samplename.smk




