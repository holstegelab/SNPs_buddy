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
4. To run this code you have to move to a directory with vcf files. (`/project/holstegelab/Share/NL_VUMC_joint_calling_splitted` for example)
5. (PREFFERED) method is to use `snakemake --jobs 20 --cluster "sbatch -n {resources.n} --mem {resources.mem_mb} -p {resources.partition} -t {resources.time_min}" --rerun-incomplete --conda-frontend conda --keep-going --use-conda --snakefile project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --config gene=MyFancyGene`
In this case you'll request nodes automatically inside the script and you won't be bothered with the amount of cores

6. (ALTernative old). You can request nodes manually. You can use `salloc` for interactive nodes or just use `srun` for example `srun -p normal -c 32 -t 12:0:0 snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --rerun-incomplete --conda-frontend conda --use-conda -c 32 --config gene=MyFancyGene`
If you are using intercative node - start snakemake script with command `snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --rerun-incomplete --conda-frontend conda --use-conda -c 2 --config gene=MyFancyGene`
 After -c you have to mention amount of cores that will be used. I recommend to use all required amount of cores. 
7. You will get a directory with a genename and some vcf-files and stat files inside


## Check samples that have a SNPs
0. Script will ask to print the path to the vcf file. Since python doesn't work good with relative paths, you have to provide the full path to the file. To check the real path for the vcf file you have to use the command `realpath MYNICEVCF.vcf` and copy the output. Then paste it as input to the script

1.  To use this script type `python Extract_samples.py`

2. The output will be a table in the shell. The script will ask for an output file. Type a path to the output file.

## Extract by sample_name
To extract SNPs by sample name you can use the script `Sample_by_samplename.smk`
usage:
`snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Sample_by_samplename.smk --rerun-incomplete --conda-frontend conda --use-conda -c 2 --config sample_names=Sample1,Sample2,SampleN dir_name=MyFancyResults`

To change core numbers change `-c 2` to actual core numbers

use `--config sample_names` to pass all sample names that you want to insclude in output. Use comma-separated list without space.
If there are more than 1 sample it's handy to pass a directory name for output with `--config dir_name=MyFancyResluts`
If dir_name is not there output dir name will be name of 1st samplename.


snakemake --jobs 20 --cluster "sbatch -n {resources.n} --mem {resources.mem_mb} -p {resources.partition} -t {resources.time_min}" --rerun-incomplete --conda-frontend conda --keep-going --use-conda --snakefile ~/Share/SNPs_buddy/Sample_by_samplename.smk

## Check missigness of SNPs in your vcf
To check the missigness of SNPs in your vcf file you can use the script `low_dp_counts.smk`
0. activate conda environment `conda activate snp_buddies_env`
Usage:
1. `snakemake -c 1 --snakefile ~/Share/SNPs_buddy/low_dp_counts.smk --config vcf={/path/to/vcf} region={/path/to/region}`
`vcf` - path to the vcf file (Obligatory)
`region` - path to the region file (Optional)
As output you will get a tsv file. 1 SNP per line. 
Columns: Chromosome, position, number of samples where this SNP is not covered, percantage of samples where this SNP is not covered

