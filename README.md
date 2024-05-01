# SNPs_buddy
Basic scripts to help with VCF files 

## Extract genes of interest
To extract genes of interest from vcf files there is the snakemake file `Extract_SNPs.smk`

As the output of the pipeline, there are 99 files, and finding the right gene could be painful. 

So here is a small script to help you with it.

### Instruction

0. install conda or miniconda [instructions here](https://docs.anaconda.com/free/miniconda/miniconda-install/)
1. clone code to the directory of choice `git clone https://github.com/holstegelab/SNPs_buddy.git` or find it on Spider `/project/holstegelab/Share/SNPs_buddy`
    If you want to update - run command `git pull main` from directory with code
3. install conda environment `conda env create -f envs/snp_buddies.yaml` and activate it
4. To run this code you have to move to a directory with vcf files. (/project/holstegelab/Share/NL_VUMC_batches1_9 for example)
5. start snakemake script with command `snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --rerun-incomplete --use-conda -c 2 --config gene=MyFancyGene`. After -c you have to mention amount of cores that will be used. I recommend to run it on nodes and use all required amount of cores. It could be interactive node or just use `srun` for example `srun -p short -c 32 -t 12:0:0 snakemake --snakefile /project/holstegelab/Share/SNPs_buddy/Extract_SNPs.smk --rerun-incomplete --use-conda -c 32 --config gene=MyFancyGene`
6. You will get a directory with a genename and some vcf-files and stat files inside


