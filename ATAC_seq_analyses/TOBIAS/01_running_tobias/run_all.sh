for i in *yaml; do
    snakemake --configfile ${i} --cores 50 --use-conda --printshellcmds --until bindetect
done