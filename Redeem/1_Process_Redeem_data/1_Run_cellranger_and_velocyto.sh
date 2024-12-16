# Run cellrancer
cellranger-arc count --id=young2_HSPC_Multi \
                       --reference=/home/liyr/public_data/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=library_csv.csv \
                       --localcores=25 \
                       --localmem=256
                       
cellranger-arc count --id=aged2_HSPC_Multi \
                       --reference=/home/liyr/public_data/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=aged2_library.csv \
                       --localcores=25 \
                       --localmem=256

# Run velocyto  

velocyto run10x -m /home/liyr/public_data/GRCh38_rmsk.gtf /home/liyr/Redeem/young2_HSPC_Multi /home/liyr/public_data/refdata-gex-GRCh38-2020-A/genes/genes.gtf
        
velocyto run10x -m /home/liyr/public_data/GRCh38_rmsk.gtf /home/liyr/Redeem/aged2_HSPC_Multi /home/liyr/public_data/refdata-gex-GRCh38-2020-A/genes/genes.gtf


