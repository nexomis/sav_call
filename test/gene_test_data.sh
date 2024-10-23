if [ ! -d "ref" ]; then
    mkdir -p ref/
    wget -O "ref/refseq.fa.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/375/GCF_002815375.1_ASM281537v1/GCF_002815375.1_ASM281537v1_genomic.fna.gz"
    wget -O "ref/refseq.gff.gz" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/815/375/GCF_002815375.1_ASM281537v1/GCF_002815375.1_ASM281537v1_genomic.gff.gz"
    gunzip ref/refseq.fa.gz ref/refseq.gff.gz
    samtools faidx ref/refseq.fa
fi

if [ ! -d "test_data" ]; then
    git clone https://github.com/nexomis/genvar.git
    cd genvar/batch2/
    python3 ../genvar.py ../../ref/refseq.fa config.yml 
    cd ../../

    bowtie2-build ref/refseq.fa ref/refseq
    for spl in P0 P1 P2 P3; do
        bowtie2 -p 6 -1 genvar/batch2/${spl}_R1.fq.gz -2 genvar/batch2/${spl}_R2.fq.gz -x ref/refseq -S genvar/batch2/${spl}.sam
        samtools view -h -@ 6 genvar/batch2/${spl}.sam | samtools sort -@ 6 -o genvar/batch2/${spl}.bam
        samtools index genvar/batch2/${spl}.bam
        docker run -w $PWD -v $PWD:$PWD -u $UID:$GID quay.io/biocontainers/abra2:2.24--hdcf5f25_3 abra2 --in genvar/batch2/${spl}.bam --ref ref/refseq.fa --threads 6 --index --out genvar/batch2/${spl}.abra2.bam
        rm genvar/batch2/${spl}.sam
    done
    mv genvar/batch2 test_data
    rm -rf genvar
fi

for spl in P0 P1 P2 P3; do
    ./sav_call --base-csv ${spl}.base.csv --indel-csv ${spl}.indel.csv --bam test_data/${spl}.abra2.bam --reference ref/refseq.fa --prefix-out ${spl}.fr --keep-read-extremities
    python check_results.py test_data/${spl}.indel.csv test_data/${spl}.snp.csv ${spl}.fr.fwd.csv ref/refseq.fa 1> check_forward_${spl}.out 2> check_forward_${spl}.err
done

for spl in P0 P1 P2 P3; do
    ./sav_call --base-csv ${spl}.base.csv --indel-csv ${spl}.indel.csv --bam test_data/${spl}.abra2.bam --reference ref/refseq.fa --prefix-out ${spl}.rf --R1-strand reverse --R2-strand forward
    python check_results.py test_data/${spl}.indel.csv test_data/${spl}.snp.csv ${spl}.rf.rev.csv ref/refseq.fa 1> check_reverse_${spl}.out 2> check_reverse_${spl}.err
done
