cd RESULTS/ITSx 

~/bin/ITSx_1.1.3/ITSx \
 -i ../all_nonchimeras.fa \
 --saveregions ITS2 \
 -o all_nochimeras \
 --cpu 8 \
 --detailed_results T \
 --preserve T
 
#get list of only Fungal ITSx sequences
cat all_nochimeras.extraction.results | awk '{if ($3 == "F") print $1;}' > Fungal_ITS_headers.txt

#subset fasta file for just Fungal ITS2 and keep only sequences longer than 100 bp
cat all_nochimeras.ITS2.fasta | seqkit grep -f Fungal_ITS_headers.txt | seqkit seq -m 100 -M 500 > all_nochimeras.Fungi.ITS2.fasta

cd ../dnabarcoder
source activate dnabarcoder-env

#search for best matches against UNITE ITS2 db
~/miniconda3/envs/dnabarcoder-env/bin/dnabarcoder/dnabarcoder.py search \
 -i ../ITSx/all_nochimeras.Fungi.ITS2.min100.fasta \
 -r ../../../../UK_SURVEY/RESULTS/dnabarcoder/unite2024ITS2.fasta -ml 50

#assign sequences to diff taxonomic groups
~/miniconda3/envs/dnabarcoder-env/bin/dnabarcoder/dnabarcoder.py classify \
 -i all_nochimeras.Fungi.ITS2.unite2024ITS2_BLAST.bestmatch \
 -c ../../../../UK_SURVEY/RESULTS/dnabarcoder/unite2024ITS2.classification \
 -cutoffs ../../../../UK_SURVEY/RESULTS/dnabarcoder/unite2024ITS2.unique.cutoffs.best.json
 
#remove size info from dynamic format files
sed 's/size.*//g' ../ITSx/all_nochimeras.Fungi.ITS2.min100.fasta > ../ITSx/all_nochimeras.Fungi.ITS2.nosize.fasta
sed 's/size[0-9A-Za-z]*//g' all_nochimeras.Fungi.ITS2.unite2024ITS2_BLAST.classification.dynamic.format > all_nochimeras.Fungi.ITS2.unite2024ITS2_BLAST.classification.dynamic.format.nosize


#RDP to get confidence scores as TBAS doesnt give them
java -jar ~/bin/rdp_classifier_2.14/dist/classifier.jar classify \
 --gene fungalits_unite \
-o all_nochimeras.Fungi.ITS2.RDP.txt \
../ITSx/all_nochimeras.Fungi.ITS2.fasta

#dynamic clustering based on dnabarcoer cutoffs
#run RStudio from conda environment
conda activate dynamic_clustering
/Applications/RStudio.app/Contents/MacOS/RStudio
