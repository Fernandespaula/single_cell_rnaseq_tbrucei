
# Create directory for EM-multimapping matrix filtering with starsolo cellfiltering parameter:

cd <PATH TO ..Gene/raw> 

mkdir raw_em
mkdir filtered_em

cp features.tsv raw_em
cp barcodes.tsv raw_em
cp UniqueAndMult-EM.mtx raw_em

mv  raw_em/UniqueAndMult-EM.mtx raw_em/matrix.mtx

echo "Done."
