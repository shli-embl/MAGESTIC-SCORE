wget http://steinmetzlab.embl.de/shiny/SCORE/s288c_chromatin_annotations.tar.gz
tar xvfz s288c_chromatin_annotations.tar.gz
rm s288c_chromatin_annotations.tar.gz
if [ ! -d inputs/data ]; then
	mkdir inputs/data
fi
cp -rf annotations/* inputs/data/
rm -r annotations