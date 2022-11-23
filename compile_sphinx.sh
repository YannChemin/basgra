# Actually compile docs with sphinx
sphinx-build -M html docs/source/ docs/build/ -W -a -j auto -n --keep-going

# express save
bash save.sh

# open docs in Firefox
#firefox docs/build/html/index.html

