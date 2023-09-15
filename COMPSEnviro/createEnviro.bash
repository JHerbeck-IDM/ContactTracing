#!/bin/bash

cd ../../phloModels
python setup.py bdist_wheel --universal
cd -
cp ../../phyloModels/dist/phyloModels-0.5-py2.py3-none-any.whl enviroDeffinition
python create_img_from_dir.py calc_tree_stats/

