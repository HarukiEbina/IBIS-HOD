#!/bin/bash -l


module load python
conda activate myenv
python make_many_forecasts_ELG.py 11.00 &
python make_many_forecasts_ELG.py 11.25 &
python make_many_forecasts_ELG.py 11.50 &
python make_many_forecasts_ELG.py 11.75 &
python make_many_forecasts_ELG.py 12.00 &
ps


