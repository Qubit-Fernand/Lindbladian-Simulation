#!/bin/bash
pip install matplotlib
pip install qutip
pip install tqdm

initial="0"

python "exp(tL)_Euler.py" --N 2 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 3 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 4 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 5 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 6 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 7 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 8 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 9 --step 100000 --initial "$initial"
python "exp(tL)_Euler.py" --N 10 --step 100000 --initial "$initial"

# python "exp(tL)_Euler.py" --N 11 --step 100000
# python "exp(tL)_Euler.py" --N 12 --step 100000

python "exp(tL)_Euler.py" --N 2 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 3 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 4 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 5 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 6 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 7 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 8 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 9 --step 100000 --initial '1'
python "exp(tL)_Euler.py" --N 10 --step 100000 --initial '1'

python "exp(tL)_Euler.py" --N 2 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 3 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 4 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 5 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 6 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 7 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 8 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 9 --step 100000 --initial '+'
python "exp(tL)_Euler.py" --N 10 --step 100000 --initial '+'

python "exp(tL)_Euler.py" --N 2 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 3 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 4 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 5 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 6 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 7 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 8 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 9 --step 100000 --initial '01'
python "exp(tL)_Euler.py" --N 10 --step 100000 --initial '01'

python "exp(tL)_Euler.py" --N 2 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 3 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 4 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 5 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 6 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 7 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 8 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 9 --step 100000 --initial 'I'
python "exp(tL)_Euler.py" --N 10 --step 100000 --initial 'I'