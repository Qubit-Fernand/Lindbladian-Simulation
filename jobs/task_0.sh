pip install matplotlib
pip install qutip
pip install tqdm

# N from 2 to 10, r in 1,2,3,7
# For Trotterization, the e^{H} is exact, 
# only the dissipative e^{D} is approximated by 1000 steps.
initial="I"

python "exp(tL)_Euler.py" --N 2 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 2 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 2 --r 2 --initial "$initial" 
python Dissipative_TFIM_rho.py --N 2 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 2 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 3 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 3 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 3 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 3 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 3 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 4 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 4 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 4 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 4 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 4 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 5 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 5 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 5 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 5 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 5 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 6 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 6 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 6 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 6 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 6 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 7 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 7 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 7 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 7 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 7 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 8 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 8 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 8 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 8 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 8 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 9 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 9 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 9 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 9 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 9 --r 7 --initial "$initial"

python "exp(tL)_Euler.py" --N 10 --step 100000 --initial "$initial"
python Dissipative_TFIM_rho.py --N 10 --r 1 --initial "$initial"
python Dissipative_TFIM_rho.py --N 10 --r 2 --initial "$initial"
python Dissipative_TFIM_rho.py --N 10 --r 3 --initial "$initial"
python Dissipative_TFIM_rho.py --N 10 --r 7 --initial "$initial"

