python circularSponge.py $1
python cyclogeostrophic_vortex.py $1
echo "starting vortex simulation"
mpiexec -n 4 python main_vortex.py $1
mpiexec -n 4 python merge.py $1/IC/
python extract_final_vortex.py $1
python array_inception.py $1

# rm data files in future

