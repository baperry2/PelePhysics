# Run from a GPU node on Kestrel

# Set up environment
module purge
module load PrgEnv-gnu/8.5.0
module load cuda/12.3
module load craype-x86-milan
module list
    
# Compile (this may take a few minutes)
make TPL USE_CUDA=TRUE && make -j USE_CUDA=TRUE

# run
./Pele3d.gnu.CUDA.ex inputs.3d_1dArray_failedcase
