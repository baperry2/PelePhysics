# Run from a CPU node on Kestrel

# Set up environment (just use default modules)
module restore

# Compile (this may take a few minutes)
make TPL && make -j

# run
./Pele3d.gnu.x86-spr.ex inputs.3d_1dArray_failedcase
