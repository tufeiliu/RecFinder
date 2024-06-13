import multiprocessing as mp
import os
import sys

# Extract the version number defined in pyproject.toml
try:
    from importlib.metadata import version, PackageNotFoundError
    __version__ = version(__name__)
except PackageNotFoundError:
    from ._version import __version__


residues = ['C', 'S', 'T', 'A', 'G', 'P', 'D', 'E', 'Q', 'N', 'H', 'R', 'K', 'M', 'I', 'L', 'V', 'F', 'Y', 'W']

# Find the total number of threads on the host machine
nThreadsDefault = mp.cpu_count()

os.environ["OPENBLAS_NUM_THREADS"] = "1"    # fix issue with numpy/openblas. Will mean that single threaded options aren't automatically parallelised 

my_env = os.environ.copy()

# Check if running inside a Conda environment
if 'CONDA_PREFIX' in my_env or 'CONDA_DEFAULT_ENV' in my_env:
    # Modify my_env to use the Conda environment
    conda_prefix = my_env.get('CONDA_PREFIX')
    if conda_prefix:
        conda_bin = os.path.join(conda_prefix, 'bin')
        my_env['PATH'] = conda_bin + os.pathsep + my_env['PATH']
else:
    if getattr(sys, 'frozen', False):
        if 'LD_LIBRARY_PATH_ORIG' in my_env:
            my_env['LD_LIBRARY_PATH'] = my_env['LD_LIBRARY_PATH_ORIG']
        else:
            my_env['LD_LIBRARY_PATH'] = ''
        if 'DYLD_LIBRARY_PATH_ORIG' in my_env:
            my_env['DYLD_LIBRARY_PATH'] = my_env['DYLD_LIBRARY_PATH_ORIG']
        else:
            my_env['DYLD_LIBRARY_PATH'] = ''