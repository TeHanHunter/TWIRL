# activate_qlp_env.sh
#
# Sets up the QLP 0.13.2 environment for qlp lctools detrend and qlp lctools hlsp.
#
# Usage (source into current shell):
#   source /pdo/users/tehan/TWIRL/scripts/activate_qlp_env.sh
#
# Then run qlp commands via:
#   qlp_run lctools detrend --help
#   qlp_run lctools hlsp --help
#
# Notes:
# - PYTHONPATH is set to the TWIRL TGLC fork only; FFITools paths are intentionally excluded
#   to prevent qlp 0.1 (FFITools) from shadowing pip-installed qlp 0.13.2
# - LD_LIBRARY_PATH includes anaconda2 (libffi.so.6) and python-3.11.9
# - For Sector < 67, pass --flag-type spoc --flag-source fits to hlsp
#   SPOC flag files are at /pdo/qlp-data/spocflags/ and confirmed present for all Sector 56 CCDs
# - All output paths must stay within /pdo/users/tehan/

export TWIRL_QLP_PYTHON=/pdo/app/qlp-environment/.venv/bin/python
export LD_LIBRARY_PATH=/pdo/app/anaconda/anaconda2-4.4.0/lib:/pdo/app/python-versions/python-3.11.9/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH=/pdo/users/tehan/tess-gaia-light-curve-twirl

qlp_run() {
    "$TWIRL_QLP_PYTHON" -m qlp "$@"
}
export -f qlp_run

echo "QLP env active: qlp $("$TWIRL_QLP_PYTHON" -c 'import qlp; print(qlp.__version__)') | $("$TWIRL_QLP_PYTHON" --version)"
