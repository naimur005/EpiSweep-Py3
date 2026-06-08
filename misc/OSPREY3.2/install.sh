#! /bin/sh
python -m pip uninstall -y osprey
python -m pip install --user osprey --find-link=wheelhouse --pre
