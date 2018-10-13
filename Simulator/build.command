#!/bin/sh
cd -- "$(dirname "$BASH_SOURCE")"
#python3 setup.py build
echo "from simulator import *" > build/startup.py
export PYTHONPATH=$(find build/lib* | head -n 1)
cp $PYTHONPATH/simulator.* ../SimulatorPy/
export PYTHONSTARTUP=build/startup.py
python3
