from lammps import lammps
import json,os,sys,re
import numpy as np
from ._version import __version__
from .reaxfit import reaxfit, default_option
if len(sys.argv) <2:
    print("useage: python3 -m reaxfit [init|fit|check]")
    pass
elif sys.argv[1] == "fit":
    reax=reaxfit()
    result=reax.fit()
    print("optmized parameters")
    print(*[f'{p:.5f}' for p in result.x])
    print("optmized energy and force")
    print(*[f'{s:.5f}' for s in reax.E])
    print(*[f'{f:.5f}' for f in reax.F])
elif sys.argv[1] == "init":
    cfile="config.json"
    if os.path.isfile(cfile): os.replace(cfile,cfile+".bk")
    with open(cfile,"w") as f:
        print("create config.json")
        json.dump(default_option,f,indent=1)
    sys.exit()
elif sys.argv[1] == "check":
    target_file=default_option["midfile"]
    if len(sys.argv) > 2: target_file = sys.argv[2]
    if not os.path.isfile(target_file):
        print(f"{target_file} not found")
        sys.exit(1)
    reax=reaxfit()
    reax.config(initfile=target_file)
    pes,fns=reax.reax()
    print(f"energy and force for {target_file} without optimization (kcal/mol)")
    print(*[f'{s:.5f}' for s in pes])
    print(*[f'{f:.5f}' for f in fns])

