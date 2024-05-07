from lammps import lammps
import json,os,sys,re
import numpy as np
from ._version import __version__
from .reaxfit import reaxfit, default_option
from argparse import ArgumentParser

if len(sys.argv) <2:
    print("useage: python3 -m reaxfit [init|fit|check|version]")
    sys.exit()
parser = ArgumentParser()
parser.add_argument('runtype')
parser.add_argument('others', nargs='*')
for k,v in default_option.items():
    parser.add_argument('--'+k ,type=type(v))
args=parser.parse_args()
newopts={k:v for k,v in vars(args).items() if v is not None}
del newopts["runtype"], newopts["others"]

if args.runtype == "fit":
    reax=reaxfit()
    if newopts: reax.config(**newopts)
    result=reax.fit()
    print("parameters after fitting")
    print(*[f'{p:.5f}' for p in result.x])
    print("fitted energy and force")
    print(*[f'{s:.5f}' for s in reax.E-reax.E[reax.base_indices]])
    fns=reax.F
    if reax.relative_force: fns-=fns[reax.base_indices]
    print(*[f'{f:.3f}' for f in fns])
elif args.runtype == "init":
    cfile="config.json"
    if os.path.isfile(cfile): os.replace(cfile,cfile+".bk")
    with open(cfile,"w") as f:
        print("create config.json")
        json.dump(default_option,f,indent=1)
    sys.exit()
elif args.runtype == "check":
    target_file=default_option["midfile"]
    if len(args.others) > 0: target_file = args.others[0]
    if not os.path.isfile(target_file):
        print(f"{target_file} not found")
        sys.exit(1)
    reax=reaxfit()
    reax.config(initfile=target_file,**newopts)
    pes,fns=reax.reax()
    print(f"reaxff absolute energy and force norm for {target_file}  (kcal/mol)")
    print('E:',*[f'{s:.5f}' for s in pes])
    print('F:',*[f'{s:.3f}' for s in fns])
    print(f"reaxff relative energy and force norm for {target_file}  (kcal/mol)")
    print('E:',*[f'{s:8.3f}' for s in pes-pes[reax.base_indices]])
    if reax.relative_force: fns-=fns[reax.base_indices]
    print('F:',*[f'{s:8.3f}' for s in fns])
    print(f"reference relative energy and force norm")
    print('E:',*[f'{s:8.3f}' for s in reax.refE_relative])
    print('F:',*[f'{s:8.3f}' for s in reax.refF_relative])
elif args.runtype == "version":
    if __version__ == '0.0.1':
        from subprocess import check_output
        dir_path=os.path.dirname(__file__)
        try:
            __version__=check_output(['git', 'describe', '--tags'],cwd=dir_path)
            __version__=__version__.decode('ascii').strip()
        except:
            pass
    print(__version__)
else:
    print(f"unknown runtype: {args.runtype}")
