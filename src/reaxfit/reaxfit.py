from lammps import lammps
import json,os,sys,re
import numpy as np
from scipy.optimize import differential_evolution

# global parameters are required for multiprocessing
def wrap_eval(): pass
def dump_best(): pass

class reaxfit():
  def __init__(self,cfile="config.json",dump_config=False):
    self.cfile=cfile
    self.option={
            "refE_file":"refE", # reference energy
            "refF_file":"refF", # reference force norm
            "datafile":"data0", # data file with initial structure
            "initfile":"ffield.temp", # initial template file
            "midfile":"ffield.currentbest", # intermediate file
            "endfile":"ffield.end", # final file
            "bound":0.1, # define range of parameters 
            "stopfile":"STOP", # file for early stopping
            "seed":None, # file for early stopping
            "workers":4, # file for early stopping
            "maxiter":1000, # maximum number of iteration
            "scrdir":"scr" # scratch directory
    }
    isconfig=os.path.isfile(cfile)
    if dump_config:
        if isconfig: os.rename(cfile,cfile+".bk")
        with open(cfile,"w") as f:
            json.dump(self.option,f,indent=1)
        sys.exit()
    #self.config()
    return

  def config(self,**kwargs):
    global dump_best
    if hasattr(kwargs,"cfile"): self.cfile=kwargs["cfile"]
    opt=self.option
    isconfig=os.path.isfile(self.cfile)
    if isconfig:
        print(f"read {self.cfile}")
        with open(self.cfile) as f:
            jwargs=json.load(f)
            optk=opt.keys()&jwargs.keys()
            opt.update(**{k:jwargs[k] for k in optk})
            if(optk != jwargs.keys()): 
                print("unrecognized options: ",*(jwargs.keys()-optk))
    # overwrite options if given
    if kwargs:
        optk=opt.keys()&kwargs.keys()
        opt.update(**{k:kwargs[k] for k in optk})
        if(optk != kwargs.keys()): 
            print("unrecognized options: ",*(kwargs.keys()-optk))
    midfile=opt["midfile"]
    stopfile=opt["stopfile"]
    for k in opt:
        setattr(self,k,opt[k])
    # create scratch directory
    os.makedirs(self.scrdir,exist_ok=True)
    with open(self.refE_file) as f:
      refE=f.read().split()
    with open(self.refF_file) as f:
      refF=f.read().split()
    refE=np.array(refE,dtype=float)
    refF=np.array(refF,dtype=float)
    # read template and x0
    with open(self.initfile) as f:
      _template=f.read()
    regex=re.compile(r"{([\d.-]+)}?")
    x0=regex.findall(_template)
    self.x0=[float(x) for x in x0]
    # define bounds
    bounds=[]
    for x in self.x0:
      delta = abs(float(self.bound)*x)
      t = (x-delta,x+delta)
      bounds.append(t)
    self.bounds=bounds
    template=regex.sub("{}",_template)
    self.template=template
    # get elements from templates
    eles=re.findall(r"\n *(:?[A-Z][a-z]?)",template)
    self.elements=" ".join(eles)
    self.refE=refE
    self.refF=refF
    def callbackF(xk,convergence=False,fout=None):
      fname=midfile
      sys.stdout.flush()
      if fout: fname=fout
      with open(fname,"w") as f:
          _x=np.round(xk,5)
          y=[f"{{{xx}" for xx in _x]
          f.write(template.format(*y))
      if os.path.isfile(stopfile):
          print("optimization will be stopped")
          os.remove(stopfile)
          return True
    dump_best=callbackF
    return
  def set_eval(self,func):
    global wrap_eval
    def wrap_eval(*args,**kwargs):
      pid=os.getpid()
      pes,fns=self.reax(*args,pid=pid)
      return func(pes,fns,self.refE,self.refF)
    return wrap_eval
  def reax(self,xk=None,pid=0):
      x=xk if xk is not None else self.x0
      ffname=f"{self.scrdir}/ffield.reax.{pid}"
      elements=self.elements
      pes,fns=[],[]
      with open(ffname,"w") as f:
          f.write(self.template.format(*np.round(x,6)))
      lmp=lammps(cmdargs=["-log", "none", "-screen", "none",  "-nocite"])
      for i in range(1):
          cmds=["units real",
              "atom_style charge",
              f"read_data {self.datafile}",
              "pair_style      reaxff NULL",
              f"pair_coeff      * * {ffname} {elements}",
              "fix             1 all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff"]
          lmp.commands_list(cmds)
          lmp.command("run 0")
          pes.append(lmp.get_thermo("pe"))
          fns.append(lmp.get_thermo("fnorm"))
      for i in range(1,len(self.refE)):
          lmp.command(f"read_dump {i}.xyz 0 x y z box no format xyz")
          lmp.command("run 0")
          pes.append(lmp.get_thermo("pe"))
          fns.append(lmp.get_thermo("fnorm"))
      lmp.close()
      pes=np.array(pes)
      fns=np.array(fns)
      return pes,fns

  def fit(self,eval_fit=None):
      global dump_best
      if not hasattr(self,"x0"): self.config()
      #refE,refF=self.refE,self.refF
      print(f"initial {len(self.x0)} parameters : {self.x0}")
      if not eval_fit:
          @self.set_eval
          def default_eval(pes,fns,refE,refF):
            fns*=0.043364124 # kcal/mol to ev
            pes-=pes[0] # initial structure is reference
            pes*=0.043364124 # kcal/mol to ev
            pes-=refE
            fns-=refF
            return pes@pes + np.linalg.norm(fns)/len(fns)/10+fns[0]*fns[0]/10
          eval_fit=default_eval
      result = differential_evolution(eval_fit, self.bounds,workers=self.workers,x0=self.x0,updating='deferred',
                                  disp=True,maxiter=self.maxiter,tol=0.01,callback=dump_best,popsize=16,seed=self.seed)
      #print(result)
      dump_best(result.x,fout=self.endfile)
      self.result=result
      self.E,self.F=self.reax(result.x)
      return result

if __name__ =='__main__':
    reax=reaxfit()
    result=reax.fit()
    pes,fns=reax.reax(result.x)
    print(pes,fns)
