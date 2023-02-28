#!/usr/bin/env python3
import reaxfit
import sys
#hoge=reaxfit.reaxfit(dump_config=True)
fftemp="ffield.currentbest"
if len(sys.argv)>1:
    fn=sys.argv[1]
    print(fn)
    fftemp=fn
hoge=reaxfit.reaxfit()
hoge.config(initfile=fftemp)
pes,fns=hoge.reax()
pes-=pes[0]
pes*=0.043364124 # kcal/mol to ev
fns*=0.043364124 # kcal/mol to ev
#pes-=hoge.refE
print(pes)
print(fns)
