#!/usr/bin/env python3
from reaxfit import reaxfit
#reaxfit(dump_config=True)
reax=reaxfit()
result=reax.fit()
print(result.x)
