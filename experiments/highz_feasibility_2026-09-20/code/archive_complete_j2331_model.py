#!/usr/bin/env python3
"""J233156 adapter: same-gas Fe1611 opacity inside the Fe1608 spectral region.

The campaign engine remains untouched. This subclass modifies only the atomic
opacity list of the 1608 region; primary reference, native pixels, errors and
integration grid are preserved. The inherited evaluator sums opacity before
exp(-tau), Gaussian convolution and finite-pixel integration. A free shift of
this region is a local spectral-region nuisance, not an isolated 1608 frequency
or alpha measurement. The target pair remains Fe2382 minus Fe2374.
"""
import argparse
import numpy as np
import archive_complete_model as engine
CoreModel=engine.Model
class NeighborModel(CoreModel):
 def __init__(self,configuration,ncomp,free=False,oversample=9):
  super().__init__(configuration,ncomp,free,oversample)
  if self.target!='J233156-090802':raise ValueError('Target-specific overlap adapter')
  if configuration.get('opacity_lines')!={'1608':[1608,1611]}:raise ValueError('Frozen neighbor policy missing')
  neighbor=engine.old.atoms([1611])[1611]
  assert neighbor.shape==(4,5)
  matched=0
  for line in self.lines:
   if line['key']==1608:
    assert line['atom'].shape==(4,5)
    line['atom']=np.vstack([line['atom'],neighbor]);matched+=1
  assert matched==1
  self.cache=None

def main():
 parser=argparse.ArgumentParser();parser.add_argument('configuration');parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
 configuration=engine.read_json(args.configuration)
 assert configuration['neighbor_opacity_adapter_sha256']==engine.sha(__file__)
 engine.Model=NeighborModel
 engine.run_campaign(args.configuration,workers=args.workers)
if __name__=='__main__':main()
