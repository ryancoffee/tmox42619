from typing import Dict, Tuple
from pydantic import BaseModel

class Config(BaseModel):
    chans: Dict[int,int] # port_0.hsd
    t0s: Dict[int,float] # port_0.t0
    logicthresh: Dict[int,int] # port_0.logicthresh
    vlsthresh: int
    vlswin: Tuple[int,int]
    l3offset: int
    inflate: int # inflate pads the DCT(FFT) with zeros, artificially over sampling the waveform
    expand: int # expand controls the fractional resolution for scanedges by scaling index values and then zero crossing round to intermediate integers.
