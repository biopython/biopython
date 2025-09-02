# Type stub for _meltingtemp_exact C extension
from typing import Optional

def tm_nn_exact(
    seq: str,
    dnac1: float = 25.0,
    dnac2: float = 25.0,
    selfcomp: bool = False,
    Na: float = 50.0,
    K: float = 0.0,
    Tris: float = 0.0,
    Mg: float = 0.0,
    dNTPs: float = 0.0,
    saltcorr: int = 5,
) -> float: ...