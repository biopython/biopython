# Copyright 2009 by Osvaldo Zagordi.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Sequencing related command line application wrappers (OBSOLETE).

We have decided to remove this module in future, and instead recommend
building your command and invoking it via the subprocess module directly.
"""

from ._Novoalign import NovoalignCommandline
from ._bwa import BwaIndexCommandline, BwaAlignCommandline, BwaSamseCommandline
from ._bwa import BwaSampeCommandline, BwaBwaswCommandline, BwaMemCommandline
from ._samtools import SamtoolsViewCommandline, SamtoolsCalmdCommandline
from ._samtools import SamtoolsCatCommandline, SamtoolsFaidxCommandline
from ._samtools import SamtoolsFixmateCommandline, SamtoolsIdxstatsCommandline
from ._samtools import SamtoolsIndexCommandline, SamtoolsMergeCommandline
from ._samtools import SamtoolsMpileupCommandline, SamtoolsPhaseCommandline
from ._samtools import SamtoolsReheaderCommandline, SamtoolsRmdupCommandline
from ._samtools import (
    SamtoolsVersion0xSortCommandline,
    SamtoolsVersion1xSortCommandline,
    SamtoolsTargetcutCommandline,
)
from ._samtools import SamtoolsVersion0xSortCommandline as SamtoolsSortCommandline


# Make this explicit, then they show up in the API docs
__all__ = (
    "BwaIndexCommandline",
    "BwaAlignCommandline",
    "BwaSamseCommandline",
    "BwaSampeCommandline",
    "BwaBwaswCommandline",
    "BwaMemCommandline",
    "NovoalignCommandline",
    "SamtoolsViewCommandline",
    "SamtoolsCalmdCommandline",
    "SamtoolsCatCommandline",
    "SamtoolsFaidxCommandline",
    "SamtoolsFixmateCommandline",
    "SamtoolsIdxstatsCommandline",
    "SamtoolsIndexCommandline",
    "SamtoolsMergeCommandline",
    "SamtoolsMpileupCommandline",
    "SamtoolsPhaseCommandline",
    "SamtoolsReheaderCommandline",
    "SamtoolsRmdupCommandline",
    "SamtoolsSortCommandline",
    "SamtoolsVersion0xSortCommandline",
    "SamtoolsVersion1xSortCommandline",
    "SamtoolsTargetcutCommandline",
)
