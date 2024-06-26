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

from ._bwa import BwaAlignCommandline
from ._bwa import BwaBwaswCommandline
from ._bwa import BwaIndexCommandline
from ._bwa import BwaMemCommandline
from ._bwa import BwaSampeCommandline
from ._bwa import BwaSamseCommandline
from ._Novoalign import NovoalignCommandline
from ._samtools import SamtoolsCalmdCommandline
from ._samtools import SamtoolsCatCommandline
from ._samtools import SamtoolsFaidxCommandline
from ._samtools import SamtoolsFixmateCommandline
from ._samtools import SamtoolsIdxstatsCommandline
from ._samtools import SamtoolsIndexCommandline
from ._samtools import SamtoolsMergeCommandline
from ._samtools import SamtoolsMpileupCommandline
from ._samtools import SamtoolsPhaseCommandline
from ._samtools import SamtoolsReheaderCommandline
from ._samtools import SamtoolsRmdupCommandline
from ._samtools import SamtoolsTargetcutCommandline
from ._samtools import SamtoolsVersion0xSortCommandline
from ._samtools import SamtoolsVersion0xSortCommandline as SamtoolsSortCommandline
from ._samtools import SamtoolsVersion1xSortCommandline
from ._samtools import SamtoolsViewCommandline

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
