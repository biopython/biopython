# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Main entry point for sequencing related command line application wrappers."""
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
