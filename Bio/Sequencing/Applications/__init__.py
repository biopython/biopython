"""Main entry point for application command line wrappers related to sequencing.
"""
from _Novoalign import NovoalignCommandline
from _bwa import BwaIndexCommandline, BwaAlignCommandline, BwaSamseCommandline
from _bwa import BwaSampeCommandline, BwaBwaswCommandline
from _samtools import SamtoolsViewCommandline,SamtoolsCalmdCommandline, SamtoolsCatCommandline, SamtoolsFaidxCommandline, SamtoolsFixmateCommandline
from _samtools import SamtoolsIdxstatsCommandline, SamtoolsIndexCommandline, SamtoolsMergeCommandline, SamtoolsMpileupCommandline
from _samtools import SamtoolsPhaseCommandline, SamtoolsReheaderCommandline, SamtoolsRmdupCommandline, SamtoolsSortCommandline, SamtoolsTargetcutCommandline
#Make this explicit, then they show up in the API docs
__all__ = ["BwaIndexCommandline",
           "BwaAlignCommandline",
           "BwaSamseCommandline",
           "BwaSampeCommandline",
           "BwaBwaswCommandline",
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
           "SamtoolsTargetcutCommandline"
           ]
