# -*- coding: utf-8 -*-
# Copyright 2014 by Carlos Pena.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to invoke the BOLD server over the internet.

This module provides code to work with the BOLD API provided by BOLDSYSTEMS
http://www.boldsystems.org/index.php/Resources

"""
import logging

from .api import call_id
from .api import call_taxon_search
from .api import call_taxon_data
from .api import call_specimen_data
from .api import call_sequence_data
from .api import call_full_data
from .api import call_trace_files


logging.basicConfig(format='[bold module]:%(levelname)s:%(message)s', level=logging.DEBUG)
