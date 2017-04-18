# -*- coding: utf-8 -*-
# Copyright 2014 by Carlos Pena.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Code to invoke the BOLD server over the internet.

This module provides code to work with the BOLD API provided by BOLDSYSTEMS
http://www.boldsystems.org/index.php/Resources

"""
from Bio._py3k import basestring


def _prepare_sequence(seq_record):
    """Outputs a DNA sequence as string.

    Args:
        seq_record: Either sequence as string or sequence object.

    Returns:
        Sequence as string.

    """
    if isinstance(seq_record, basestring):
        return seq_record
    else:
        try:
            return str(seq_record.seq)
        except AttributeError:
            raise AttributeError("No valid sequence was found for %s." % seq_record)
