"""Useful utilities for helping in parsing GenBank files.
"""
class FeatureValueCleaner(object):
    """Provide specialized capabilities for cleaning up values in features.

    This class is designed to provide a mechanism to clean up and process
    values in the key/value pairs of GenBank features. This is useful 
    because in cases like:
        
         /translation="MED
         YDPWNLRFQSKYKSRDA"

    you'll end up with a value with \012s and spaces in it like:
        "MED\012 YDPWEL..."

    which you probably don't want. 
    
    This cleaning needs to be done on a case by case basis since it is
    impossible to interpret whether you should be concatenating everything
    (as in translations), or combining things with spaces (as might be
    the case with /notes).
    """
    keys_to_process = ["translation"]
    def __init__(self, to_process = keys_to_process):
        """Initialize with the keys we should deal with.
        """
        self._to_process = to_process

    def clean_value(self, key_name, value):
        """Clean the specified value and return it.

        If the value is not specified to be dealt with, the original value
        will be returned.
        """
        if key_name in self._to_process:
            try:
                cleaner = getattr(self, "_clean_%s" % key_name)
                value = cleaner(value)
            except AttributeError:
                raise AssertionError("No function to clean key: %s" 
                                     % key_name)
        return value

    def _clean_translation(self, value):
        """Concatenate a translation value to one long protein string.
        """
        translation_parts = value.split()
        return "".join(translation_parts)
