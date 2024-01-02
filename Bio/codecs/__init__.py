"""Registers a search function for codecs in this module."""


import codecs
import importlib

from typing import Union


def _search_function(encoding_name: str) -> Union[codecs.CodecInfo, None]:
    try:
        module = importlib.import_module(f"Bio.codecs.{encoding_name}")
        codec = module.Codec()

        return codecs.CodecInfo(
            name=encoding_name,
            encode=codec.encode,
            decode=codec.decode,
        )
    except ModuleNotFoundError:
        return None


codecs.register(_search_function)
