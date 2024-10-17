"""
A module to interact with BinaryCIF-formatted files.
"""

import gzip
from collections import deque
from typing import Optional

import numpy as np

try:
    import msgpack
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install msgpack to use Bio.PDB.binaryCIF (e.g. pip install msgpack)"
    ) from None

import Bio.PDB._bcif_helper as _bcif_helper
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder

# https://github.com/ihmwg/python-ihm/blob/main/ihm/format_bcif.py
# https://numpy.org/doc/stable/reference/arrays.dtypes.html#
# The "<" tells NumPy to use little endian representation.
# BinaryCIF always uses little endian.
_dtypes = {
    1: np.dtype("<i1"),  # Int8
    2: np.dtype("<i2"),  # Int16
    3: np.dtype("<i4"),  # Int32
    4: np.dtype("<u1"),  # UInt8
    5: np.dtype("<u2"),  # UInt16
    6: np.dtype("<u4"),  # UInt32
    32: np.dtype("<f4"),  # Float32
    33: np.dtype("<f8"),  # Float64
}


def _byte_array_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "ByteArray"

    dtype = _dtypes[encoding["type"]]
    column["data"]["data"] = np.frombuffer(column["data"]["data"], dtype)
    column["data"]["encoding"].pop()


def _fixed_point_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "FixedPoint"

    dtype = _dtypes[encoding["srcType"]]
    factor = encoding["factor"]
    data = column["data"]["data"]
    assert data.dtype.type in (np.int32, np.uint32)
    decoded_data = np.divide(data, factor, dtype=dtype)

    column["data"]["data"] = decoded_data
    column["data"]["encoding"].pop()


def _interval_quantization_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "IntervalQuantization"

    min_val = encoding["min"]
    max_val = encoding["max"]
    num_steps = encoding["num_steps"]
    delta = max_val - min_val / (num_steps - 1)
    data = column["data"]["data"]
    dtype = _dtypes[encoding["srcType"]]
    decoded_data = np.add(min_val, np.multiply(data, delta, dtype=dtype), dtype=dtype)

    column["data"]["data"] = decoded_data
    column["data"]["encoding"].pop()


def _run_length_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "RunLength"

    data = column["data"]["data"]
    dtype = _dtypes[encoding["srcType"]]
    decoded_data = np.repeat(data[::2].astype(dtype), data[1::2])

    assert len(decoded_data) == encoding["srcSize"]
    column["data"]["data"] = decoded_data
    column["data"]["encoding"].pop()


def _delta_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "Delta"

    dtype = _dtypes[encoding["srcType"]]
    data = column["data"]["data"]
    decoded_data = data.astype(dtype, copy=False)
    decoded_data[0] += encoding["origin"]
    decoded_data.cumsum(out=decoded_data)

    column["data"]["data"] = decoded_data
    column["data"]["encoding"].pop()


def _integer_packing_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "IntegerPacking"

    byte_count = encoding["byteCount"]
    src_size = encoding["srcSize"]
    is_unsigned = encoding["isUnsigned"]

    if is_unsigned:
        dtype = np.dtype("<u4")
    else:
        dtype = np.dtype("<i4")

    data = column["data"]["data"]
    assert byte_count == data.dtype.itemsize
    assert np.issubdtype(data.dtype, np.unsignedinteger) == is_unsigned
    decoded_data = np.empty((src_size,), dtype)
    _bcif_helper.integer_unpack(data, decoded_data)

    column["data"]["data"] = decoded_data
    column["data"]["encoding"].pop()


def _string_array_decoder(column):
    encoding = column["data"]["encoding"][-1]
    assert encoding["kind"] == "StringArray"

    offsets_column = {
        "data": {
            "data": encoding["offsets"],
            "encoding": encoding["offsetEncoding"],
        }
    }
    lookup_column = {
        "data": {
            "data": column["data"]["data"],
            "encoding": encoding["dataEncoding"],
        }
    }

    string_data = encoding["stringData"]
    offsets = _decode(offsets_column)
    unique_strings = np.empty((len(offsets) - 1,), dtype=object)

    for index in range(len(unique_strings)):
        unique_string = string_data[offsets[index] : offsets[index + 1]]
        unique_strings[index] = unique_string

    lookups = _decode(lookup_column)
    column["data"]["data"] = unique_strings[lookups]
    column["data"]["encoding"].pop()


_decoders = {
    "ByteArray": _byte_array_decoder,
    "FixedPoint": _fixed_point_decoder,
    "IntervalQuantization": _interval_quantization_decoder,
    "RunLength": _run_length_decoder,
    "Delta": _delta_decoder,
    "IntegerPacking": _integer_packing_decoder,
    "StringArray": _string_array_decoder,
}


def _decode(column):
    # Note that decode modifies the column.
    encodings = deque(column["data"]["encoding"])
    column["data"]["encoding"] = encodings

    while encodings:
        encoding = encodings[-1]
        _decoders[encoding["kind"]](column)

    return column["data"]["data"]


class BinaryCIFParser:
    """A parser for BinaryCIF files.

    See the `BinaryCIF specification <https://github.com/molstar/BinaryCIF>`_.
    """

    def __init__(self):
        """Initialize a BinaryCIF parser."""
        self._structure_builder = StructureBuilder()

    def _get_hetero_field(self, atom_group: str, component_id: str) -> str:
        if atom_group == "HETATM":
            hetero_field = "W" if component_id in ("HOH", "WAT") else "H"
        else:
            hetero_field = " "

        return hetero_field

    def _get_residue_ids(self, columns):
        atom_groups = _decode(columns["_atom_site.group_PDB"])
        component_ids = _decode(columns["_atom_site.label_comp_id"])
        hetero_fields = [
            self._get_hetero_field(atom_group, component_id)
            for atom_group, component_id in zip(atom_groups, component_ids)
        ]
        insertion_codes = [
            code or " " for code in _decode(columns["_atom_site.pdbx_PDB_ins_code"])
        ]
        sequence_ids = _decode(columns["_atom_site.auth_seq_id"])

        return list(zip(hetero_fields, sequence_ids, insertion_codes))

    def _get_atoms(self, columns):
        names = _decode(columns["_atom_site.label_atom_id"])
        x_list = _decode(columns["_atom_site.Cartn_x"])
        y_list = _decode(columns["_atom_site.Cartn_y"])
        z_list = _decode(columns["_atom_site.Cartn_z"])
        coordinates_list = np.stack((x_list, y_list, z_list), axis=1)
        b_factors = _decode(columns["_atom_site.B_iso_or_equiv"])
        occupancies = _decode(columns["_atom_site.occupancy"])
        alt_ids = [
            str(alt_id or " ") for alt_id in _decode(columns["_atom_site.label_alt_id"])
        ]
        serial_numbers = _decode(columns["_atom_site.id"])
        type_symbols = _decode(columns["_atom_site.type_symbol"])

        return [
            {
                "name": names[index],
                "coord": coordinates_list[index],
                "b_factor": b_factors[index],
                "occupancy": occupancies[index],
                "altloc": alt_ids[index],
                "fullname": names[index],
                "serial_number": serial_numbers[index],
                "element": type_symbols[index],
            }
            for index in range(len(serial_numbers))
        ]

    def get_structure(self, id: Optional[str], source: str) -> Structure:
        """Parse and return the PDB structure from a BinaryCIF file.

        :param str id: the PDB code for this structure
        :param str source: the path to the BinaryCIF file
        :return: the PDB structure
        :rtype: Bio.PDB.Structure.Structure
        """
        if hasattr(source, "seek"):
            # This resets the source if source is a file handle.
            source.seek(0)

        with (
            gzip.open(source, mode="rb")
            if source.endswith(".gz")
            else open(source, mode="rb")
        ) as file:
            result = msgpack.unpack(file, use_list=True)

        columns = {
            f"{category['name']}.{column['name']}": column
            for data_block in result["dataBlocks"]
            for category in data_block["categories"]
            for column in category["columns"]
        }

        atom_model_numbers = _decode(columns["_atom_site.pdbx_PDB_model_num"])
        atom_chain_ids = _decode(columns["_atom_site.label_asym_id"])
        atom_residue_ids = self._get_residue_ids(columns)
        atom_component_ids = _decode(columns["_atom_site.label_comp_id"])
        atoms = self._get_atoms(columns)

        entry_id = _decode(columns["_entry.id"])[0]
        self._structure_builder.init_structure(id or entry_id)
        builder_model_count = 0
        builder_model_number = None
        builder_chain_id = None
        builder_residue_id = None
        builder_component_id = None

        for index in range(len(atom_model_numbers)):
            model_number = atom_model_numbers[index]
            chain_id = atom_chain_ids[index]
            residue_id = atom_residue_ids[index]
            component_id = atom_component_ids[index]

            if model_number != builder_model_number:
                self._structure_builder.init_model(builder_model_count, model_number)
                builder_model_count += 1
                builder_model_number = model_number
                builder_chain_id = None
                builder_residue_id = None
            if chain_id != builder_chain_id:
                self._structure_builder.init_chain(chain_id)
                builder_chain_id = chain_id
                builder_residue_id = None
            if residue_id != builder_residue_id or component_id != builder_component_id:
                self._structure_builder.init_residue(component_id, *residue_id)
                builder_residue_id = residue_id
                builder_component_id = component_id

            self._structure_builder.init_atom(**atoms[index])

        return self._structure_builder.get_structure()
