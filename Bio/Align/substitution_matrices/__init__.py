# Copyright 2019 by Michiel de Hoon.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Substitution matrices."""

import os
import string
import numpy as np

from Bio.File import as_handle


class Array(np.ndarray):
    """numpy array subclass indexed by integers and by letters."""

    def __new__(cls, alphabet=None, dims=None, data=None, dtype=float):
        """Create a new Array instance."""
        if isinstance(data, dict):
            if alphabet is not None:
                raise ValueError("alphabet should be None if data is a dict")
            if dims is not None:
                raise ValueError("dims should be None if data is a dict")
            alphabet = []
            single_letters = True
            for key in data:
                if isinstance(key, str):
                    if dims is None:
                        dims = 1
                    elif dims != 1:
                        raise ValueError("inconsistent dimensions in data")
                    alphabet.append(key)
                elif isinstance(key, tuple):
                    if dims is None:
                        dims = len(key)
                    elif dims != len(key):
                        raise ValueError("inconsistent dimensions in data")
                    if dims == 1:
                        if not isinstance(key, str):
                            raise ValueError("expected string")
                        if len(key) > 1:
                            single_letters = False
                        alphabet.append(key)
                    elif dims == 2:
                        for letter in key:
                            if not isinstance(letter, str):
                                raise ValueError("expected string")
                            if len(letter) > 1:
                                single_letters = False
                            alphabet.append(letter)
                    else:
                        raise ValueError(
                            "data array should be 1- or 2- dimensional "
                            "(found %d dimensions) in key" % dims
                        )
            alphabet = sorted(set(alphabet))
            if single_letters:
                alphabet = "".join(alphabet)
            else:
                alphabet = tuple(alphabet)
            n = len(alphabet)
            if dims == 1:
                shape = (n,)
            elif dims == 2:
                shape = (n, n)
            else:  # dims is None
                raise ValueError("data is an empty dictionary")
            obj = super().__new__(cls, shape, dtype)
            if dims == 1:
                for i, key in enumerate(alphabet):
                    obj[i] = data.get(letter, 0.0)
            elif dims == 2:
                for i1, letter1 in enumerate(alphabet):
                    for i2, letter2 in enumerate(alphabet):
                        key = (letter1, letter2)
                        value = data.get(key, 0.0)
                        obj[i1, i2] = value
            obj._alphabet = alphabet
            return obj
        if alphabet is None:
            alphabet = string.ascii_uppercase
        elif not (isinstance(alphabet, (str, tuple))):
            raise ValueError("alphabet should be a string or a tuple")
        n = len(alphabet)
        if data is None:
            if dims is None:
                dims = 1
            elif dims not in (1, 2):
                raise ValueError("dims should be 1 or 2 (found %s)" % dims)
            shape = (n,) * dims
        else:
            if dims is None:
                shape = data.shape
                dims = len(shape)
                if dims == 1:
                    pass
                elif dims == 2:
                    if shape[0] != shape[1]:
                        raise ValueError("data array is not square")
                else:
                    raise ValueError(
                        "data array should be 1- or 2- dimensional "
                        "(found %d dimensions) " % dims
                    )
            else:
                shape = (n,) * dims
                if data.shape != shape:
                    raise ValueError(
                        "data shape has inconsistent shape (expected (%s), found (%s))"
                        % (shape, data.shape)
                    )
        obj = super().__new__(cls, shape, dtype)
        if data is None:
            obj[:] = 0.0
        else:
            obj[:] = data
        obj._alphabet = alphabet
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._alphabet = getattr(obj, "_alphabet", None)

    def _convert_key(self, key):
        if isinstance(key, tuple):
            indices = []
            for index in key:
                if isinstance(index, str):
                    try:
                        index = self._alphabet.index(index)
                    except ValueError:
                        raise IndexError("'%s'" % index) from None
                indices.append(index)
            key = tuple(indices)
        elif isinstance(key, str):
            try:
                key = self._alphabet.index(key)
            except ValueError:
                raise IndexError("'%s'" % key) from None
        return key

    def __getitem__(self, key):
        key = self._convert_key(key)
        value = np.ndarray.__getitem__(self, key)
        if value.ndim == 2:
            if self.ndim == 2:
                if value.shape != self.shape:
                    raise IndexError("Requesting truncated array")
            elif self.ndim == 1:
                length = self.shape[0]
                if value.shape[0] == length and value.shape[1] == 1:
                    pass
                elif value.shape[0] == 1 and value.shape[1] == length:
                    pass
                else:
                    raise IndexError("Requesting truncated array")
        elif value.ndim == 1:
            if value.shape[0] != self.shape[0]:
                value._alphabet = self.alphabet[key]
        return value.view(Array)

    def __setitem__(self, key, value):
        key = self._convert_key(key)
        np.ndarray.__setitem__(self, key, value)

    def __contains__(self, key):
        # Follow dict definition of __contains__
        return key in self.keys()

    def __array_prepare__(self, out_arr, context=None):
        # needed for numpy older than 1.13.0
        ufunc, inputs, i = context
        alphabet = self.alphabet
        for arg in inputs:
            if isinstance(arg, Array):
                if arg.alphabet != alphabet:
                    raise ValueError("alphabets are inconsistent")
        return np.ndarray.__array_prepare__(self, out_arr, context)

    def __array_wrap__(self, out_arr, context=None):
        if len(out_arr) == 1:
            return out_arr[0]
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        args = []
        alphabet = self._alphabet
        for arg in inputs:
            if isinstance(arg, Array):
                if arg.alphabet != alphabet:
                    raise ValueError("alphabets are inconsistent")
                args.append(arg.view(np.ndarray))
            else:
                args.append(arg)

        outputs = kwargs.pop("out", None)
        if outputs:
            out_args = []
            for arg in outputs:
                if isinstance(arg, Array):
                    if arg.alphabet != alphabet:
                        raise ValueError("alphabets are inconsistent")
                    out_args.append(arg.view(np.ndarray))
                else:
                    out_args.append(arg)
            kwargs["out"] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        raw_results = super().__array_ufunc__(ufunc, method, *args, **kwargs)
        if raw_results is NotImplemented:
            return NotImplemented

        if method == "at":
            return

        if ufunc.nout == 1:
            raw_results = (raw_results,)

        results = []
        for raw_result, output in zip(raw_results, outputs):
            if raw_result.ndim == 0:
                result = raw_result
            elif output is None:
                result = np.asarray(raw_result).view(Array)
                result._alphabet = self._alphabet
            else:
                result = output
                result._alphabet = self._alphabet
            results.append(result)

        return results[0] if len(results) == 1 else results

    def __reduce__(self):
        import pickle

        values = np.array(self)
        state = pickle.dumps(values)
        alphabet = self._alphabet
        dims = len(self.shape)
        dtype = self.dtype
        arguments = (Array, alphabet, dims, None, dtype)
        return (Array.__new__, arguments, state)

    def __setstate__(self, state):
        import pickle

        self[:, :] = pickle.loads(state)

    def transpose(self, axes=None):
        """Transpose the array."""
        other = np.ndarray.transpose(self, axes)
        other._alphabet = self._alphabet
        return other

    @property
    def alphabet(self):
        """Return the alphabet property."""
        return self._alphabet

    def copy(self):
        """Create and return a copy of the array."""
        other = Array(alphabet=self._alphabet, data=self)
        return other

    def get(self, key, value=None):
        """Return the value of the key if found; return value otherwise."""
        try:
            return self[key]
        except IndexError:
            return value

    def items(self):
        """Return an iterator  of (key, value) pairs in the array."""
        dims = len(self.shape)
        if dims == 1:
            for index, key in enumerate(self._alphabet):
                value = np.ndarray.__getitem__(self, index)
                yield key, value
        elif dims == 2:
            for i1, c1 in enumerate(self._alphabet):
                for i2, c2 in enumerate(self._alphabet):
                    key = (c1, c2)
                    value = np.ndarray.__getitem__(self, (i1, i2))
                    yield key, value
        else:
            raise RuntimeError("array has unexpected shape %s" % self.shape)

    def keys(self):
        """Return a tuple with the keys associated with the array."""
        dims = len(self.shape)
        alphabet = self._alphabet
        if dims == 1:
            return tuple(alphabet)
        elif dims == 2:
            return tuple((c1, c2) for c2 in alphabet for c1 in alphabet)
        else:
            raise RuntimeError("array has unexpected shape %s" % self.shape)

    def values(self):
        """Return a tuple with the values stored in the array."""
        dims = len(self.shape)
        alphabet = self._alphabet
        if dims == 1:
            return tuple(self)
        elif dims == 2:
            n1, n2 = self.shape
            return tuple(
                np.ndarray.__getitem__(self, (i1, i2))
                for i2 in range(n2)
                for i1 in range(n1)
            )
        else:
            raise RuntimeError("array has unexpected shape %s" % self.shape)

    def update(self, E=None, **F):
        """Update the array from dict/iterable E and F."""
        if E is not None:
            try:
                alphabet = E.keys()
            except AttributeError:
                for key, value in E:
                    self[key] = value
            else:
                for key in E:
                    self[key] = E[key]
        for key in F:
            self[key] = F[key]

    def select(self, alphabet):
        """Subset the array by selecting the letters from the specified alphabet."""
        ii = []
        jj = []
        for i, key in enumerate(alphabet):
            try:
                j = self._alphabet.index(key)
            except ValueError:
                continue
            ii.append(i)
            jj.append(j)
        dims = len(self.shape)
        a = Array(alphabet, dims=dims)
        ii = np.ix_(*[ii] * dims)
        jj = np.ix_(*[jj] * dims)
        a[ii] = np.ndarray.__getitem__(self, jj)
        return a

    def _format_1D(self, fmt):
        _alphabet = self._alphabet
        n = len(_alphabet)
        words = [None] * n
        lines = []
        try:
            header = self.header
        except AttributeError:
            pass
        else:
            for line in header:
                line = "#  %s\n" % line
                lines.append(line)
        maxwidth = 0
        for i, key in enumerate(_alphabet):
            value = self[key]
            word = fmt % value
            width = len(word)
            if width > maxwidth:
                maxwidth = width
            words[i] = word
        fmt2 = " %" + str(maxwidth) + "s"
        for letter, word in zip(_alphabet, words):
            word = fmt2 % word
            line = letter + word + "\n"
            lines.append(line)
        text = "".join(lines)
        return text

    def _format_2D(self, fmt):
        alphabet = self.alphabet
        n = len(alphabet)
        words = [[None] * n for _ in range(n)]
        lines = []
        try:
            header = self.header
        except AttributeError:
            pass
        else:
            for line in header:
                line = "#  %s\n" % line
                lines.append(line)
        keywidth = max(len(c) for c in alphabet)
        keyfmt = "%" + str(keywidth) + "s"
        line = " " * keywidth
        for j, c2 in enumerate(alphabet):
            maxwidth = 0
            for i, c1 in enumerate(alphabet):
                key = (c1, c2)
                value = self[key]
                word = fmt % value
                width = len(word)
                if width > maxwidth:
                    maxwidth = width
                words[i][j] = word
            fmt2 = " %" + str(maxwidth) + "s"
            word = fmt2 % c2
            line += word
            for i, c1 in enumerate(alphabet):
                word = words[i][j]
                words[i][j] = fmt2 % word
        line = line.rstrip() + "\n"
        lines.append(line)
        for letter, row in zip(alphabet, words):
            key = keyfmt % letter
            line = key + "".join(row) + "\n"
            lines.append(line)
        text = "".join(lines)
        return text

    def __format__(self, fmt):
        return self.format(fmt)

    def format(self, fmt=""):
        """Return a string representation of the array.

        The argument ``fmt`` specifies the number format to be used.
        By default, the number format is "%i" if the array contains integer
        numbers, and "%.1f" otherwise.

        """
        if fmt == "":
            if np.issubdtype(self.dtype, np.integer):
                fmt = "%i"
            else:
                fmt = "%.1f"
        n = len(self.shape)
        if n == 1:
            return self._format_1D(fmt)
        elif n == 2:
            return self._format_2D(fmt)
        else:
            raise RuntimeError("Array has unexpected rank %d" % n)

    def __str__(self):
        return self.format()

    def __repr__(self):
        text = np.ndarray.__repr__(self)
        alphabet = self._alphabet
        if isinstance(alphabet, str):
            assert text.endswith(")")
            text = text[:-1] + ",\n         alphabet='%s')" % self._alphabet
        return text


def read(handle, dtype=float):
    """Parse the file and return an Array object."""
    with as_handle(handle) as fp:
        lines = fp.readlines()

    header = []
    for i, line in enumerate(lines):
        if not line.startswith("#"):
            break
        header.append(line[1:].strip())
    rows = [line.split() for line in lines[i:]]
    if len(rows[0]) == len(rows[1]) == 2:
        alphabet = [key for key, value in rows]
        for key in alphabet:
            if len(key) > 1:
                alphabet = tuple(alphabet)
                break
        else:
            alphabet = "".join(alphabet)
        matrix = Array(alphabet=alphabet, dims=1, dtype=dtype)
        matrix.update(rows)
    else:
        alphabet = rows.pop(0)
        for key in alphabet:
            if len(key) > 1:
                alphabet = tuple(alphabet)
                break
        else:
            alphabet = "".join(alphabet)
        matrix = Array(alphabet=alphabet, dims=2, dtype=dtype)
        for letter1, row in zip(alphabet, rows):
            assert letter1 == row.pop(0)
            for letter2, word in zip(alphabet, row):
                matrix[letter1, letter2] = float(word)
    matrix.header = header
    return matrix


def load(name=None):
    """Load and return a precalculated substitution matrix.

    >>> from Bio.Align import substitution_matrices
    >>> names = substitution_matrices.load()
    """
    path = os.path.realpath(__file__)
    directory = os.path.dirname(path)
    subdirectory = os.path.join(directory, "data")
    if name is None:
        filenames = os.listdir(subdirectory)
        try:
            filenames.remove("README.txt")
            # The README.txt file is not present in usual Biopython
            # installations, but is included in a development install.
        except ValueError:
            pass
        return sorted(filenames)
    path = os.path.join(subdirectory, name)
    matrix = read(path)
    return matrix
