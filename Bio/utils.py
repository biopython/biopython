import string
import Seq

from PropertyManager import default_manager

def translate(seq, id = None):
    if id is None:
        s = "translator"
    else:
        s = "translator.id.%d" % id
    translator = default_manager.resolve(seq.alphabet, s)
    return translator.translate(seq)

def translate_to_stop(seq, id = None):
    if id is None:
        s = "translator"
    else:
        s = "translator.id.%d" % id
    translator = default_manager.resolve(seq.alphabet, s)
    return translator.translate_to_stop(seq)

def back_translate(seq, id = None):
    if id is None:
        s = "translator"
    else:
        s = "translator.id.%d" % id
    translator = default_manager.resolve(seq.alphabet, s)
    return translator.back_translate(seq)


def transcribe(seq):
    transcriber = default_manager.resolve(seq.alphabet, "transcriber")
    return transcriber.transcribe(seq)

def back_transcribe(seq):
    transcriber = default_manager.resolve(seq.alphabet, "transcriber")
    return transcriber.back_transcribe(seq)

def ungap(seq):
    """given a sequence with gap encoding, return the ungapped sequence"""
    gap = seq.gap_char
    letters = []
    for c in seq.data:
        if c != gap:
            letters.append(c)
    return Seq.Seq(string.join(letters, ""), seq.alphabet.alphabet)

def verify_alphabet(seq):
    letters = {}
    for c in seq.alphabet.letters:
        letters[c] = 1
    try:
        for c in seq.data:
            letters[c]
    except KeyError:
        return 0
    return 1

def count_monomers(seq):
    dict = {}
    s = buffer(seq.data)  # works for strings and array.arrays
    for c in seq.alphabet.letters:
        dict[c] = string.count(s, c)
    return dict

def sum(seq, table, zero = 0.0):
    total = zero
    for c in getattr(seq, "data", seq):
        total = total + table[c]
    return total

# For ranged addition
def sum_2ple(seq, table, zero = (0.0, 0.0)):
    x, y = zero
    data = getattr(seq, "data", seq)
    for c in data:
        x2, y2 = table[c]
        x = x + x2
        y = y + y2
    return (x, y)

def total_weight(seq, weight_table = None):
    if weight_table is None:
        weight_table = default_manager.resolve(seq.alphabet, "weight_table")
    return sum(seq, weight_table)

def total_weight_range(seq, weight_table = None):
    if weight_table is None:
        weight_table = default_manager.resolve(seq.alphabet, "weight_range_table")
    return sum_2ple(seq, weight_table)
