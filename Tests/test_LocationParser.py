#!/usr/bin/env python
"""Test the spark based location parser for parsing GenBank locations.
"""

import sys
if sys.version_info[0] >= 3:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "This deprecated module doesn't work on Python 3.")

# test the feature table parser
from Bio.GenBank import LocationParser

# a list of strings to test
test_data = (
    "467",
    "23..400",
    "join(544..589,688..1032)",
    "1..1000",
    "<345..500",
    "<1..888",
    "(102.110)",
    "(23.45)..600",
    "(122.133)..(204.221)",
    "123^124",
    "145^177",
    "join(12..78,134..202)",
    "complement(join(2691..4571,4918..5163))",
    "join(complement(4918..5163),complement(2691..4571))",
    "complement(34..(122.126))",
    # The doc example allows "J00194:(100..202)" but not the BNF
    "J00194:100..202",
    "1..1509",
    "<1..9",
    "join(10..567,789..1320)",
    "join(54..567,789..1254)",
    "10..567",
    "join(complement(<1..799),complement(5080..5120))",
    "complement(1697..2512)",
    "complement(4170..4829)",
    # added a comma from the documentation
    "join(complement(567..795),complement(21..349))",
    "join(2004..2195,3..20)",
    "<1..>336",
    "394..>402",
    "join(10000,10200..10450)",

    # a few examples from from hum1
    "join(AB001090.1:1669..1713)",
    "join(AB001090.1:1669..1713,AB001091.1:85..196)",
    "join(AB001090.1:1669..1713,AB001091.1:85..196,AB001092.1:40..248,AB001093.1:96..212,AB001094.1:71..223,AB001095.1:87..231,AB001096.1:33..211,AB001097.1:35..175,AB001098.1:213..395,AB001099.1:56..309,AB001100.1:54..196,AB001101.1:171..404,AB001102.1:160..378,210..217)",
    "join(9106..9239,9843..9993,11889..11960,16575..16650)",
    "join(<1..109,620..>674)",
    "join(AB003599.1:<61..315,AB003599.1:587..874,47..325,425..>556)",
    "join(<85..194,296..458,547..>653)",

    # brad's examples
    "5.12",
    # this next location is in the BioPerl test records, but not in GenBank.
    # I'm not sure if we really need to worry about it, but the parser
    # won't handle the ugly thing right now.
    # "J00194:(100..202),1..245,300..422"
    # one-of syntax
    "one-of(1888,1901)..2200",
    "39..one-of(1956,2843)"
    )

# test all of the data using Andrew's parser
for s in test_data:
    print "--> Trying", s
    print repr(LocationParser.parse(LocationParser.scan(s)))



