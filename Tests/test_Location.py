#!/usr/bin/env python
"""Test the Location code located in SeqFeature.py

This checks to be sure fuzzy and non-fuzzy representations of locations
are working properly.
"""
from Bio import SeqFeature

# --- test fuzzy representations
print "Testing fuzzy representations..."

# check the positions alone
exact_pos = SeqFeature.ExactPosition(5)
within_pos = SeqFeature.WithinPosition(10, 3)
between_pos = SeqFeature.BetweenPosition(20, 4)
before_pos = SeqFeature.BeforePosition(15)
after_pos = SeqFeature.AfterPosition(40)

print "Exact:", exact_pos
print "Within:", within_pos
print "Between:", between_pos
print "Before:", before_pos
print "After:", after_pos

# put these into Locations
location1 = SeqFeature.FeatureLocation(exact_pos, within_pos)
location2 = SeqFeature.FeatureLocation(before_pos, between_pos)
location3 = SeqFeature.FeatureLocation(within_pos, after_pos)

for location in [location1, location2, location3]:
    print "Location:", location
    print "   Start:", location.start
    print "   End  :", location.end

# --- test non-fuzzy represenations
print "Testing non-fuzzy representations..."
for location in [location1, location2, location3]:
    print "Location:", location
    print "  Non-Fuzzy Start:", location.nofuzzy_start
    print "  Non-Fuzzy End:", location.nofuzzy_end

