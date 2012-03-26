import sys
import os
import time
from reportlab.lib import colors
from reportlab.lib.units import cm

from Bio.Blast.Applications import NcbitblastxCommandline
from Bio.Graphics.GenomeDiagram import Diagram, CrossLink
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

#Modify this line to point at the Artemis/ACT example data which is online at:
# https://github.com/sanger-pathogens/Artemis/tree/master/etc
input_folder = "/Applications/Artemis/Artemis.app/Contents/artemis/etc"

name = "af063097_v_b132222"
file_a = "af063097.embl"
file_b = "b132222.embl"
format_a = "embl"
format_b = "embl"
file_a_vs_b = "af063097_v_b132222.crunch"

for f in [file_a, file_b, file_a_vs_b]:
    if not os.path.isfile(os.path.join(input_folder, f)):
        print "Missing input file %s.fna" % f
        sys.exit(1)

#Only doing a_vs_b here, could also have b_vs_c and c_vs_d etc
genomes = [
    (os.path.join(input_folder, file_a), format_a),
    (os.path.join(input_folder, file_b), format_b),
]
comparisons = [os.path.join(input_folder, file_a_vs_b)]

#Create diagram with tracks, each with a feature set
assert len(genomes) >= 2 and len(genomes) == len(comparisons)+1
gd_diagram = Diagram(name, track_size=0.35, circular=False)
tracks = dict()
feature_sets = dict()
records = dict()
for f, format in genomes:
    records[f] = SeqIO.read(f, format)
    tracks[f] = gd_diagram.new_track(1, name=f, start=0, end=len(records[f]),
                                     scale_smalltick_interval=1000,
                                     scale_largetick_interval=10000,
                                     greytrack=True, greytrack_labels=0)
    feature_sets[f] = tracks[f].new_set()

print "Drawing matches..."
for i, crunch_file in enumerate(comparisons):
    q = genomes[i+1][0] #query file
    s = genomes[i][0] #subject file
    q_set = feature_sets[q]
    s_set = feature_sets[s]
    handle = open(crunch_file)
    for line in handle:
        if line[0]=="#":
            continue
        parts = line.rstrip("\n").split(None,7)
        #0 = score
        #1 = id
        #2 = S1
        #3 = E1
        #4 = seq1
        #5 = S2
        #6 = E2
        #7 = seq2
        try:
            q_start, q_end = int(parts[2]), int(parts[3])
            s_start, s_end = int(parts[5]), int(parts[6])
        except IndexError:
            sys.stderr.write(repr(line) + "\n")
            sys.stderr.write(repr(parts) + "\n")
            raise
        flip = False
        if q_start > q_end:
            flip = not flip
            q_start, q_end = q_end, q_start
        if s_start > s_end:
            flip = not flip
            s_start, s_end = s_end, s_start
        if flip:
            c = colors.Color(0, 0, 1, alpha=0.25)
            b = False
        else:
            c = colors.Color(1, 0, 0, alpha=0.25)
            b = False
        q_feature = q_set.add_feature(SeqFeature(FeatureLocation(q_start-1, q_end)),
                                                 color=c, border=b)
        s_feature = s_set.add_feature(SeqFeature(FeatureLocation(s_start-1, s_end)),
                                                 color=c, border=b)
        gd_diagram.cross_track_links.append(CrossLink(q_feature, s_feature, c, b))
        #NOTE: We are using the same colour for all the matches,
        #with transparency. This means overlayed matches will appear darker.
        #It also means the drawing order not very important.
        #Note ACT puts long hits at the back, and colours by hit score
    handle.close()

print "Drawing CDS features..."
for f, format in genomes:
    record = records[f]
    feature_set = feature_sets[f]
    #Mark the CDS features
    for cds in record.features:
        if cds.type != "CDS":
            continue
        feature_set.add_feature(cds, sigil="ARROW",
                                color=colors.lightblue,
                                border=colors.blue)

gd_diagram.draw(format="linear", fragments=3,
                orientation="landscape", pagesize=(20*cm,10*cm))
gd_diagram.write(name + ".pdf", "PDF")

gd_diagram.draw(format="circular",
                orientation="landscape", pagesize=(20*cm,20*cm))
gd_diagram.write(name + "_c.pdf", "PDF")
