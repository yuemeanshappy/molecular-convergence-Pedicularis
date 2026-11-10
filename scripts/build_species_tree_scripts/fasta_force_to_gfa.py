from GetOrganelleLib.seq_parser import *
import sys
seqs = SequenceList(sys.argv[1])
with open(sys.argv[2], "w") as output_h, open(sys.argv[2] + ".name_map", "w") as name_h:
    counter = 0
    for record in seqs:
        counter += 1
        output_h.write("S\t{}\t{}\tLN:i:{}\tRC:i:{}\n".format(counter, record.seq, len(record.seq), len(record.seq) * 100))
        name_h.write("{}\t{}\n".format(record.label, counter))

