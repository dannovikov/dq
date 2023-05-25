import json
import sys
import os
from tn93 import tn93
from Bio.Seq import Seq, SeqRecord
from Bealign import bealign, references
from emit import emit

if __name__ == "__main__":
    # stdin will contain the sequence to be aligned and other metadata, all in json format
    in_text = sys.stdin.read().strip()

    seq_dict = json.loads(in_text)
    seq = seq_dict["clean_seq"]

    seq_dict["flags"]["pr_uni"] = 0
    seq_dict["flags"]["rt_uni"] = 0
    seq_dict["flags"]["int_uni"] = 0
    seq_dict["flags"]["pr_bi"] = 0
    seq_dict["flags"]["rt_bi"] = 0
    seq_dict["flags"]["int_bi"] = 0

    if seq_dict["flags"]["unknownChars"] == 1:
        print(json.dumps(seq_dict))
        sys.exit(0)

    tn93 = tn93.TN93(max_ambig_fraction=1.0, minimum_overlap=0)
    regions = ["HXB2_pr", "HXB2_rt", "HXB2_int"]
    for region in regions:
        if region == "HXB2_rt":
            ref = "CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACA"
            ref = SeqRecord(Seq(ref), id="HXB2_rt")
        else:
            ref = references[region].load()

        aln_uni = bealign(seq, ref, do_reverse_complement=False)
        aln_bi = bealign(seq, ref, do_reverse_complement=True)
        seq_uni = SeqRecord(Seq(aln_uni), id="seq_uni")
        seq_bi = SeqRecord(Seq(aln_bi), id="seq_bi")
        for i, aln in enumerate([seq_uni, seq_bi]):
            if aln.count("-") / len(aln) <= 1 and len(aln) == len(ref):
                tn93_output = tn93.tn93_distance(aln, ref, "resolve")
                distance = float(tn93_output.split(",")[2])
                if distance < 0.2:
                    if region == "HXB2_pr":
                        if i == 0:
                            seq_dict["flags"]["pr_uni"] = 1
                        else:
                            seq_dict["flags"]["pr_bi"] = 1
                    elif region == "HXB2_rt":
                        if i == 0:
                            seq_dict["flags"]["rt_uni"] = 1
                        else:
                            seq_dict["flags"]["rt_bi"] = 1
                    elif region == "HXB2_int":
                        if i == 0:
                            seq_dict["flags"]["int_uni"] = 1
                        else:
                            seq_dict["flags"]["int_bi"] = 1
    emit(json.dumps(seq_dict))
