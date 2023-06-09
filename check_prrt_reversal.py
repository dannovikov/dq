from Bealign import bealign, references
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tn93 import tn93 as TN93
from emit import emit
import json
import sys


def check_prrt_reversal(seq_dict):
    tn93 = TN93.TN93(max_ambig_fraction=1.0, minimum_overlap=0)
    HXB2_PRRT = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACA"
    ref = SeqRecord(Seq(HXB2_PRRT), id="ref")
    # seq = SeqRecord(Seq(seq_dict["clean_seq"]), id="seq")
    seq = seq_dict["clean_seq"]
    RTPR = 0
    PRRT_HXB2_LINK = 0

    aln_uni = bealign(seq, ref, do_reverse_complement=False)

    if len(aln_uni) >= 300:
        if "-" * 250 in aln_uni[0:300]:
            RTPR = 1

    aln_bi = SeqRecord(bealign(seq, ref, do_reverse_complement=True), id="seq_bi")
    tn93_output = tn93.tn93_distance(aln_bi, ref, "resolve")

    distance = tn93_output.split(",")[2]
    if distance != "-" and float(distance) < 0.015:
        PRRT_HXB2_LINK = 1

    seq_dict["flags"]["pr_forward"] = 1 if (seq_dict["flags"]["pr_uni"] == 1 and seq_dict["flags"]["pr_bi"] == 1) else 0
    seq_dict["flags"]["pr_reverse"] = 1 if (seq_dict["flags"]["pr_uni"] == 0 and seq_dict["flags"]["pr_bi"] == 1) else 0
    seq_dict["flags"]["rt_forward"] = 1 if (seq_dict["flags"]["rt_uni"] == 1 and seq_dict["flags"]["rt_bi"] == 1) else 0
    seq_dict["flags"]["rt_reverse"] = 1 if (seq_dict["flags"]["rt_uni"] == 0 and seq_dict["flags"]["rt_bi"] == 1) else 0
    seq_dict["flags"]["int_forward"] = 1 if (seq_dict["flags"]["int_uni"] == 1 and seq_dict["flags"]["int_bi"] == 1) else 0
    seq_dict["flags"]["int_reverse"] = 1 if (seq_dict["flags"]["int_uni"] == 0 and seq_dict["flags"]["int_bi"] == 1) else 0
    seq_dict["flags"]["RTPR"] = 1 if (RTPR == 1 and seq_dict["flags"]["pr_uni"] == 1 and seq_dict["flags"]["rt_uni"] == 1) else 0

    return seq_dict


if __name__ == "__main__":
    # stdin will contain the sequence to be aligned and other metadata, all in json format
    in_text = sys.stdin.read().strip()
    seq_dict = json.loads(in_text)
    if isinstance(seq_dict, str):
        seq_dict = json.loads(seq_dict)

    seq_dict = check_prrt_reversal(seq_dict)
    emit(json.dumps(seq_dict))
