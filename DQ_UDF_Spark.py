# Databricks notebook source
# MAGIC %run "./MTNAP bealign and tn93 defs"

# COMMAND ----------

import pyspark.sql.functions as F
import pyspark.sql.types as T
from Bio.Seq import Seq
import re

# set up intial seq_dict dictionary
# test bealign2 and df transform
# DHAP/HIV-Trace/ForDQScript/test/part-00000-tid-3928563708143819176-c99c4a64-7908-4a36-8105-510a9c2f34dd-812-1-c000.txt

iupac_code = set(["A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", "-", "*"])
ambig_code = set(["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"])


# return 1 if it contains unknown chars, 0 otherwise
def checkUnknownChars(seq):
    if set(seq).issubset(iupac_code):
        returnCode = 0
    else:
        returnCode = 1

    return returnCode


def checkAmbiguities(seq):
    ambig_count = 0
    for nt in seq:
        if nt in ambig_code:
            ambig_count += 1
    if float(ambig_count) / len(seq) >= 0.05:
        ambig5percent = 1
    else:
        ambig5percent = 0

    return ambig5percent


def checkNString(seq):
    if "NNNNN" in seq:
        nString = 1
    else:
        nString = 0

    return nString


def removeInvalidChars(seq):
    newSeq = seq
    for char in set(seq).difference(iupac_code):
        # print(seq_dict[seq_id]['seq'])
        newSeq = re.sub(re.escape(char), "", seq)
    return newSeq


# seqId - (String) sequence ID
# seq - (String) sequence
def describeSeqs_udf(seqId, seq):
    # print(seqId)
    seq = seq.upper()
    unkChars = checkUnknownChars(seq)
    seq = removeInvalidChars(seq)
    skip = 0  ##'SeqID\tUnknownCharacters\tNString\tAmbig5%\tPRRT_HXB2_Link\tPR_Forward\tPR_Reverse\tRT_Forward\tRT_Reverse\tINT_Forward\tINT_Reverse\tRTPR_Swap\n')
    # return a Sequence of values that represent the output.
    ambigs = checkAmbiguities(seq)
    nString = checkNString(seq)
    pr_Uni = 0
    pr_Bi = 0
    rt_Uni = 0
    rt_Bi = 0
    int_Uni = 0
    int_Bi = 0
    RTPR = 0
    PRRT_HXB2_Link = 0

    if len(seq) < 100 or unkChars == 1:
        # make everything 0
        return (seqId, unkChars, nString, ambigs, 0, 0, 0, 0, 0, 0, 0, 0)

    # check PR, RT, IN presence and direction
    regions = ["pr", "rt", "int"]
    # TODO use custom RT region reference.
    for region in regions:
        if skip != 1:
            # print ('PERFORMING ALIGNMENT AND TN93 DISTANCE CALCULATIONS TO HXB2 REFERENCE FOR '+region+' REGION\n')

            args_bealign = "-r HXB2_" + region + " -m HIV_BETWEEN_F"
            args_bealign_R = "-r HXB2_" + region + " -m HIV_BETWEEN_F -R"
            args = bealignParser.parse_args(args_bealign.split(" "))
            args_R = bealignParser.parse_args(args_bealign_R.split(" "))

            if region == "rt":  # use custom reference
                refString = "CCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACA"
                args.reference = refString
                args_R.reference = refString
            else:
                refString = references["HXB2_" + region].load().seq

            # print(ref)
            bealign_Uni = bealign_MTNAP(
                seqId,
                Seq(seq),
                args.reference,
                args.expected_identity,
                args.alphabet,
                args.reverse_complement,
                args.score_matrix,
                args.sort,
                args.quiet,
                args.keep_reference,
            )

            bealign_Bi = bealign_MTNAP(
                seqId,
                Seq(seq),
                args_R.reference,
                args_R.expected_identity,
                args_R.alphabet,
                args_R.reverse_complement,
                args_R.score_matrix,
                args_R.sort,
                args_R.quiet,
                args_R.keep_reference,
            )

            # configure TN93
            args_tn93 = "-a resolve -l 250"  # -t 1.0 -l 250
            args = tn93Parser.parse_args(args_tn93.split(" "))
            # print(args)
            # print('in udf!')
            match_mode = args.ambigs
            # print(match_mode)
            tn93_1 = tn93.TN93(
                verbose=args.verbose,
                ignore_gaps=args.ignore_terminal_gaps,
                json_output=args.json_output,
                max_ambig_fraction=args.max_ambig_fraction,
                minimum_overlap=args.minimum_overlap,
            )
            # if the sequence has fewer than %75 gaps, then do tn93.
            if bealign_Uni.count("-") / len(bealign_Uni) <= 0.75 and len(bealign_Uni) == len(refString):
                seqRec1 = SeqRecord(bealign_Uni, id=seqId, name=seqId, description=seqId)
                refId = "HXB2_" + region
                seqRec2 = SeqRecord(refString, id=refId, name=refId, description=refId)
                tn93_Uni = float(tn93_1.tn93_distance(seqRec1, seqRec2, match_mode).split(",")[2])
                if tn93_Uni < 0.2:
                    if region == "pr":
                        pr_Uni = 1
                    elif region == "rt":
                        rt_Uni = 1
                    elif region == "int":
                        int_Uni = 1
                    else:
                        print("unkown region " + region)
            if bealign_Bi.count("-") / len(bealign_Bi) <= 0.75 and len(bealign_Bi) == len(refString):
                seqRec1 = SeqRecord(bealign_Bi, id=seqId, name=seqId, description=seqId)
                refId = "HXB2_" + region
                seqRec2 = SeqRecord(refString, id=refId, name=refId, description=refId)

                tn93_Bi = float(tn93_1.tn93_distance(seqRec1, seqRec2, match_mode).split(",")[2])
                if tn93_Bi < 0.2:
                    if region == "pr":
                        pr_Bi = 1
                    elif region == "rt":
                        rt_Bi = 1
                    elif region == "int":
                        int_Bi = 1
                    else:
                        print("unkown region " + region)

    # CHECKING FOR INVERTED SEQUENCES IN PR/RT REGION
    HXB2_PRRT = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACA"

    args_PRRT = "-r " + HXB2_PRRT + " -m HIV_BETWEEN_F"
    args_PRRT_R = "-r " + HXB2_PRRT + " -m HIV_BETWEEN_F -R"

    args = bealignParser.parse_args(args_PRRT.split(" "))
    args_R = bealignParser.parse_args(args_PRRT_R.split(" "))

    bealign_Uni = bealign_MTNAP(
        seqId,
        Seq(seq),
        args.reference,
        args.expected_identity,
        args.alphabet,
        args.reverse_complement,
        args.score_matrix,
        args.sort,
        args.quiet,
        args.keep_reference,
    )

    bealign_Bi = bealign_MTNAP(
        seqId,
        Seq(seq),
        args_R.reference,
        args_R.expected_identity,
        args_R.alphabet,
        args_R.reverse_complement,
        args_R.score_matrix,
        args_R.sort,
        args_R.quiet,
        args_R.keep_reference,
    )

    # If the data has more than 250 gaps in the Unidirectional output, then the PR and RT are reversed.
    #    if len(prrt_dict[seq_id]) >= 300:
    #        if '-'*250 in prrt_dict[seq_id][0:300]:
    #             seq_dict[seq_id]['RTPR'] = 1

    if len(bealign_Uni) >= 300:
        if "-" * 250 in bealign_Uni[0:300]:
            RTPR = 1

    # check for HXB2 links in the PRRT region
    #'tn93 -t 0.015 -a resolve -g 1.0 -f csv -l 250
    # threshhold of .015 will be filtered after we get output.
    args = "-a resolve -g 1.0 -l 250"
    args = tn93Parser.parse_args(args.split(" "))
    match_mode = args.ambigs

    tn93_1 = tn93.TN93(
        verbose=args.verbose,
        ignore_gaps=args.ignore_terminal_gaps,
        json_output=args.json_output,
        max_ambig_fraction=args.max_ambig_fraction,
        minimum_overlap=args.minimum_overlap,
    )
    # check or PRRT links in bidirectional alignment.
    seqRec1 = SeqRecord(bealign_Bi, id=seqId, name=seqId, description=seqId)
    refId = "HXB2_PRRT"
    seqRec2 = SeqRecord(HXB2_PRRT, id=refId, name=refId, description=refId)
    tn93_temp = tn93_1.tn93_distance(seqRec1, seqRec2, match_mode)
    # print(tn93_temp)
    tn93_Uni = tn93_temp.split(",")[2]

    if tn93_Uni != "-" and (float(tn93_Uni) < 0.015):
        PRRT_HXB2_Link = 1

    PR_Forward = 1 if (pr_Uni == 1 and pr_Bi == 1) else 0
    PR_Reverse = 1 if (pr_Uni == 0 and pr_Bi == 1) else 0
    RT_Forward = 1 if (rt_Uni == 1 and rt_Bi == 1) else 0
    RT_Reverse = 1 if (rt_Uni == 0 and rt_Bi == 1) else 0
    INT_Forward = 1 if (int_Uni == 1 and int_Bi == 1) else 0
    INT_Reverse = 1 if (int_Uni == 0 and int_Bi == 1) else 0
    RTPR = 1 if (RTPR == 1 and pr_Uni == 1 and rt_Uni == 1) else 0

    # PR_Forward + PR_Reverse
    # if seq_dict[seq_id]['PR_Uni'] == 1 and seq_dict[seq_id]['PR_Bi'] == 1: outfile.write('1\t0\t')
    # elif seq_dict[seq_id]['PR_Bi'] == 1: outfile.write('0\t1\t')
    # else: outfile.write('0\t0\t')
    # RT_Forward + RT_Reverse
    # if seq_dict[seq_id]['RT_Uni'] == 1 and seq_dict[seq_id]['RT_Bi'] == 1: outfile.write('1\t0\t')
    # elif seq_dict[seq_id]['RT_Bi'] == 1: outfile.write('0\t1\t')
    # else: outfile.write('0\t0\t')
    # INT_forward + INT_Reverse
    # if seq_dict[seq_id]['INT_Uni'] == 1 and seq_dict[seq_id]['INT_Bi'] == 1: outfile.write('1\t0\t')
    # elif seq_dict[seq_id]['INT_Bi'] == 1: outfile.write('0\t1\t')
    # else: outfile.write('0\t0\t')
    # if seq_dict[seq_id]['RTPR'] == 1 and seq_dict[seq_id]['PR_Uni'] == 1 and seq_dict[seq_id]['RT_Uni'] == 1: outfile.write('1\n')
    # else: outfile.write('0\n')
    # seqID, sequence, unknown chars plus 10 dash fields.

    ##'SeqID\tUnknownCharacters\tNString\tAmbig5%\tPRRT_HXB2_Link\tPR_Forward\tPR_Reverse\tRT_Forward\tRT_Reverse\tINT_Forward\tINT_Reverse\tRTPR_Swap\n')
    return (
        seqId,
        unkChars,
        nString,
        ambigs,
        PRRT_HXB2_Link,
        PR_Forward,
        PR_Reverse,
        RT_Forward,
        RT_Reverse,
        INT_Forward,
        INT_Reverse,
        RTPR,
    )
