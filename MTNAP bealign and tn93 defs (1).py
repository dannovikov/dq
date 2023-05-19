# Databricks notebook source
# MAGIC %pip install Bio
# MAGIC %pip install BioExt
# MAGIC %pip install --upgrade tn93

# COMMAND ----------

# tn93 UDF
from tn93 import tn93
import argparse


class FASTA:
    def __init__(self, seqID, sequence):
        self.seqID = seqID
        self.sequence = sequence


def setup_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--ambigs",
        action="store",
        required=False,
        help="""
        How to handle ambiguities. This can be one of four options:
        average - Averages the possible nucleotide values for each ambiguity in a sequence;
        resolve - Tries to resolve ambiguities;
        skip - Ignores gaps and ambiguities;
        gapmm - Treats gaps in only one sequence as 'N's;
        """,
    )
    parser.add_argument(
        "-g",
        "--max_ambig_fraction",
        action="store",
        required=False,
        default=1.0,
        help="Sequences that have proportions of ambiguities lower than this value will be resolved, otherwise they will be averaged (RESOLVE only) (Default: 1.0)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity, One copy prints intermediate values and final counts, two copies produces a CSV file with pairwise counts for each non-gap nucleotide",
    )
    parser.add_argument(
        "-l",
        "--minimum_overlap",
        action="store",
        default=500,
        help="What's the minimum amount of overlapping sequence to make a comparison? (Default: 500)",
    )
    parser.add_argument(
        "-n",
        "--ignore_terminal_gaps",
        action="store_true",
        default=False,
        help="Should gaps at the beginning and end of a sequence be ignored (GAPMM only)? (Default: False)",
    )
    parser.add_argument(
        "-j",
        "--json_output",
        action="store_true",
        default=False,
        help="Should the output be in JSON format? (Default: False)",
    )
    return parser


#'tn93 -t 1.0 -f csv -l 250 -o "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.csv" -s "'+str(in_file.split(f)[0])+'.'+str(ref)+'_HXB2.fa" "'+str(in_file.split(f)[0])+'.NoUnkChar.'+str(ref)+'_HXB2.Unidirectional.fa"'
# -t threshhold 1.0 - report all distances less than this
# -f format csv output
# -l overlap
# -s second fasta...  this is the reference we pass in to check if seqs are identical to reference.
# -o output location
# need to recreate parser without input file as a required parameter.
tn93Parser = setup_parser()


# seqId1 - (String) sequence ID of the first sequence
# seq1 - (String) sequence
# seqId2 - (String) sequence ID of the second sequence
# seq2 - (String) sequence
# args - (String) string of args passed in as if on CLI, and then parsed by CLI parser.
# this is an annotation built into databricks to define a Spark UDF
@udf
def tn93_udf(seqId1, seq1, seqId2, seq2, args):
    # print(args)
    args = tn93Parser.parse_args(args.split(" "))
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
    # need to be SeqRecord objects.
    seqRec1 = SeqRecord(seq1, id=seqId1, name=seqId1, description=seqId1)
    seqRec2 = SeqRecord(seq2, id=seqId2, name=seqId2, description=seqId2)
    # result of tn93 is in this format: 'AL00I001279552-9,HXB2_pr,0.0825914'
    # we only need the distance portion.
    return tn93_1.tn93_distance(seqRec1, seqRec2, match_mode).split(",")[2]


# COMMAND ----------

# bealign UDF parameters
import argparse
from functools import partial

from BioExt.args import add_alphabet, add_reference, add_scorematrix

from BioExt.references import hxb2, nl4_3, cov2
from Bio.Seq import Seq, reverse_complement as rc
from BioExt.uds import _align_par
from BioExt import __version__ as version

references = {
    "HXB2_env": hxb2.env,
    "HXB2_gag": hxb2.gag,
    "HXB2_int": hxb2.int,
    "HXB2_nef": hxb2.nef,
    "HXB2_pol": hxb2.pol,
    "HXB2_pr": hxb2.pr,
    "HXB2_prrt": hxb2.prrt,
    "HXB2_rev": hxb2.rev,
    "HXB2_rt": hxb2.rt,
    "HXB2_tat": hxb2.tat,
    "HXB2_vif": hxb2.vif,
    "HXB2_vpr": hxb2.vpr,
    "HXB2_vpu": hxb2.vpu,
    "NL4-3_prrt": nl4_3.prrt,
    "CoV2-3C": cov2.threeC,
    "CoV2-E": cov2.E,
    "CoV2-endornase": cov2.endornase,
    "CoV2-exonuclease": cov2.exonuclease,
    "CoV2-helicase": cov2.helicase,
    "CoV2-leader": cov2.leader,
    "CoV2-methyltransferase": cov2.methyltransferase,
    "CoV2-M": cov2.M,
    "CoV2-N": cov2.N,
    "CoV2-nsp10": cov2.nsp10,
    "CoV2-nsp2": cov2.nsp2,
    "CoV2-nsp3": cov2.nsp3,
    "CoV2-nsp4": cov2.nsp4,
    "CoV2-nsp6": cov2.nsp6,
    "CoV2-nsp7": cov2.nsp7,
    "CoV2-nsp8": cov2.nsp8,
    "CoV2-nsp9": cov2.nsp9,
    "CoV2-ORF10": cov2.ORF10,
    "CoV2-ORF1a": cov2.ORF1a,
    "CoV2-ORF1b": cov2.ORF1b,
    "CoV2-ORF3a": cov2.ORF3a,
    "CoV2-ORF5": cov2.ORF5,
    "CoV2-ORF6": cov2.ORF6,
    "CoV2-ORF7a": cov2.ORF7a,
    "CoV2-ORF7b": cov2.ORF7b,
    "CoV2-ORF8": cov2.ORF8,
    "CoV2-RdRp": cov2.RdRp,
    "CoV2-S": cov2.S,
}


def probability(string):
    try:
        print(string)
        p = float(string)
        print(p)
        if p < 0 or p > 1:
            raise ValueError()
        return p
    except ValueError:
        msg = "'{0}' is not a probability in [0, 1]".format(string)
        raise argparse.ArgumentTypeError(msg)


def add_reference_string(parser, *args):
    from argparse import ArgumentTypeError
    from Bio import SeqIO
    from BioExt.references import hxb2, nl4_3, cov2

    def reference(string):
        if string in references:
            return references[string].load().seq
        else:
            # do a check if it's a valid ref sequence?
            return Seq(string)
        # elif isinstance(string, Seq):
        #    return Seq(string)
        # else:
        #    msg = "'{0}' is not a REFERENCE sequence or string in '{1}' ".format(string, ', '.join(references.keys()))
        #    raise argparse.ArgumentTypeError(msg)

    kwargs = dict(
        metavar="REFERENCE", type=reference, help="REFERENCE sequence or {{{0}}}".format(", ".join(references.keys()))
    )

    parser.add_argument(*args, **kwargs)

    return parser


bealignParser = argparse.ArgumentParser(
    description=("align sequences to a reference using " "a codon alignment algorithm and output to a BAM file")
)
bealignParser.add_argument(
    "-e",
    "--expected-identity",
    type=probability,
    default=None,
    help="discard sequences that are insufficiently identical to the reference",
)
bealignParser.add_argument(
    "-R",
    "--reverse-complement",
    action="store_true",
    help=("also align the reverse complement of each query sequence, " "returning it if the alignment is superior"),
)
bealignParser.add_argument(
    "-K",
    "--keep-reference",
    dest="keep_reference",
    action="store_true",
    help="include the reference sequence as the first one in the resulting BAM file [the default is to strip it]",
)
bealignParser.add_argument(
    "-S",
    "--no-sort",
    dest="sort",
    action="store_false",
    help="do NOT sort the resulting BAM file [the default is to sort]",
)
bealignParser.add_argument("-q", "--quiet", action="store_true", help="do not print status update messages")
# Discard is not supported yet.  We can send this to a new column with the name passed in here.
# bealignParser.add_argument(
#    '-D', '--discard',
#    metavar='DISCARD',
#    type=argparse.FileType('w'),
#    help='discarded sequences are sent to DISCARD'
#    )
add_alphabet(bealignParser, "-a", "--alphabet")
add_scorematrix(bealignParser, "-m", "--score-matrix")
add_reference_string(bealignParser, "-r", "--reference")

# COMMAND ----------

# bealign UDF definition
from BioExt.align import Aligner
from BioExt.misc import compute_cigar, gapful, gapless
from BioExt.io._SamBamIO import _from_seqrecord, _to_cigarstring
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
import builtins
import pysam


def bealign_udf(title, seq, args):
    # print("udf!" + seq)
    # args has to be passed in as a string and split into an array.  I can't seem to pass an array as a column literal for the UDF.
    args = bealignParser.parse_args(args.split(" "))
    # return str(args.reference.seq)
    # if we are keeping the reference, it needs to be added to the sequences before here.
    return bealign_MTNAP(
        title,
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


# needs to return an algined sequence object!
def bealign_MTNAP(
    title,  # input sequence header row from FASTA
    input_sequence,  # Seq Object
    reference,  # seqObject
    expected_identity,  # probability number between 0 and 1 inclusive
    alphabet,  # choices=('amino', 'dna', 'codon'),
    reverse_complement,
    score_matrix,
    # discard_handle,
    do_sort,
    quiet,
    keep_reference,
):
    # print('hi!')
    # print(reference)
    # print(expected_identity)
    # print(reverse_complement)
    # print(score_matrix)
    # print(input_sequence)
    # create SeqRecord - from Bio/SeqIO/FastaIO.py line 239
    try:
        first_word = title.split(None, 1)[0]
    except IndexError:
        assert not title, repr(title)
        # Should we use SeqRecord default for no ID?
        first_word = ""
    input_sequence = SeqRecord(input_sequence, id=first_word, name=first_word, description=title)

    # from bealign script
    try:
        score_matrix_ = score_matrix.load()
    except:
        raise RuntimeError("could not load the score matrix")

    if (alphabet == "dna" and not isinstance(score_matrix, DNAScoreMatrix)) and not isinstance(
        score_matrix, ProteinScoreMatrix
    ):
        raise ValueError(
            "DNA alphabet requires a DNA score matrix, "
            "while amino and codon alphabets require a protein score matrix"
        )

    do_codon = alphabet == "codon"

    # code from _align_par here...  this does the actual alignment.
    aln = Aligner(score_matrix_, do_codon=do_codon, expected_identity=expected_identity)

    if isinstance(reference, str):
        refstr = reference
    elif isinstance(reference, Seq):
        refstr = str(reference)
    elif isinstance(reference, SeqRecord):
        refstr = str(reference.seq)
    else:
        raise ValueError("reference must be one of str, Bio.Seq, Bio.SeqRecord")

    def keep(score, record):
        if aln.expected(score):
            return True
        # elif discard:
        #    discard(record)
        else:
            return False

    def _rc(record):
        if isinstance(record, str):
            return rc(record)
        elif isinstance(record, Seq):
            return record.reverse_complement()
        elif isinstance(record, SeqRecord):
            return SeqRecord(
                record.seq.reverse_complement(), id=record.id, name=record.name, description=record.description
            )
        else:
            raise ValueError("record must be one of str, Bio.Seq, or Bio.SeqRecord")

    # aln, ref, ref_name, and do_revcomp are set by set_globals below
    def _align(record, aln, ref, ref_name, do_revcomp):
        records = (record, _rc(record)) if do_revcomp else (record,)
        score, ref_, record = builtins.max((aln(ref, record) for record in records), key=itemgetter(0))
        record_ = compute_cigar(ref_, record, ref_name)
        return score, record_

    def _to_seqrecord(refname, read):
        seq = Seq(read.seq)

        qname = read.qname

        annotations = {}
        annotations["sam_flag"] = read.flag
        annotations["reference_name"] = refname
        annotations["position"] = read.pos
        annotations["mapping_quality"] = read.mapq
        annotations["CIGAR"] = _to_cigarstring(read.cigar)
        annotations["reference_next"] = read.rnext
        annotations["position_next"] = read.pnext
        annotations["template_length"] = read.tlen

        if read.qual:
            letter_annotations = {}
            letter_annotations["phred_quality"] = [ord(q) - 33 for q in read.qual]
        else:
            letter_annotations = None

        record = SeqRecord(
            seq, id=qname, name=qname, description=qname, annotations=annotations, letter_annotations=letter_annotations
        )

        return record

    (score, record) = _align(input_sequence, aln, refstr, reference, reverse_complement)
    header = None
    samread_from_record = _from_seqrecord(header, record)
    # print(samread_from_record)
    # need the length of the reference sequence.  Normally this is in the bam file header, but we are not using the bam file.
    # we can get it directly.
    # (seq + ('-' * (length - len(seq))))
    # _to_seq_record(samfile, read): takes a AlignedSegment and record and a  pysam.AlignedSegment()
    # I imported the code and changed the parameters.  It only needs the samfile to read the reference name,
    # see https://github.com/veg/BioExt/blob/master/scripts/bam2msa
    # lines 52 and 53
    seq = gapful(_to_seqrecord(record.annotations["reference_name"], samread_from_record), insertions=False)
    seq = seq + (
        "-" * (len(reference) - len(seq))
    )  # pad with gaps on the end of the sequence to match reference length.

    return str(seq.seq)
