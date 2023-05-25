from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, reverse_complement as rc
from BioExt.align import Aligner
from BioExt.misc import compute_cigar, gapful, gapless
from BioExt.io._SamBamIO import _from_seqrecord, _to_cigarstring
from BioExt.args import add_alphabet, add_reference, add_scorematrix
from BioExt.scorematrices import DNAScoreMatrix, FrequenciesError, ProteinScoreMatrix
from BioExt.references import hxb2, nl4_3, cov2
from BioExt.uds import _align_par
from BioExt import __version__ as version
import pysam
from functools import partial
from operator import itemgetter
import builtins
import argparse

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
add_alphabet(bealignParser, "-a", "--alphabet")
add_scorematrix(bealignParser, "-m", "--score-matrix")
add_reference_string(bealignParser, "-r", "--reference")


def bealign(seq, ref, do_reverse_complement: bool = False, args: str = "-m HIV_BETWEEN_F") -> str:
    args = bealignParser.parse_args(args.split(" "))
    # return str(args.reference.seq)
    # if we are keeping the reference, it needs to be added to the sequences before here.
    return bealign_MTNAP(
        "NA",
        Seq(seq),
        ref,
        args.expected_identity,
        args.alphabet,
        do_reverse_complement,
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

    # do_codon = alphabet == "codon"
    do_codon = False

    # code from _align_par here...  this does the actual alignment.
    aln = Aligner(score_matrix_, do_codon=do_codon, expected_identity=expected_identity)

    if isinstance(reference, str):
        refstr = reference
    elif isinstance(reference, Seq):
        refstr = str(reference)
    elif isinstance(reference, SeqRecord):
        refstr = str(reference.seq)
    else:
        print(type(reference))
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
    seq = gapful(_to_seqrecord(record.annotations["reference_name"], samread_from_record), insertions=False)
    seq = seq + (
        "-" * (len(reference) - len(seq))
    )  # pad with gaps on the end of the sequence to match reference length.

    return str(seq.seq)
