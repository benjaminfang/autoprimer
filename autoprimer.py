#! /usr/bin/env python3
import argparse
import sys
from biobrary.bioparse.gtf_parse import GTF
from biobrary.bioparse.fasta_parse import FASTA
from biobrary.misc import change_coordinate
from biobrary.seq import Seq

def getargs_info(cmd_args):
    parser = argparse.ArgumentParser(prog="autoprimer info",
                    description="print information of gene")

    parser.add_argument("--fasta", "-F", help="genome fasta file", required=True)
    parser.add_argument("--gtf", "-G", help="genome gtf file", required=True)
    parser.add_argument("--extend_up", type=int, default=0, help="extend gene upstream")
    parser.add_argument("--extend_down", type=int, default=0, help="extend gene downstream")
    parser.add_argument("--add_geneseq", "-A", action='store_true', help="add gene sequence information")
    parser.add_argument("--add_transseq", "-B", action='store_true', help="add transcript sequence information")
    parser.add_argument("--add_cdsseq", "-D", action='store_true', help="add CDS sequence information")

    parser.add_argument("--outfile", "-O", help="name of output file.",
                        default= None)
    parser.add_argument("geneid", nargs="+", help="ncbi sequence id")

    args = parser.parse_args(cmd_args)
    return  (args.gtf, args.fasta, args.outfile, args.geneid,
             args.extend_up, args.extend_down,
             args.add_geneseq, args.add_transseq, args.add_cdsseq)


def getargs_primer(cmd_args):
    pass


def get_info(args):
    (gtf, fasta, outfile, geneids,
     extend_up, extend_down,
     add_geneseq, add_transseq, add_cdsseq) = getargs_info(args)
    fasta = FASTA(fasta)
    gtf = GTF(gtf)

    if outfile:
        outfile_fout = open(outfile, "w")
    else:
        outfile_fout = sys.stdout

    for geneid in geneids:
        gene = gtf.get_gene(geneid)
        if not gene:
            continue
        geneid = gene.get_geneid()
        geneinfo = gene.get_info()
        seqname = geneinfo[0]
        pos = geneinfo[2][0]
        ori = geneinfo[3]

        seq_len = str(pos[1] - pos[0] + 1)
        print("\t".join([">" + geneid, seqname,
                ",".join([str(ele) for ele in pos]), ori, seq_len]),
                file=outfile_fout)

        seq = fasta.get_seq(seqname)
        if ori == "+":
            if extend_up > 0:
                up_seq = seq.substr(pos[0] - extend_up, pos[0] - 1)
            else:
                up_seq = ""
            if extend_down > 0:
                down_seq = seq.substr(pos[1] + 1, pos[1] + extend_down)
            else:
                down_seq = ""
            gene_seq = seq.substr(pos[0], pos[1])
        else:
            if extend_up > 0:
                up_seq = seq.substr(pos[1] + 1, pos[1] + extend_up)
                up_seq = Seq(up_seq).reverse_comp().seq
            else:
                up_seq = ""
            if extend_down > 0:
                down_seq = seq.substr(pos[0] - extend_up, pos[0] - 1)
                down_seq = Seq(down_seq).reverse_comp().seq
            else:
                down_seq = ""
            gene_seq = Seq(seq.substr(pos[0], pos[1])).reverse_comp().seq

        if add_geneseq:
            print("\t".join(["&UP", up_seq]), file=outfile_fout)
            print("\t".join(["&GENESEQ", gene_seq]), file=outfile_fout)
            print("\t".join(["&DOWN", down_seq]), file=outfile_fout)

        if ori == "+":
            ref_pos = geneinfo[2][0][0]
        else:
            ref_pos = geneinfo[2][-1][-1]

        trans = gene.get_transcript()
        for tran in trans:
            if tran.get_gbkey() == "mRNA":
                transid = tran.get_transid()
                traninfo = tran.get_info()
                pos = traninfo[2]
                ori = traninfo[3]
                pos = change_coordinate(ref_pos, pos, coor="relative", ori=ori)
                trans_len = pos[0][1] - pos[0][0] + 1
                print("\t".join(["@" + transid,
                    ",".join([str(ele) for ele in pos[0]]), ori, str(trans_len)]),
                    file=outfile_fout)
                if add_transseq:
                    trans_seq = gene_seq[pos[0] - 1: pos[1]]
                    print("\t".join(["&TRANSEQ", trans_seq]), file=outfile_fout)

                exon, cds, start_codon, stop_codon = tran.get_child()

                if start_codon:
                    start_info = start_codon.get_info()
                    pos = start_info[2]
                    ori = start_info[3]
                    pos = change_coordinate(ref_pos, pos, coor="relative", ori=ori)
                    start_up = pos[0][0]
                    start_len = pos[-1][-1] - pos[0][0] + 1
                    pos = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    start_line = "\t".join(["$START", ";".join(pos), ori, str(start_len)])

                if stop_codon:
                    stop_info = stop_codon.get_info()
                    pos = stop_info[2]
                    ori = stop_info[3]
                    pos = change_coordinate(ref_pos, pos, coor="relative", ori=ori)
                    stop_up = pos[0][0]
                    stop_len = pos[-1][-1] - pos[0][0] + 1
                    pos = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    stop_line = "\t".join(["$STOP", ";".join(pos), ori], str(stop_len))

                if cds:
                    cds_info = cds.get_info()
                    pos = cds_info[2]
                    ori = cds_info[3]
                    pos = change_coordinate(ref_pos, pos, coor="relative", ori=ori)
                    cds_len = 0
                    for pp in pos:
                        pp += pp[1] - pp[0] + 1
                    pos_str = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    cds_line = "\t".join(["$CDS", ";".join(pos_str), ori, str(cds_len)])
                    cds_seq = ""
                    if add_cdsseq:
                        for pp in pos:
                            cds_seq += gene_seq[pp[0] - 1: pp[1]]

                if exon:
                    exoninfo = exon.get_info()
                    pos = exoninfo[2]
                    ori = exoninfo[3]
                    pos = change_coordinate(ref_pos, pos, coor="relative", ori=ori)
                    exon_len = 0
                    for pp in pos:
                        exon_len += pp[1] - pp[0] + 1
                    pos_str = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    exon_line = "\t".join(["$EXON", ";".join(pos_str), ori, str(exon_len)])

                    exon_seq = []
                    exon_seq_anno = []
                    for pp in pos:
                        seq_piece = gene_seq[pp[0] - 1: pp[1]]
                        exon_seq.append(seq_piece)
                        if start_up <= pp[1] and start_up >= pp[0]:
                            seq_piece = list(seq_piece)
                            seq_piece.insert(start_up - pp[0], "|>" + str(start_up) + "|")
                            seq_piece = "".join(seq_piece)
                        if stop_up <= pp[1] and stop_up >= pp[0]:
                            seq_piece = list(seq_piece)
                            seq_piece.insert(stop_up - pp[0], "|<" + str(stop_up) + "|")
                            seq_piece = "".join(seq_piece)
                        exon_seq_anno.append("".join(["[" + str(pp[0]) + "]", seq_piece, "[" + str(pp[1]) + "]"]))

                    exon_seq = "".join(exon_seq)
                    exon_seq_anno = "".join(exon_seq_anno)


                print(exon_line, file=outfile_fout)
                print(exon_seq_anno, file=outfile_fout)

                print(cds_line, file=outfile_fout)
                if add_cdsseq:
                    print(cds_seq, file=outfile_fout)

                print(start_line, file=outfile_fout)
                print(stop_line, file=outfile_fout)

            if outfile:
                outfile_fout.close()

def get_primer():
    pass


def print_usage():
    print("Usage:")
    print("autoprimer subcommand -h/args\n")
    print("subcommand:")
    print("    info: print information of gene.")
    print("    primer: design primer for gene.")
    print()


def main():
    args = sys.argv
    if len(args) < 2 or args[1] == "-h" or args[1] == "--help":
        print_usage()
        exit()

    if args[1] == "info":
        get_info(args[2:])

    elif args[1] == "primer":
        get_primer(args[2:])
    else:
        print_usage()
        exit()


if __name__ == "__main__":
    main()
