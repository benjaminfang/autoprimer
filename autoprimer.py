#! /usr/bin/env python3
import argparse
import sys
from biobrary.bioparse.gtf_parse import GTF
from biobrary.bioparse.fasta_parse import FASTA
from biobrary.misc import change_coordinate
from biobrary.seq import Seq

def getargs_info(cmd_args):
    parser = argparse.ArgumentParser(prog="autoprimer",
                    description="A primer3 designer based on primer3.")

    parser.add_argument("--fasta", "-F", help="genome fasta file", required=True)
    parser.add_argument("--gtf", "-G", help="genome gtf file", required=True)
    parser.add_argument("--extend_up", type=int, default=0, help="extend gene upstream")
    parser.add_argument("--extend_down", type=int, default=0, help="extend gene downstream")
    parser.add_argument("--add_geneseq", "-A", action='store_true', help="add gene sequence information")
    parser.add_argument("--add_transseq", "-B", action='store_true', help="add transcript sequence information")
    parser.add_argument("--add_exonseq", "-C", action='store_true', help="add exon sequence informamtion")
    parser.add_argument("--anno_exonseq", action='store_true', help="add exon sequence informamtion")
    parser.add_argument("--add_cdsseq", "-D", action='store_true', help="add CDS sequence information")
    
    parser.add_argument("--outfile", "-O", help="name of output file.",
                        default= None)
    parser.add_argument("geneid", nargs="+", help="ncbi sequence id")

    args = parser.parse_args(cmd_args)
    return  (args.gtf, args.fasta, args.outfile, args.geneid,
             args.extend_up, args.extend_down,
             args.add_geneseq, args.add_transseq, args.add_exonseq, args.add_cdsseq,
             args.anno_exonseq)


def getargs_primer():
    pass


def get_info(args):
    (gtf, fasta, outfile, geneids,
     extend_up, extend_down,
     add_geneseq, add_transseq, add_exonseq, add_cdsseq,
     anno_exonseq) \
        = getargs_info(args)
    fasta = FASTA(fasta)
    gtf = GTF(gtf)
    for geneid in geneids:
        gene = gtf.get_gene(geneid)
        geneid = gene.get_geneid()
        geneinfo = gene.get_info()
        seqname = geneinfo[0]
        pos = geneinfo[2][0]
        ori = geneinfo[3]
        print("\t".join([">" + geneid, seqname, ",".join([str(ele) for ele in pos]), ori]))
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
            print("\t".join(["&UP", up_seq]))
            print("\t".join(["&GENESEQ", gene_seq]))
            print("\t".join(["&DOWN", down_seq]))
                
        gene_left = geneinfo[2][0][0]
        gene_right = geneinfo[2][-1][-1]
        trans = gene.get_transcript()
        for tran in trans:
            if tran.get_gbkey() == "mRNA":
                transid = tran.get_transid()
                traninfo = tran.get_info()
                pos = traninfo[2]
                ori = traninfo[3]
                pos = change_coordinate(gene_left, pos, coor="relative", ori=ori)
                print("\t".join(["@" + transid, ",".join([str(ele) for ele in pos[0]]), ori]))
                if add_transseq:
                    trans_seq = gene_seq[pos[0] - 1: pos[1]]
                    print("\t".join(["&TRANSEQ", trans_seq]))

                exon, cds, start_codon, stop_codon = tran.get_child()

                if start_codon:
                    start_info = start_codon.get_info()
                    pos = start_info[2]
                    ori = start_info[3]
                    pos = change_coordinate(gene_left, pos, coor="relative", ori=ori)
                    start_up = pos[0][0]
                    pos = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    start_line = "\t".join(["$START", ";".join(pos), ori])

                if stop_codon:
                    stop_info = stop_codon.get_info()
                    pos = stop_info[2]
                    ori = stop_info[3]
                    pos = change_coordinate(gene_left, pos, coor="relative", ori=ori)
                    stop_up = pos[0][0]
                    pos = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    stop_line = "\t".join(["$STOP", ";".join(pos), ori])

                if cds:
                    cds_info = cds.get_info()
                    pos = cds_info[2]
                    ori = cds_info[3]
                    pos = change_coordinate(gene_left, pos, coor="relative", ori=ori)
                    pos_str = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    cds_line = "\t".join(["$CDS", ";".join(pos_str), ori])
                    cds_seq = ""
                    if add_cdsseq:
                        for pp in pos:
                            cds_seq += gene_seq[pp[0] - 1: pp[1]]
                 
                if exon:
                    exoninfo = exon.get_info()
                    pos = exoninfo[2]
                    ori = exoninfo[3]
                    pos = change_coordinate(gene_left, pos, coor="relative", ori=ori)
                    pos_str = [str(ele[0]) + "," + str(ele[1]) for ele in pos]
                    exon_line = "\t".join(["$EXON", ";".join(pos_str), ori])
                    exon_seq = []
                    exon_seq_anno = []
                    if add_exonseq:
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

                         
                print(exon_line)
                print(exon_seq_anno)
                print(exon_seq)

                print(cds_line)
                print(cds_seq)

                print(start_line)
                print(stop_line)

              



def get_primer():
    pass


def main():
    args = sys.argv
    if len(args) < 2:
        print("Usage")
        exit()

    if args[1] == "info":
        get_info(args[2:])
        
    elif args[1] == "primer":
        get_primer(args[2:])
    else:
        print("Error, subcommand not recognized.", file=sys.stderr)




if __name__ == "__main__":
    main()
