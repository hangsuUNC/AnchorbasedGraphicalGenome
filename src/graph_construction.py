import sys
sys.path.append("/home/suhang/jupyter/jupyter_env/lib/python3.11/site-packages")
import gzip
import pyfrost
import gzip
import argparse
from collections import defaultdict
import json
import numpy


class GraphicalGenome: 
        
    def read_reference(self, filename):
        
        if filename.endswith(".gz"):
            with gzip.open(filename, "rb") as fp:
                data =fp.readlines()
        else:
            with open(filename, 'rb') as fp:
                data = fp.readlines()
                
        kmers = [km.decode()[:-1] for km in data[1:]]
        contig = "".join(kmers)
        contig = "+" + contig # 1-based coordinates
        return contig
        
        
    def revised_find_shared_kmers(self, g, k, colorsnum):
        anchorlist = []
        no_colors = 0
        for n, data in g.nodes(data=True):
            unitig = data['unitig_sequence']
            for i in range(len(unitig)- k + 1):
                n = unitig[i:i+k]
                if "N" in n:
                    continue
                try:
                    colors = g.nodes[pyfrost.Kmer(n)]['colors'] 
                    if len(set(colors)) == colorsnum:
                        anchorlist.append(n)
                        break
                except:
                    no_colors += 1

        return anchorlist, no_colors

    def map_to_genome(self, contig, anchorlist, k):
        PositionDict = defaultdict(list)
        for anchor in anchorlist:
            PositionDict[anchor]
        for i in range(1, len(contig) - k + 1):
            kmer = contig[i:i+k]
            if kmer in PositionDict:
                PositionDict[kmer] = PositionDict.get(kmer, []) + [i]
        return PositionDict

    def create_anchors(self, PositionDict, k):
        anchor_updated_list = []

        for kmer in PositionDict:
            if len(PositionDict[kmer]) == 1:
                anchor_updated_list.append(kmer)
        AnchorInfo = {}
        for kmer in anchor_updated_list:
            pos = PositionDict[kmer][0]
            anchorname = "A%06d" %  (pos//k + 1)
            AnchorInfo[anchorname] = {}
            AnchorInfo[anchorname]['seq'] = kmer
            AnchorInfo[anchorname]['pos'] = pos

        anchornames = sorted(AnchorInfo.keys())
        anchor_unadjacent_list = []
        index = 0
        sanchor = anchornames[index]
        while sanchor < anchornames[-1]:
            for danchor in anchornames[index+1:]:
                if AnchorInfo[danchor]['pos'] - AnchorInfo[sanchor]['pos'] > k+1:
                    break
                index += 1
            anchor_unadjacent_list += [sanchor, danchor]
            sanchor = danchor
        anchor_unadjacent_list = sorted(set(anchor_unadjacent_list))[:-1]
        return AnchorInfo, anchor_unadjacent_list


    def loadFasta(self, filename):
        """ Parses a classically formatted and possibly 
            compressed FASTA file into a list of headers 
            and fragment sequences for each sequence contained.
            The resulting sequences are 0-indexed! """
        if (filename.endswith(".gz")):
            fp = gzip.open(filename, 'rb')
        else:
            fp = open(filename, 'rb')
        # split at headers
        data = fp.read().decode().split('>')
        fp.close()
        # ignore whatever appears before the 1st header
        data.pop(0)     
        headers = []
        sequences = []
        for sequence in data:
            lines = sequence.split('\n')
            headers.append(lines.pop(0))
            sequences.append(''.join(lines))
        return (headers, sequences)



    def mapping_info(self, AnchorInfo, contig, k):
        anchor_list = list(AnchorInfo.keys())
        Anchorseq = {AnchorInfo[anchor]['seq']:anchor for anchor in anchor_list}

        seqlist = Anchorseq.keys()
        PositionDict = defaultdict(list)
        for anchor_seq in seqlist:
            anchor_rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(anchor_seq)])
            PositionDict[anchor_seq]
            PositionDict[anchor_rev]

        for i in range(1, len(contig) - k + 1):
            kmer = contig[i:i+k]
            if kmer in PositionDict:
                PositionDict[kmer] = PositionDict.get(kmer, []) + [i]

        A = {}
        SVs = {}
        for anchor, D in AnchorInfo.items():
            anchor_seq = D['seq']
            anchor_rev = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(anchor_seq)])
            poslist = PositionDict[anchor_seq] + PositionDict[anchor_rev]
            if len(poslist) == 1:
                A[anchor] = poslist[0]
            else:
                SVs[anchor] = poslist

        return A, SVs

    def find_edge_info(self, src_pos, dst_pos, k, contig, contigname, sample, Anchorseq):
        E = {} # edgeinfo
        # get source infomation
        if src_pos == 0:
            src = "SOURCE"
            src_seq = ""
            pr = False
        else:
            src_seq = contig[src_pos:src_pos + k]
            try:
                src = Anchorseq[src_seq]
                pr = False
            except:
                src_seq = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(src_seq)])
                src = Anchorseq[src_seq]
                pr = True

        dst_seq = contig[dst_pos:dst_pos+k]

        if dst_pos == len(contig): # fix sink edge issue 08/29/23
            dst = "SINK"
            dst_seq = ""
            sr = True
        else:
            try:
                dst = Anchorseq[dst_seq]
                sr = False
            except:
                dst_seq = ''.join([{'A':'T','C':'G','G':'C','T':'A'}[base] for base in reversed(dst_seq)])
                dst = Anchorseq[dst_seq]
                sr = True

        if src_pos == 0:
            edge_seq = contig[src_pos:dst_pos] # first edge fix bug 08/28/2023
        else:
            edge_seq = contig[src_pos+k:dst_pos] # fix bug 08/23/2023

        if pr and sr:
            edge_seq = ''.join([{'A':'T','C':'G','G':'C','T':'A', 'N':'N'}[base] for base in reversed(edge_seq)])
            node = src
            src = dst
            dst = node


        E = {}
        E['seq'] = edge_seq
        E['src'] = src
        E['dst'] = dst
        E['reads'] = [contigname]
        E['strain'] = [sample]


        return E

    def create_edgefile(self, headers, sequences, AnchorInfo, k):
        Edge_info = {}
        Outgoing = {}

        edgenum_perread = []
        anchor_list = list(AnchorInfo.keys())
        Anchorseq = {AnchorInfo[anchor]['seq']:anchor for anchor in anchor_list}

        for contig_index in range(len(headers)):
            contig = sequences[contig_index]
            contig_name = headers[contig_index]
            sample_name = contig_name.split('_')[-1]
            A, SVs = self.mapping_info(AnchorInfo, contig, k)
            splitposlist = sorted(A.values())
            edgeindex = 0 # 
            src_pos = 0
            for dst_pos in splitposlist:
                E = self.find_edge_info(src_pos, dst_pos, k, contig, contig_name, sample_name, Anchorseq)
                src = E['src']
                edgelist = Outgoing.get(src, [])
                for edge in edgelist:
                    if Edge_info[edge]['dst'] != E['dst']:
                        continue
                    if Edge_info[edge]['seq'] == E['seq']:
                        Edge_info[edge]['reads'] += E['reads']
                        Edge_info[edge]['strain'] = sorted(set(Edge_info[edge]['strain']) | set(E['strain']))
                        break
                else:
                    edgename = "E%05d.%04d" % (contig_index, edgeindex)
                    Edge_info[edgename] = E
                    Outgoing[src] = Outgoing.get(src,[]) + [edgename]
                    edgeindex += 1
                # update
                src_pos = dst_pos

            dst_pos = len(contig) # fix bug 08/29/2023 
            E = self.find_edge_info(src_pos, dst_pos, k, contig, contig_name, sample_name, Anchorseq)
            src = E['src']
            edgelist = Outgoing.get(src, [])
            for edge in edgelist:
                if Edge_info[edge]['dst'] != E['dst']:
                    continue
                if Edge_info[edge]['seq'] == E['seq']:
                    Edge_info[edge]['reads'] += E['reads']
                    Edge_info[edge]['strain'] = sorted(set(Edge_info[edge]['strain']) | set(E['strain']))
                    break
            else:
                edgename = "E%05d.%04d" % (contig_index, edgeindex)
                Edge_info[edgename] = E
                Outgoing[src] = Outgoing.get(src,[]) + [edgename]
                edgeindex += 1

            edgenum_perread.append(edgeindex)   

        return Edge_info, Outgoing

    def filter_undersupported_edges(self, Edge_info, ref, threshold = 4):
        # filter edges that supported by > 3 reads
        deleted_edge = []
        for edge in Edge_info:
            if ref in Edge_info[edge]['reads']:
                continue
            if len(Edge_info[edge]['reads']) < 4:
                deleted_edge.append(edge)
        print(len(deleted_edge))
        for edge in deleted_edge:
            del Edge_info[edge]
        return 

    def write_gfa(self, AnchorInfo, Edge_info, outputfilename):
        header = ['H\tVN:Z:1.0\n']
        # add json annotation to edge and anchor sequences
        AnchorS = []
        for anchor,D in AnchorInfo.items():
            seq = D.pop('seq')
            json_string = json.dumps(D)
            AnchorS += ['S\t%s\t%s\t%s\n' % (anchor, seq, "PG:J:" + json_string)]

        EdgeS = []   
        Link = [] 
        for edge,edge_dict in Edge_info.items():
            seq = edge_dict.pop('seq')
            src = edge_dict.pop('src')
            dst = edge_dict.pop('dst')

            json_string = json.dumps(edge_dict)
            EdgeS += ['S\t%s\t%s\t%s\t%s\n' % (edge, seq, "PG:J:" + json_string, "RC:i:" + str(len(edge_dict['reads'])))]
            Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (src, "+", edge, "+", "0M"))
            Link.append('L\t%s\t%s\t%s\t%s\t%s\n'% (edge, "+", dst, "+", "0M"))     

        with open(outputfilename, 'w') as fp:
            for h in header:
                fp.write(h)
            for s in AnchorS:
                fp.write(s)
            for s in EdgeS:
                fp.write(s)
            for l in Link:
                fp.write(l)

def main():
    # create parser object
    parser = argparse.ArgumentParser(description = "Construct Anchor-based Graphical Genome")

    
    parser.add_argument("-i", "--Bifrost_graph_input", type = str, nargs = 1,
                        metavar = "input bifrost graph", default = None,
                        help = "a path for bifrost graph")
    parser.add_argument("-f", "--ReferenceFa", type = str, nargs = 1,
                    metavar = "Reference fasta file", help = "path for linear reference fasta file")
    parser.add_argument("-n", "--referencename", type = str, nargs = 1,
                    metavar = "a string for reference accession", help = "a reference id")
    parser.add_argument("-r", "--readFa", type = str, nargs = 1,
                    metavar = "a fasta file for all input reads with readname_sample information in header", help = "a fasta file for all input reads with readname_sample information in header")
    
    parser.add_argument("-t", "--ReadcountThreshold", type = int, nargs = 1,
                    metavar = "an integer, edges with less than t reads supported are purged, unless its on reference assemblies")
    parser.add_argument("-o", "--Output", type = str, nargs = 1,
                    metavar = "output filename", help = "output gfa filename")
    
    args = parser.parse_args()
    bifrost_graph = args.Bifrost_graph_input[0]
    ref_fa = args.ReferenceFa[0]
    ref_name = args.referencename[0]
    read_fa = args.readFa[0]
    outputgfa = args.Output[0]
    threshold = args.ReadcountThreshold[0]
    
    g = GraphicalGenome()
    
    bifrostg = pyfrost.load(bifrost_graph)
    colorsnum = len(bifrostg.graph['color_names'])
    k = bifrostg.graph['k']
    anchorlist, no_colors = g.revised_find_shared_kmers(bifrostg, k, colorsnum)
    ref_contig = g.read_reference(ref_fa)
    PositionDict = g.map_to_genome(ref_contig, anchorlist, k)
    AnchorInfo, anchor_unadjacent_list = g.create_anchors(PositionDict, k)
    Final_anchor = {anchor:AnchorInfo[anchor] for anchor in anchor_unadjacent_list}
    
    pos1 = [AnchorInfo[anchor]['pos'] for anchor in anchor_unadjacent_list[1:]]
    pos = [AnchorInfo[anchor]['pos'] for anchor in anchor_unadjacent_list[:-1]]
    assert len(numpy.where((numpy.array(pos1) - numpy.array(pos))<=k)[0]) == 0

    headers, sequences = g.loadFasta(read_fa)
    headers.append(ref_name)
    sequences.append(ref_contig[1:])

    Edge_info, Outgoing = g.create_edgefile(headers, sequences, Final_anchor, k)
    g.filter_undersupported_edges(Edge_info, ref_name, threshold)
    g.write_gfa(Final_anchor, Edge_info, outputgfa)
    
   

if __name__ == "__main__":
    # calling the main function
    main()

