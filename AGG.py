import gzip
import json

class GraphicalGenome:
    def __init__(self, filename):
        self.anchor, self.edges, self.outgoing, self.incoming = self.loadgraph(filename)

    def loadgraph(self, filename):
        if filename.endswith(".gz"):
            with gzip.open(filename, 'rb') as fp:
                data = fp.readlines()
        else:
            with open(filename, 'rb') as fp:
                data = fp.readlines()
        Anchor_dict = {}
        Edge_dict = {}
        Outgoing = {}
        Incoming = {}

        for line in data:
            line = line.decode()[:-1]
            if line.startswith("S"):
                itemlist = line.split('\t')
                name = itemlist[1]
                seq = itemlist[2]
                annotation = itemlist[3][5:]
                if name.startswith("A"):
                    Anchor_dict[name] = json.loads(annotation)
                    Anchor_dict[name]['seq'] = seq
                elif name.startswith("E"):
                    Edge_dict[name] = json.loads(annotation)
                    Edge_dict[name]['seq'] = seq
            elif line.startswith("L"):
                itemlist = line.split('\t')
                src = itemlist[1]
                dst = itemlist[3]
                Outgoing[src] = Outgoing.get(src, []) + [dst]
                Incoming[dst] = Incoming.get(dst, []) + [src]
                    
        return Anchor_dict, Edge_dict, Outgoing, Incoming
    
def reconstruct_refpath(graph, ref = 'NC_012920.1'):
    src = "SOURCE"
    Path = []

    visited = set()

    while src != "SINK":
        edgelist = graph.outgoing[src]
        ref_edge = []
        for edge in edgelist:
            if ref in graph.edges[edge]['reads']:
                ref_edge += [edge]
        if len(ref_edge) == 1:
            edge = ref_edge[0]
        else:
            for edge in ref_edge:
                if graph.outgoing[edge] == "SINK":
                    break
        Path += [src, edge]
        visited.add(src)
        src = graph.outgoing[edge][0]
        assert src.startswith("A") or src == "SINK"
        if src in visited:
            break

    reconstruct = ''
    for item in Path:
        if item == "SOURCE":
            reconstruct += ""
        elif item.startswith('A'):
            reconstruct += graph.anchor[item]['seq']
        elif item.startswith("E"):
            reconstruct += graph.edges[item]['seq']
        else:
            print(item)
    return reconstruct, Path

def write_gfa(AnchorInfo, Edge_info, outputfilename):
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

def reconstruct_sample_subgraph(graph, samplelist, outputfile):
    Anchor_dict = {}
    Edge_dict = {}
    edgelist = list(graph.edges.keys())
    Anchor_dict = {}
    Edge_dict = {}
    nodelist = []
    for edge in edgelist:
        if len(set(samplelist) & set(graph.edges[edge]['strain'])) > 0:
            Edge_dict[edge] = graph.edges[edge]
            Edge_dict[edge]['src'] = graph.incoming[edge][0]
            Edge_dict[edge]['dst'] = graph.outgoing[edge][0]
        nodelist += graph.incoming[edge]
        nodelist += graph.outgoing[edge]

    nodelist = list(set(nodelist))
    for node in nodelist:
        if node.startswith("A"):
            Anchor_dict[node] = graph.anchor[node]

    write_gfa(Anchor_dict, Edge_dict, outputfile)

    return Anchor_dict, Edge_dict

def find_most_supported_path(graph, sample, ref = '012920.1'):
    src = "SOURCE"
    Path = []

    visited = set()

    while src != "SINK":
        edgelist = graph.outgoing[src]
        sample_edge = []
        for edge in edgelist:
            if sample in set(graph.edges[edge]['strain']):
                sample_edge += [edge]
        if len(sample_edge)<1: # no haplotype for this sample at this place the go to the reference one
            for edge in edgelist:
                if ref in set(graph.edges[edge]['strain']):
                    break

        elif len(sample_edge) == 1:
            edge = sample_edge[0]
        else:
            # find the most supported edge
            read_count = [len(graph.edges[e]['reads']) for e in sample_edge]

            index = numpy.argmax(read_count)

            edge = sample_edge[index]


        Path += [src, edge]

        visited.add(src)

        src = graph.outgoing[edge][0]

        assert src.startswith("A") or src == "SINK"

        if src in visited:
            break

    reconstruct = ''
    for item in Path:
        if item == "SOURCE":
            reconstruct += ""
        elif item.startswith('A'):
            reconstruct += graph.anchor[item]['seq']
        elif item.startswith("E"):
            reconstruct += graph.edges[item]['seq']
        else:
            print(item)

    return reconstruct, Path