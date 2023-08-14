package algo;

import io.GFAWriter;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.KmerUtils;
import structures.ConnectedComponent;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static io.GFAWriter.normalizeDna;

public class Comp2Graph implements Runnable {
    private final int k;
    private final int compID;
    private final Logger logger;

    private final ByteArrayOutputStream outputStream;

    private final ConnectedComponent comp;
    private final BigLong2ShortHashMap all_kmers;
    private final BigLong2ShortHashMap comp_kmers;
    private final Map<String, Short> subgraph;

    private int size;
    private SingleNode[] nodes;

    public Comp2Graph(int k, int icomp, ConnectedComponent comp, BigLong2ShortHashMap all_kmers, ByteArrayOutputStream outputStream, Logger logger) {
        this.k = k;
        this.compID = icomp;
        this.logger = logger;

        this.comp = comp;
        this.all_kmers = all_kmers;
        this.outputStream = outputStream;
        this.comp_kmers = new BigLong2ShortHashMap(8, 12, true);
        this.subgraph = new HashMap<>();
    }


    @Override
    public void run() {
        buildComponent();
        initStructures();
        mergePaths();
        createGFA();
    }

    private void buildComponent() {
        logger.debug("Building graph component...");
        if (all_kmers != null) {
            for (long kmer: comp.kmers) {
                comp_kmers.put(kmer, all_kmers.getWithZero(kmer));
            }
        } else {
            for (long kmer: comp.kmers) {
                comp_kmers.put(kmer, (short)1);
            }
        }

        for (long kmer: comp.kmers) {
            subgraph.put(normalizeDna(KmerUtils.kmer2String(kmer, k)), comp_kmers.getWithZero(kmer));
        }
    }

    private void initStructures() {
        logger.debug("Initializing structures for graph traversal...");
        Map<String, List<SingleNode>> nodeByKmer = new HashMap<>();
        size = subgraph.size() * 2;
        nodes = new SingleNode[size];
        {
            int id = 0;
            for (Map.Entry<String, Short> stringShortEntry : subgraph.entrySet()) {
                String seq = stringShortEntry.getKey();
                String rc = DnaTools.reverseComplement(seq);
                nodes[id] = new SingleNode(seq, id);
                nodes[id + 1] = new SingleNode(rc, id);
                nodes[id].rc = nodes[id + 1];
                nodes[id + 1].rc = nodes[id];
                id += 2;
            }
        }
        for (int i = 0; i < size; i++) {
            String key = nodes[i].sequence.substring(0, k - 1);
            if (!nodeByKmer.containsKey(key)) {
                nodeByKmer.put(key, new ArrayList<>());
            }
            nodeByKmer.get(key).add(nodes[i]);
        }
        for (int i = 0; i < size; i++) {
            String lastK = nodes[i].sequence.substring(1);
            if (nodeByKmer.containsKey(lastK)) {
                nodes[i].rc.neighbors.addAll(nodeByKmer.get(lastK));
            }
        }
    }

    private void mergePaths() {
        logger.debug("Merging graph paths...");
        while (true) {
            boolean acted = false;
            for (int i = 0; i < size; i++) {
                if (!nodes[i].deleted && nodes[i].neighbors.size() == 1) {
                    SingleNode other = nodes[i].neighbors.get(0);
                    if (other.neighbors.size() != 1) {
                        continue;
                    }
                    mergeNodes(nodes[i], other);
                    acted = true;
                }
            }
            if (!acted) {
                break;
            }
        }
    }

    private void mergeNodes(SingleNode firstPlus, SingleNode secondMinus) {
        // first k-1 symbols of firstPlus coincide with complement of first k-1 symbols of secondMinus
        SingleNode firstMinus = firstPlus.rc, secondPlus = secondMinus.rc;
        String newSeq = mergeLabels(secondPlus.sequence, firstPlus.sequence);
        String newSeqRC = mergeLabels(firstMinus.sequence, secondMinus.sequence);

        secondPlus.sequence = newSeq;
        firstMinus.sequence = newSeqRC;
        secondPlus.rc = firstMinus;
        firstMinus.rc = secondPlus;

        firstPlus.deleted = secondMinus.deleted = true;
    }

    private void checkLabels(String a, String b) {
        if (!a.substring(a.length() - (k - 1)).equals(b.substring(0, k - 1))) {
            throw new AssertionError("Labels should be merged, but can not: " + a + " and " + b);
        }
    }

    private String mergeLabels(String a, String b) {
        checkLabels(a, b);
        return a + b.substring(k - 1);
    }

    private void createGFA() {
        logger.debug("Preparing GFA output...");
        GFAWriter writer = new GFAWriter(k, compID, outputStream, nodes, subgraph);
        writer.generateOutput();
    }
}
