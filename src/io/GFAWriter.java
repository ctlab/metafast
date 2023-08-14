package io;

import algo.SingleNode;


import java.io.ByteArrayOutputStream;
import java.io.PrintWriter;
import java.util.Map;

import static ru.ifmo.genetics.dna.DnaTools.reverseComplement;

public class GFAWriter {
    private final int k;
    private final int compID;
    private final ByteArrayOutputStream outputStream;
    private final SingleNode[] nodes;
    private final Map<String, Short> subgraph;
    private final int size;


    public GFAWriter(int k, int compID, ByteArrayOutputStream outputStream, SingleNode[] nodes, Map<String, Short> subgraph) {
        this.k = k;
        this.compID = compID;
        this.outputStream = outputStream;
        this.nodes = nodes;
        this.subgraph = subgraph;
        this.size = nodes.length;
    }

    public void generateOutput() {
        PrintWriter out = new PrintWriter(outputStream, false);
        for (int i = 0; i < size; i++) {
            if (!nodes[i].deleted && nodes[i].sequence.compareTo(nodes[i].rc.sequence) <= 0) {
                printLabel(out, nodes[i]);
            }
        }
        for (SingleNode i : nodes) {
            if (!i.deleted) {
                for (SingleNode j : i.neighbors) {
                    if (!j.deleted) {
                        printEdge(out, i, j);
                    }
                }
            }
        }
        out.flush();
        out.close();
    }

    private void printEdge(PrintWriter out, SingleNode first, SingleNode second) {
        StringBuilder edge = new StringBuilder("L\t");
        edge.append(getNodeId(first)).append('\t')
            .append((first.sequence.compareTo(first.rc.sequence) >= 0 ? "+" : "-")).append('\t')
            .append(getNodeId(second)).append('\t')
            .append((second.sequence.compareTo(second.rc.sequence) <= 0 ? "+" : "-")).append('\t')
            .append(k - 1).append("M");
        out.println(edge);
    }

    private String getNodeId(SingleNode node) {
        return (Math.min(node.rc.id, node.id) + 1) + "_i" + compID;
    }

    private void printLabel(PrintWriter out, SingleNode node) {
        out.print("S\t" + getNodeId(node) + "\t" + node.sequence);
        long coverage = 0;
        for (int i = 0; i + k <= node.sequence.length(); i++) {
            String kmer = node.sequence.substring(i, i + k);
            coverage += subgraph.get(normalizeDna(kmer));
        }
        coverage += subgraph.get(normalizeDna(node.sequence.substring(node.sequence.length()-k)))*(k-1);
        //coverage += coverage / (node.sequence.length() - k + 1) * (k-1);

        out.println("\tLN:i:" + (node.sequence.length()) + "\tKC:i:" + coverage);
    }

    public static String normalizeDna(String s) {
        String rc = reverseComplement(s);
        if (s.compareTo(rc) < 0) {
            return s;
        } else {
            return rc;
        }
    }
}

