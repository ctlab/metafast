package structures;


import ru.ifmo.genetics.dna.kmers.ShortKmer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

public class Sequence {
    String repr;

    long totalWeight, minWeight, maxWeight;

    public Sequence(String repr, long totalWeight, long minWeight, long maxWeight) {
        this.repr = repr;
        this.totalWeight = totalWeight;
        this.minWeight = minWeight;
        this.maxWeight = maxWeight;
    }

    public int length() {
        return repr.length();
    }

    @Override
    public String toString() {
        return repr;
    }

    public int averageWeight() {
        return (int) (totalWeight / length());
    }

    public ShortKmer startKmer(int k) {
        return new ShortKmer(repr.substring(0, k));
    }

    public ShortKmer endKmer(int k) {
        return new ShortKmer(repr.substring(repr.length() - k, repr.length()));
    }

    public static void printSequences(Iterable<Sequence> sequences, File destination) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(destination);
        int sequenceId = 0;
        for (Sequence seq : sequences) {
            sequenceId++;
            String seqInfo = String.format(">%d length=%d av_weight=%d min_weight=%d max_weight=%d",
                    sequenceId, seq.length(), seq.averageWeight(), seq.minWeight, seq.maxWeight);
            pw.println(seqInfo);
            pw.println(seq);
        }

    }
}
