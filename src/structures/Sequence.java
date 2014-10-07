package structures;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.io.writers.WritersUtils;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;

public class Sequence extends Dna {
    int avgWeight, minWeight, maxWeight;

    public Sequence(Dna dna, int avgWeight, int minWeight, int maxWeight) {
        super(dna);
        this.avgWeight = avgWeight;
        this.minWeight = minWeight;
        this.maxWeight = maxWeight;
    }


    public int averageWeight() {
        return avgWeight;
    }


    public static void printSequences(Iterable<Sequence> sequences, File file) throws IOException {
        LinkedList<String> comments = new LinkedList<String>();

        int sequenceId = 0;
        for (Sequence seq : sequences) {
            sequenceId++;
            String seqInfo = String.format("%d length=%d av_weight=%d min_weight=%d max_weight=%d",
                    sequenceId, seq.length(), seq.avgWeight, seq.minWeight, seq.maxWeight);
            comments.add(seqInfo);
        }
        WritersUtils.writeDnasToFastaFile(comments, sequences, file, false);
    }
}
