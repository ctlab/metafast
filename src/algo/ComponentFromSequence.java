package algo;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import structures.SequenceComponent;

import java.util.List;

/**
 * Created by -- on 03.02.2020.
 */
public class ComponentFromSequence implements Runnable {
    private final int k;
    private final Dna dna;
    private final List<SequenceComponent> components;

    public ComponentFromSequence(List<SequenceComponent> components, Dna dna, int k) {
        this.k = k;
        this.dna = dna;
        this.components = components;
    }

    @Override
    public void run() {
        SequenceComponent comp = new SequenceComponent();
        for (ShortKmer kmer : ShortKmer.kmersOf(dna, k)) {
            comp.add(kmer.toLong());
        }
        components.add(comp);
    }
}
