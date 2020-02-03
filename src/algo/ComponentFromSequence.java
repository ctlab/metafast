package algo;

import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import structures.ConnectedComponent;

import java.util.Iterator;
import java.util.List;

/**
 * Created by -- on 03.02.2020.
 */
public class ComponentFromSequence implements Runnable {
    private final int k;
    private final Dna dna;
    private final List<ConnectedComponent> components;

    public ComponentFromSequence(List<ConnectedComponent> components, Dna dna, int k) {
        this.k = k;
        this.dna = dna;
        this.components = components;
    }

    @Override
    public void run() {
        ConnectedComponent comp = new ConnectedComponent();
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(5, 12, true);
        for (ShortKmer kmer : ShortKmer.kmersOf(dna, k)) {
            hm.addAndBound(kmer.toLong(), (short) 1);
        }
        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            short value = entry.getValue();
            comp.add(key, value);
        }
        components.add(comp);
    }
}
