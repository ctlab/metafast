package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.statistics.QuantitativeStatistics;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.structures.set.BigLongHashSet;
import structures.ConnectedComponent;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class ComponentsBuilder {
    public static List<ConnectedComponent> splitStrategy(BigLong2ShortHashMap hm,
                                                         int k,
                                                         int b1, int b2,
                                                         String statFP,
                                                         Logger logger)
            throws FileNotFoundException {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();
        BigLongHashSet processedKmers = new BigLongHashSet(hm.capacity());

        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.size\tcomponent.weight\tfreqThreshold");
        for (int freqThreshold = 0; ; freqThreshold++) {
            Timer t = new Timer();
            List<ConnectedComponent> components = getComponents(hm, k, freqThreshold, b2, processedKmers);
            if (components.size() == 0) {
                break;
            }
            logger.debug("Time to find components = " + t);
            int added = 0;
            for (ConnectedComponent comp : components) {
                if (comp.size() < b1) {
                    banComponent(hm, comp);
                } else if (comp.size() < b2) {
                    ans.add(comp);
                    added++;
                    statPW.println(comp.size() + "\t" + comp.getWeight() + "\t" + freqThreshold);
                    banComponent(hm, comp);
                }
            }
            logger.debug("FreqThreshold = " + freqThreshold + ", " +
                    "components added = " + added + ", total components = " + ans.size());
        }
        statPW.close();

        return ans;
    }

    private static void banComponent(BigLong2ShortHashMap hm, ConnectedComponent component) {
        for (long kmer : component.kmers) {
            hm.put(kmer, (short) 0);
        }
    }

    private static List<ConnectedComponent> getComponents(BigLong2ShortHashMap hm,
                                                   int k, int freqThreshold, int b2,
                                                   BigLongHashSet processedKmers) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();

        processedKmers.reset();

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long kmer = entry.getKey();
            short value = entry.getValue();
            if (value > freqThreshold && !processedKmers.contains(kmer)) {
                ConnectedComponent comp = getComponent(hm, k, kmer, freqThreshold, b2, processedKmers);
                ans.add(comp);
            }
        }

        return ans;
    }

    private static LongArrayFIFOQueue queue = new LongArrayFIFOQueue(1 << 16);

    private static ConnectedComponent getComponent(BigLong2ShortHashMap hm,
                                                   int k,
                                                   long startKmer,
                                                   int freqThreshold, int b2,
                                                   BigLongHashSet processedKmers) {

        ConnectedComponent ans = new ConnectedComponent();
        long weight = 0;

        queue.clear();

        queue.enqueue(startKmer);
        processedKmers.add(startKmer);
        ans.add(startKmer);
        weight += hm.get(startKmer);

        while (queue.size() > 0) {
            long kmer = queue.dequeue();

            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                short value = hm.get(neighbour);
                if (value > freqThreshold && !processedKmers.contains(neighbour)) {
                    queue.enqueue(neighbour);
                    processedKmers.add(neighbour);
                    if (ans.size() < b2) {
                        ans.add(neighbour);
                        weight += value;
                    }
                }
            }
        }

        ans.setWeight(weight);
        return ans;
    }

}
