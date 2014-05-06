package algo;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import structures.ConnectedComponent;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ulyantsev on 06.05.14.
 *
 */
public class ComponentsBuilder {
    public static List<ConnectedComponent> splitStrategy(ArrayLong2IntHashMap hm,
                                                         int k,
                                                         int b1,
                                                         int b2,
                                                         String statFP,
                                                         Logger logger,
                                                         int availableProcessors) throws FileNotFoundException {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();

        PrintWriter statPW = new PrintWriter(statFP);
        for (int freqThreshold = 0; ; freqThreshold++) {
            List<ConnectedComponent> components = getComponents(hm, k, freqThreshold, availableProcessors);
            if (components.size() == 0) {
                break;
            }
            for (ConnectedComponent comp : components) {
                if (comp.size() < b1) {
                    banComponent(hm, comp);
                } else if (comp.size() < b2) {
                    ans.add(comp);
                    statPW.println(comp.size() + " " + comp.getWeight() + " " + freqThreshold);
                    banComponent(hm, comp);
                }
            }
            logger.debug("Freq = " + freqThreshold + ", components count = " + ans.size());
        }
        statPW.close();

        return ans;
    }

    private static void banComponent(ArrayLong2IntHashMap hm, ConnectedComponent component) {
        for (long kmer : component.kmers) {
            hm.add(kmer, -hm.get(kmer));
        }
    }

    private static List<ConnectedComponent> getComponents(ArrayLong2IntHashMap hm,
                                                   int k,
                                                   int freqThreshold,
                                                   int availableProcessors) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();

        ArrayLong2IntHashMap processedKmers =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors) / Math.log(2)) + 4);
        for (int i = 0; i < hm.hm.length; ++i) {
            for (Long2IntMap.Entry entry : hm.hm[i].long2IntEntrySet()) {
                long kmer = entry.getLongKey();
                if (processedKmers.get(kmer) > 0) {
                    continue;
                }

                int value = entry.getIntValue();
                if (value > freqThreshold) {
                    ConnectedComponent comp = getComponent(hm, k, kmer, freqThreshold, processedKmers);
                    ans.add(comp);
                }
            }
        }

        return ans;
    }

    private static ConnectedComponent getComponent(ArrayLong2IntHashMap hm,
                                                   int kValue,
                                                   long kmer,
                                                   int freqThreshold,
                                                   ArrayLong2IntHashMap processedKmers) {
        if (hm.get(kmer) <= freqThreshold || processedKmers.get(kmer) > 0) {
            return null;
        }
        ConnectedComponent ans = new ConnectedComponent();

        long weight = hm.get(kmer);
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue();

        queue.enqueue(kmer);
        processedKmers.add(kmer, 1);

        while (queue.size() > 0) {
            long kmerRepr = queue.dequeue();
            ans.add(kmerRepr);

            for (long neighbour : KmerOperations.possibleNeighbours(kmerRepr, kValue)) {
                int value = hm.get(neighbour);
                if (value <= freqThreshold || processedKmers.get(neighbour) > 0) {
                    continue;
                }
                weight += value;
                processedKmers.add(neighbour, 1);
                queue.enqueue(neighbour);
            }
        }

        ans.setWeight(weight);
        return ans;
    }

}
