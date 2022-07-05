package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.Long2LongHashMapInterface;
import ru.ifmo.genetics.structures.map.MutableLongLongEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;
import structures.SequenceComponent;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by -- on 28.06.2022.
 */
public class ColoredComponentsBuilder {

    final private HashMap<Integer, List<SequenceComponent>> ans;
    final int k;
    final String statFP;
    final private Logger logger;
    final int n_groups;
    final boolean separate;
    final boolean linear;
    final int n_comps;
    final double perc;


    public ColoredComponentsBuilder(int k, String statFP, Logger logger, int n_groups, boolean separate, boolean linear, int n_comps, double perc) {
        this.ans = new HashMap<>();
        this.k = k;
        this.statFP = statFP;
        this.logger = logger;
        this.n_groups = n_groups;
        this.separate = separate;
        this.linear = linear;
        this.n_comps = n_comps;
        this.perc = perc;
    }


    public static HashMap<Integer, List<SequenceComponent>> splitStrategy(
                        BigLong2LongHashMap hm, int k, String statFP, Logger logger,
                        int n_groups, boolean separate, boolean linear, int n_comps, double perc) throws FileNotFoundException {
        ColoredComponentsBuilder builder = new ColoredComponentsBuilder(k, statFP, logger, n_groups, separate, linear, n_comps, perc);
        builder.run(hm);
        return builder.ans;
    }

    private void run(BigLong2LongHashMap hm) throws FileNotFoundException {
        Timer t = new Timer();

        HashMap<Integer, List<SequenceComponent>> comps = findAllComponents(hm);

        ans.putAll(comps);

        int ok = ans.values().stream().mapToInt(List::size).sum();

        Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " components");
        Tool.info(logger, "Iteration was finished in " + t);

        Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());


        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.no\tcomponent.size\tcomponent.weight\tcomponent.color");

        int total = 0;
        for (Map.Entry<Integer, List<SequenceComponent>> entry: ans.entrySet()) {
            for (SequenceComponent comp: entry.getValue()) {
                statPW.println((total + 1) + "\t" + comp.size + "\t" + comp.weight + "\t" + entry.getKey());
                total++;
            }
        }

        statPW.close();
    }

    private HashMap<Integer, List<SequenceComponent>> findAllComponents(BigLong2LongHashMap hm) {
        HashMap<Integer, List<SequenceComponent>> ans = new HashMap<>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));

        Iterator<MutableLongLongEntry> iterator = hm.entryIterator();

        ArrayList<Integer> compsPerGroup = new ArrayList<>();
        for (int group = 0; group < n_groups; group++) {
            compsPerGroup.add(0);
            ans.put(group, new ArrayList<>());
        }

        while (iterator.hasNext() && (n_comps == -1 || compsPerGroup.stream().mapToInt(Integer::intValue).sum() < n_groups*n_comps)) {
            MutableLongLongEntry startKmer = iterator.next();
            if (hm.get(startKmer.getKey()) > 0) {
                int color = ColoredKmerOperations.getColor(startKmer.getValue(), perc);
                if (color == -1) {
                    continue;
                }

                if (n_comps == -1 || compsPerGroup.get(color) < n_comps) {
                    SequenceComponent comp;

                    if (linear) {
                        comp = bfsLinear(hm, startKmer.getKey(), queue);
                    } else {
                        comp = bfs(hm, startKmer.getKey(), queue);
                    }

                    if (comp.size > 0) {
                        compsPerGroup.set(color, compsPerGroup.get(color) + 1);
                        ans.get(color).add(comp);
                    }
                }
            }

        }

        return ans;
    }

    private SequenceComponent bfsLinear(Long2LongHashMapInterface hm, long startKmer,
                                               LongArrayFIFOQueue queue) {
        SequenceComponent comp = new SequenceComponent();
        queue.clear();

        long value = hm.get(startKmer);
        int startColor = ColoredKmerOperations.getColor(value, perc);

        queue.enqueue(startKmer);
        hm.put(startKmer, -value);
        comp.add(startKmer);

        while (queue.size() > 0) {
            long kmer = queue.dequeue();

            long[] neighbours = Arrays.stream(KmerOperations.possibleNeighbours(kmer, k)).filter(x -> hm.get(x) > 0).toArray();
            if (neighbours.length > 1) {
                long bestNeighbour = 0;
                int maxGood = -1;
                for (long neighbour : neighbours) {
                    if (hm.get(neighbour) > 0) {
                        int goodOnPath = countColorOnPath(neighbour, kmer, startColor, hm, comp);
                        if (goodOnPath > maxGood) {
                            maxGood = goodOnPath;
                            bestNeighbour = neighbour;
                        }
                    }
                }
                if (maxGood > 0) {
                    List<Long> kmersOnPath = getKmersOnPath(bestNeighbour, kmer, hm, comp);
                    for (long v : kmersOnPath) {
                        value = hm.get(v);
                        int color = ColoredKmerOperations.getColor(value, perc);

                        if (value > 0) {
                            if (color == startColor) {
                                hm.put(v, -value);
                                comp.add(v);
                            } else if (color == -1 && !comp.contains(v)) {
                                comp.add(v);
                            }
                        }
                   }
                    queue.enqueue(kmersOnPath.get(kmersOnPath.size() - 1));
                }
            } else if (neighbours.length == 1) {
                long neighbour = neighbours[0];
                value = hm.get(neighbour);
                int color = ColoredKmerOperations.getColor(value, perc);

                if (value > 0) {
                    if (color == startColor) {
                        queue.enqueue(neighbour);
                        hm.put(neighbour, -value);
                        comp.add(neighbour);
                    } else if (color == -1 && !comp.contains(neighbour)) {
                        queue.enqueue(neighbour);
                        comp.add(neighbour);
                    }
                }
            }
        }

        return comp;
    }

    private int countColorOnPath(long startKmer, long prevKmer, int startColor, Long2LongHashMapInterface hm, SequenceComponent comp) {
        int count = 0;
        long curV = startKmer;
        long prevV = prevKmer;
        while (true) {
            long value = hm.get(curV);
            int color = ColoredKmerOperations.getColor(value, perc);
            if (value < 0) {
                return -1;
            }
            if (color == startColor) {
                count += 1;
            }

            List<Long> next = new ArrayList<>();
            for (long neighbour: KmerOperations.possibleNeighbours(curV, k)) {
                if (neighbour != prevV && hm.get(neighbour) > 0) {
                    next.add(neighbour);
                }
            }
            if (next.size() == 1) {
                prevV = curV;
                curV = next.get(0);
            } else {
                break;
            }
        }
        return count;
    }

    private List<Long> getKmersOnPath(long startKmer, long prevKmer, Long2LongHashMapInterface hm, SequenceComponent comp) {
        List<Long> kmers = new ArrayList<>();
        long curV = startKmer;
        long prevV = prevKmer;
        while (true) {
            long value = hm.get(curV);
            if (value < 0) {
                return null;
            }
            kmers.add(curV);

            List<Long> next = new ArrayList<>();
            for (long neighbour: KmerOperations.possibleNeighbours(curV, k)) {
                if (neighbour != prevV && hm.get(neighbour) > 0) {
                    next.add(neighbour);
                }
            }
            if (next.size() == 1) {
                prevV = curV;
                curV = next.get(0);
            } else {
                break;
            }
        }
        return kmers;
    }


    private SequenceComponent bfs(Long2LongHashMapInterface hm, long startKmer, LongArrayFIFOQueue queue) {
        SequenceComponent comp = new SequenceComponent();
        queue.clear();

        long value = hm.get(startKmer);
        int startColor = ColoredKmerOperations.getColor(value, perc);

        queue.enqueue(startKmer);
        hm.put(startKmer, -value);
        comp.add(startKmer);

        while (queue.size() > 0) {
            long kmer = queue.dequeue();

            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                int color = ColoredKmerOperations.getColor(value, perc);
                if (value > 0) {
                    if (color == startColor) {
                        queue.enqueue(neighbour);
                        hm.put(neighbour, -value);
                        comp.add(neighbour);
                    } else if (!separate && color == -1 && !comp.contains(neighbour)) {
                        queue.enqueue(neighbour);
                        comp.add(neighbour);
                    }
                }
            }
        }
        return comp;
    }

}
