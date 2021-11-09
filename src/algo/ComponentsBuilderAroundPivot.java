package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.Long2ShortHashMapInterface;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;
import structures.ConnectedComponent;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by -- on 05.02.2021.
 */
public class ComponentsBuilderAroundPivot {

    public static List<ConnectedComponent> splitStrategy(BigLong2ShortHashMap hm,
                                                         int k, BigLong2ShortHashMap pivot,
                                                         String statFP, Logger logger) throws FileNotFoundException {

        ComponentsBuilderAroundPivot builder = new ComponentsBuilderAroundPivot(k, statFP, logger);
        builder.run(hm, pivot);
        return builder.ans;
    }

    final private List<ConnectedComponent> ans;
    final int k;
    final String statFP;
    final private Logger logger;


    public ComponentsBuilderAroundPivot(int k, String statFP, Logger logger) {
        this.ans = new ArrayList<ConnectedComponent>();
        this.k = k;
        this.statFP = statFP;
        this.logger = logger;
    }

    private void run(BigLong2ShortHashMap hm, BigLong2ShortHashMap pivot) throws FileNotFoundException {
        Timer t = new Timer();

        // current component is formed of k-mers with frequency >= 1
        List<ConnectedComponent> newComps = findAllComponents(hm, k, pivot);

        int ok = 0;


        for (ConnectedComponent comp : newComps) {
            ok++;
            ans.add(comp);
        }

        Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " components");
        Tool.info(logger, "Iteration was finished in " + t);

        Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());

        hm = null;  // for cleaning
        newComps = null;
        Tool.debug(logger, "Memory used after cleaning = " + Misc.usedMemoryAsString() + ", final time = " + t);


        // post processing...
        Tool.debug(logger, "ans.size = " + ans.size());


        Collections.sort(ans);

        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.no\tcomponent.size\tcomponent.weight\tusedFreqThreshold");
        for (int i = 0; i < ans.size(); i++) {
            ConnectedComponent comp = ans.get(i);
            statPW.println((i + 1) + "\t" + comp.size + "\t" + comp.weight + "\t" + comp.usedFreqThreshold);
        }
        statPW.close();
    }


    /**
     * Assuming running in one thread for current hm!
     */
    private static List<ConnectedComponent> findAllComponents(Long2ShortHashMapInterface hm,
                                                              int k, BigLong2ShortHashMap pivot) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size()/2));
        LongArrayFIFOQueue parent = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size()/2));

        Iterator<MutableLongShortEntry> iterator = pivot.entryIterator();
        while (iterator.hasNext()) {
            MutableLongShortEntry startKmer = iterator.next();
            if (startKmer.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponent comp = bfs(hm, startKmer.getKey(), queue, parent, k, pivot);
                ans.add(comp);
            }
        }

        return ans;
    }


    /**
     * Breadth-first search to make the traversal of the component.
     * All its kmers are saved to ConnectedComponent.kmers.
     */
    private static ConnectedComponent bfs(Long2ShortHashMapInterface hm, long startKmer,
                                          LongArrayFIFOQueue queue,
                                          LongArrayFIFOQueue parent, int k, BigLong2ShortHashMap pivot) {
        ConnectedComponent comp = new ConnectedComponent();

        queue.clear();
        parent.clear();

        short value = hm.get(startKmer);
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        pivot.put(startKmer, (short) -value);
        comp.add(startKmer, value);

        // extend to right
        {
            int n_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour: KmerOperations.rightNeighbours(startKmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    n_neighbours++;
                }
            }
            if (n_neighbours == 0) {
                // do nothing
            } else {
                // if single path =>  extend
                if (n_neighbours == 1) {
                    long neighbour = rightNeighbours.get(0);
                    value = hm.get(neighbour);
                    queue.enqueue(neighbour);
                    parent.enqueue(startKmer);
                    hm.put(neighbour, (short) -value);
                    if (pivot.get(neighbour) > 0) {
                        pivot.put(neighbour, (short) -value);
                    }
                    comp.add(neighbour, value);
                }
                // if branching path, dfs into each branch to find another pivot or fail
                else {
                    for (long neighbour : rightNeighbours) {
                        List<Long> kmersOnPath = new ArrayList<Long>();
                        boolean goodPath = dfs(neighbour, startKmer, hm, pivot, k, kmersOnPath);
                        if (goodPath) {
                            value = hm.get(neighbour);
                            hm.put(neighbour, (short) -value);
                            if (pivot.get(neighbour) > 0) {
                                pivot.put(neighbour, (short) -value);
                            }
                            comp.add(neighbour, value);
                            int pathLength = kmersOnPath.size();
                            for (long foundKmer: kmersOnPath) {
                                value = hm.get(foundKmer);
                                comp.add(foundKmer, (short) -value);
                            }

                            if (kmersOnPath.size() >= 2) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(kmersOnPath.get(pathLength - 2));
                            }
                            else {
                                if (kmersOnPath.size() == 1) {
                                    queue.enqueue(kmersOnPath.get(pathLength - 1));
                                    parent.enqueue(neighbour);
                                } else {
                                    queue.enqueue(neighbour);
                                    parent.enqueue(startKmer);
                                }
                            }
                        }
                    }
                }
            }
        }

        // extend to left
        {
            int n_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour: KmerOperations.leftNeighbours(startKmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    leftNeighbours.add(neighbour);
                    n_neighbours++;
                }
            }
            if (n_neighbours == 0) {
                // do nothing
            } else {
                // if single path =>  extend
                if (n_neighbours == 1) {
                    long neighbour = leftNeighbours.get(0);
                    value = hm.get(neighbour);
                    queue.enqueue(neighbour);
                    parent.enqueue(startKmer);
                    hm.put(neighbour, (short) -value);
                    if (pivot.get(neighbour) > 0) {
                        pivot.put(neighbour, (short) -value);
                    }
                    comp.add(neighbour, value);
                }
                // if branching path, dfs into each branch to find another pivot or fail
                else {
                    for (long neighbour : leftNeighbours) {
                        List<Long> kmersOnPath = new ArrayList<Long>();
                        boolean goodPath = dfs(neighbour, startKmer, hm, pivot, k, kmersOnPath);
                        if (goodPath) {
                            value = hm.get(neighbour);
                            hm.put(neighbour, (short) -value);
                            if (pivot.get(neighbour) > 0) {
                                pivot.put(neighbour, (short) -value);
                            }
                            comp.add(neighbour, value);
                            int pathLength = kmersOnPath.size();
                            for (long foundKmer: kmersOnPath) {
                                value = hm.get(foundKmer);
                                comp.add(foundKmer, (short) -value);
                            }

                            if (kmersOnPath.size() >= 2) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(kmersOnPath.get(pathLength - 2));
                            }
                            else {
                                if (kmersOnPath.size() == 1) {
                                    queue.enqueue(kmersOnPath.get(pathLength - 1));
                                    parent.enqueue(neighbour);
                                } else {
                                    queue.enqueue(neighbour);
                                    parent.enqueue(startKmer);
                                }
                            }
                        }
                    }
                }
            }
        }


        while (queue.size() > 0) {
            long kmer = queue.dequeue();
            long prev = parent.dequeue();

            int right_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour: KmerOperations.rightNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    right_neighbours++;
                }
            }
            int left_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour: KmerOperations.leftNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    leftNeighbours.add(neighbour);
                    left_neighbours++;
                }
            }

            int n_neighbours = 0;
            List<Long> neighbours = null;
            for (long val : KmerOperations.leftNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = right_neighbours;
                    neighbours = rightNeighbours;
                }
            }
            for (long val : KmerOperations.rightNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = left_neighbours;
                    neighbours = leftNeighbours;
                }
            }

            if (n_neighbours == 0) {
                continue;
                // do nothing
            }
            // if single path =>  extend
            if (n_neighbours == 1) {
                long neighbour = neighbours.get(0);
                value = hm.get(neighbour);
                queue.enqueue(neighbour);
                parent.enqueue(kmer);
                hm.put(neighbour, (short) -value);
                if (pivot.get(neighbour) > 0) {
                    pivot.put(neighbour, (short) -value);
                }
                comp.add(neighbour, value);
            }
            // if branching path, dfs into each branch to find another pivot or fail
            else
            {
                for (long neighbour: neighbours) {
                    List<Long> kmersOnPath = new ArrayList<Long>();
                    boolean goodPath = dfs(neighbour, kmer, hm, pivot, k, kmersOnPath);
                    if (goodPath) {
                        value = hm.get(neighbour);
                        hm.put(neighbour, (short) -value);
                        if (pivot.get(neighbour) > 0) {
                            pivot.put(neighbour, (short) -value);
                        }
                        comp.add(neighbour, value);
                        int pathLength = kmersOnPath.size();
                        for (long foundKmer: kmersOnPath) {
                            value = hm.get(foundKmer);
                            comp.add(foundKmer, (short) -value);
                        }
                        if (kmersOnPath.size() >= 2) {
                            queue.enqueue(kmersOnPath.get(pathLength - 1));
                            parent.enqueue(kmersOnPath.get(pathLength - 2));
                        }
                        else {
                            if (kmersOnPath.size() == 1) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(neighbour);
                            } else {
                                queue.enqueue(neighbour);
                                parent.enqueue(kmer);
                            }
                        }
                    }
                }
            }
        }



        return comp;
    }


    private static boolean dfs(long startKmer, long parentKmer, Long2ShortHashMapInterface hm,
                               BigLong2ShortHashMap pivot, int k, List<Long> kmersOnPath) {
        boolean foundPivot = false;
        long kmer = startKmer;
        long prev = parentKmer;

        while (true) {
            int right_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour: KmerOperations.rightNeighbours(kmer, k)) {
                short value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    right_neighbours++;
                }
            }
            int left_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour: KmerOperations.leftNeighbours(kmer, k)) {
                short value = hm.get(neighbour);
                if (value > 0) {
                    leftNeighbours.add(neighbour);
                    left_neighbours++;
                }
            }

            int n_neighbours = 0;
            List<Long> neighbours = null;
            for (long val : KmerOperations.leftNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = right_neighbours;
                    neighbours = rightNeighbours;
                }
            }
            for (long val : KmerOperations.rightNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = left_neighbours;
                    neighbours = leftNeighbours;
                }
            }


            // if single path =>  extend
            if (n_neighbours == 1) {
                long neighbour = neighbours.get(0);
                kmersOnPath.add(neighbour);

                short value = hm.get(neighbour);
                hm.put(neighbour, (short) -value);
                if (pivot.get(neighbour) > 0) {
                    foundPivot = true;
                    pivot.put(neighbour, (short) -value);
                }
                prev = kmer;
                kmer = neighbour;
            }
            // if branching path or no path, stop and return
            else
            {
                break;
            }

        }

        return foundPivot;
    }


}
