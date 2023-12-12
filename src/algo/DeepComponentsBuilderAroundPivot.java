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
 * Created by -- on 15.10.2023.
 */
public class DeepComponentsBuilderAroundPivot {

    public static List<ConnectedComponent> splitStrategy(BigLong2ShortHashMap hm,
                                                         int k, int depth, BigLong2ShortHashMap pivot,
                                                         String statFP, Logger logger) throws FileNotFoundException {

        DeepComponentsBuilderAroundPivot builder = new DeepComponentsBuilderAroundPivot(k, depth, statFP, logger);
        builder.run(hm, pivot);
        return builder.ans;
    }

    final private List<ConnectedComponent> ans;
    final int k;
    final int depth;
    final String statFP;
    final private Logger logger;


    public DeepComponentsBuilderAroundPivot(int k, int depth, String statFP, Logger logger) {
        this.ans = new ArrayList<ConnectedComponent>();
        this.k = k;
        this.depth = depth;
        this.statFP = statFP;
        this.logger = logger;
    }

    private void run(BigLong2ShortHashMap hm, BigLong2ShortHashMap pivot) throws FileNotFoundException {
        Timer t = new Timer();

        // current component is formed of k-mers with frequency >= 1
        List<ConnectedComponent> newComps = findAllComponents(hm, k, depth, pivot);

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
        statPW.println("# component.no\tcomponent.size\tcomponent.weight\tcomponent.nPivotKmers\tusedFreqThreshold");
        for (int i = 0; i < ans.size(); i++) {
            ConnectedComponent comp = ans.get(i);
            statPW.println((i + 1) + "\t" + comp.size + "\t" + comp.weight + "\t" + comp.n_pivot + "\t" + comp.usedFreqThreshold);
        }
        statPW.close();
    }


    /**
     * Assuming running in one thread for current hm!
     */
    private static List<ConnectedComponent> findAllComponents(Long2ShortHashMapInterface hm,
                                                              int k, int depth, BigLong2ShortHashMap pivot) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));
        LongArrayFIFOQueue parent = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));

        Iterator<MutableLongShortEntry> iterator = pivot.entryIterator();
        while (iterator.hasNext()) {
            MutableLongShortEntry startKmer = iterator.next();
            if (startKmer.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponent comp = bfs(hm, startKmer.getKey(), queue, parent, k, depth, pivot);
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
                                          LongArrayFIFOQueue parent, int k, int depth, BigLong2ShortHashMap pivot) {
        ConnectedComponent comp = new ConnectedComponent();

        queue.clear();
        parent.clear();

        short value = hm.get(startKmer);
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        pivot.put(startKmer, (short) -value);
        comp.add(startKmer, value);
        comp.add_pivot();

        // extend to right
        {
            int n_neighbours = 0;
            List<Long> rightNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.rightNeighbours(startKmer, k)) {
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
                        comp.add_pivot();
                    }
                    comp.add(neighbour, value);
                }
                // if branching path, dfs into each branch to find another pivot or fail
                else {
                    for (long neighbour : rightNeighbours) {
                        List<Long> kmersOnPath = new ArrayList<Long>();
                        HashSet<Long> bestPath = new HashSet<Long>();

                        int n_piv = dfs(neighbour, startKmer, hm, pivot, k, bestPath, kmersOnPath, 0, depth, 0);
                        assert bestPath.size() <= depth;

                        if (n_piv > 0) {
                            value = hm.get(neighbour);
                            hm.put(neighbour, (short) -value);
                            value = pivot.get(neighbour);
                            if (value > 0) {
                                comp.add_pivot();
                                pivot.put(neighbour, (short) -value);
                            }
                            comp.add(neighbour, value);
                            int pathLength = kmersOnPath.size();
                            for (long foundKmer : kmersOnPath) {
                                value = hm.get(foundKmer);
                                comp.add(foundKmer, value);
                                hm.put(foundKmer, (short) -value);
                                value = pivot.get(foundKmer);
                                comp.add_pivot(n_piv);

                                if (value > 0) {
                                    pivot.put(foundKmer, (short) -value);
                                }
                            }
                            if (pathLength >= 2) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(kmersOnPath.get(pathLength - 2));
                            } else {
                                if (pathLength == 1) {
                                    queue.enqueue(kmersOnPath.get(0));
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
            for (long neighbour : KmerOperations.leftNeighbours(startKmer, k)) {
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
                        comp.add_pivot();
                    }
                    comp.add(neighbour, value);
                }
                // if branching path, dfs into each branch to find another pivot or fail
                else {
                    for (long neighbour : leftNeighbours) {
                        List<Long> kmersOnPath = new ArrayList<Long>();
                        HashSet<Long> bestPath = new HashSet<Long>();

                        int n_piv = dfs(neighbour, startKmer, hm, pivot, k, bestPath, kmersOnPath, 0, depth, 0);
                        assert bestPath.size() <= depth;

                        if (n_piv > 0) {
                            value = hm.get(neighbour);
                            hm.put(neighbour, (short) -value);
                            if (pivot.get(neighbour) > 0) {
                                comp.add_pivot();
                                pivot.put(neighbour, (short) -pivot.get(neighbour));
                            }
                            comp.add(neighbour, value);
                            int pathLength = kmersOnPath.size();
                            for (long foundKmer : kmersOnPath) {
                                value = hm.get(foundKmer);
                                comp.add(foundKmer, value);
                                hm.put(foundKmer, (short) -value);
                                comp.add_pivot(n_piv);

                                if (pivot.get(foundKmer) > 0) {
                                    pivot.put(foundKmer, (short) -pivot.get(foundKmer));
                                }
                            }
                            if (pathLength >= 2) {
                                queue.enqueue(kmersOnPath.get(pathLength - 1));
                                parent.enqueue(kmersOnPath.get(pathLength - 2));
                            } else {
                                if (pathLength == 1) {
                                    queue.enqueue(kmersOnPath.get(0));
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
            for (long neighbour : KmerOperations.rightNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    rightNeighbours.add(neighbour);
                    right_neighbours++;
                }
            }
            int left_neighbours = 0;
            List<Long> leftNeighbours = new ArrayList<Long>();
            for (long neighbour : KmerOperations.leftNeighbours(kmer, k)) {
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
                    break;
                }
            }
            for (long val : KmerOperations.rightNeighbours(kmer, k)) {
                if (val == prev) {
                    n_neighbours = left_neighbours;
                    neighbours = leftNeighbours;
                    break;
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
                    comp.add_pivot();
                }
                comp.add(neighbour, value);
            }
            // if branching path, dfs into each branch to find another pivot or fail
            else {
                for (long neighbour : neighbours) {
                    List<Long> kmersOnPath = new ArrayList<Long>();
                    HashSet<Long> bestPath = new HashSet<Long>();
                    int n_piv = dfs(neighbour, kmer, hm, pivot, k, bestPath, kmersOnPath, 0, depth, 0);
                    assert bestPath.size() <= depth;

                    if (n_piv > 0) {
                        value = hm.get(neighbour);
                        hm.put(neighbour, (short) -value);
                        if (pivot.get(neighbour) > 0) {
                            comp.add_pivot();
                            pivot.put(neighbour, (short) -pivot.get(neighbour));
                        }
                        comp.add(neighbour, value);
                        int pathLength = kmersOnPath.size();
                        for (long foundKmer : kmersOnPath) {
                            value = hm.get(foundKmer);
                            comp.add(foundKmer, value);
                            hm.put(foundKmer, (short) -value);
                            comp.add_pivot(n_piv);

                            if (pivot.get(foundKmer) > 0) {
                                pivot.put(foundKmer, (short) -pivot.get(foundKmer));
                            }
                        }
                        if (pathLength >= 2) {
                            queue.enqueue(kmersOnPath.get(pathLength - 1));
                            parent.enqueue(kmersOnPath.get(pathLength - 2));
                        } else {
                            if (pathLength == 1) {
                                queue.enqueue(kmersOnPath.get(0));
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


        return comp;
    }


    /**
     *
     * @param startKmer – k-mer to start DFS
     * @param parentKmer – previous for `startKmer`
     * @param hm – hashmap with all k-mers
     * @param pivot – hashmap with all pivot k-mers
     * @param k – k-mers length
     * @param kmersOnPath – hashset of all k-mers in current path
     * @param bestPath – list of vertices in current path
     * @param pivotKmersOnPath – amount of pivot k-mers on current path
     * @param depthAvailable – available depth in k-mers
     * @param globalBest – amount of pivot k-mers on the best path
     * @return amount of pivot k-mers on best selected path
     */
    private static int dfs(long startKmer, long parentKmer, Long2ShortHashMapInterface hm,
                           BigLong2ShortHashMap pivot, int k, HashSet<Long> kmersOnPath,
                           List<Long> bestPath, int pivotKmersOnPath, int depthAvailable, int globalBest) {
        if (depthAvailable == 0) {
            return pivotKmersOnPath;
        }

        List<Long> rightNeighbours = new ArrayList<Long>();
        for (long neighbour : KmerOperations.rightNeighbours(startKmer, k)) {
            short value = hm.get(neighbour);
            if (value > 0) {
                rightNeighbours.add(neighbour);
            }
        }
        List<Long> leftNeighbours = new ArrayList<Long>();
        for (long neighbour : KmerOperations.leftNeighbours(startKmer, k)) {
            short value = hm.get(neighbour);
            if (value > 0) {
                leftNeighbours.add(neighbour);
            }
        }

        List<Long> neighbours = null;
        for (long val : KmerOperations.leftNeighbours(startKmer, k)) {
            if (val == parentKmer) {
                neighbours = rightNeighbours;
                break;
            }
        }
        for (long val : KmerOperations.rightNeighbours(startKmer, k)) {
            if (val == parentKmer) {
                neighbours = leftNeighbours;
                break;
            }
        }

        assert neighbours != null;
        /*
        Explore each neighbour and update amount of pivot k-mers on each path.
        Update `bestPath` if it is better, than path of maximal depth.
         */
        boolean hasNeigbours = false;
        for (long neighbour : neighbours) {
            // make true copy
            HashSet<Long> kmersOnPathBack = new HashSet<Long>(kmersOnPath);

            // explore each neighbour which is available and not on current path
            short value = hm.get(neighbour);
            if (value > 0 && !kmersOnPath.contains(neighbour)) {
                kmersOnPathBack.add(neighbour);
                hasNeigbours = true;
            } else {
                continue;
            }

            int n_pivot;
            if (pivot.get(neighbour) > 0) {
                n_pivot = dfs(neighbour, startKmer, hm, pivot, k, kmersOnPathBack, bestPath, pivotKmersOnPath + 1, depthAvailable - 1, globalBest);
            } else {
                n_pivot = dfs(neighbour, startKmer, hm, pivot, k, kmersOnPathBack, bestPath, pivotKmersOnPath, depthAvailable - 1, globalBest);
            }

            // update best path only for max depth
            if (n_pivot > globalBest) {
                globalBest = n_pivot;
                if (depthAvailable == 1) {
                    bestPath.clear();
                    bestPath.addAll(kmersOnPathBack);
                }
            }
        }

        // update best path if no more neighbours
        if (!hasNeigbours && pivotKmersOnPath > globalBest) {
            globalBest = pivotKmersOnPath;
            bestPath.clear();
            bestPath.addAll(kmersOnPath);
        }

        return globalBest;
    }

}
