package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.executors.NonBlockingQueueExecutor;
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

import static io.IOUtils.withP;

public class ComponentsBuilder {


    public static List<ConnectedComponent> splitStrategy(BigLong2ShortHashMap hm,
                                                         int k, int b1, int b2,
                                                         String statFP, Logger logger,
                                                         int availableProcessors) throws FileNotFoundException {

        ComponentsBuilder builder = new ComponentsBuilder(k, b1, b2, availableProcessors, statFP, logger);
        builder.run(hm);
        return builder.ans;
    }




    final private List<ConnectedComponent> ans;
    final private NonBlockingQueueExecutor executor;
    final int k;
    final int b1, b2;
    final String statFP;
    final private Logger logger;


    private ComponentsBuilder(int k, int b1, int b2, int availableProcessors, String statFP, Logger logger) {
        this.ans = new ArrayList<ConnectedComponent>();
        this.executor = new NonBlockingQueueExecutor(availableProcessors, comparator);
        this.k = k;
        this.b1 = b1;
        this.b2 = b2;
        this.statFP = statFP;
        this.logger = logger;
    }




    private void run(BigLong2ShortHashMap hm) throws FileNotFoundException {
        Tool.info(logger, "First iteration...");
        Timer t = new Timer();

        long hmSize = hm.size();
        int curFreqThreshold = 1;  // current component is formed of k-mers with frequency >= 1
        List<ConnectedComponent> newComps = findAllComponents(hm, k, b2, curFreqThreshold);

        int small = 0, ok = 0, big = 0;
        long smallK = 0, okK = 0;


        List<ConnectedComponent> toProcess = new ArrayList<ConnectedComponent>();
        for (ConnectedComponent comp : newComps) {
            if (comp.size < b1) {
                small++;
                smallK += comp.size;
            } else if (comp.size <= b2) {
                ok++;
                okK += comp.size;
                ans.add(comp);
            } else {
                big++;
                toProcess.add(comp);
            }
        }

        int ansFirst = ans.size();
        Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " good components, " +
                "and " + NumUtils.groupDigits(big) + " big ones");
        Tool.info(logger, "First iteration was finished in " + t);

        Tool.debug(logger, "Total components found = " + NumUtils.groupDigits(newComps.size()) + ", " +
                "kmers = " + NumUtils.groupDigits(hmSize));
        Tool.debug(logger, "Components count: small = " + withP(small, newComps.size()) + ", " +
                "ok = " + withP(ok, newComps.size()) + ", " +
                "big = " + withP(big, newComps.size()));
        Tool.debug(logger, "Components kmers: small = " + withP(smallK, hmSize) + ", " +
                "ok = " + withP(okK, hmSize) + ", " +
                "big = " + withP(hmSize - smallK - okK, hmSize));
        Tool.debug(logger, "FreqThreshold = " + curFreqThreshold + ", " +
                "components added = " + ok + ", total components added = " + ans.size());

        Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());

        hm = null;  // for cleaning
        newComps = null;

        Tool.debug(logger, "Memory used after cleaning = " + Misc.usedMemoryAsString() + ", final time = " + t);


        if (big != 0) {
            Tool.info(logger, "Following iterations...");
            t.start();

            for (ConnectedComponent comp : toProcess) {
                executor.addTask(new Task(comp));
            }
            toProcess = null;

            ConnectedComponent biggest = ((Task) executor.tasks.peek()).component;
            Tool.debug(logger, "Biggest component has " +
                    withP(biggest.size, hmSize, "kmers", "of initial hm size"));
            Tool.debug(logger, "Saved to new hm from it = " +
                    withP(biggest.nextHM.size(), biggest.size, "kmers", "of its size"));
            biggest = null;

            executor.startWorkers();

            executor.waitForTasksToFinish();
            executor.shutdownAndAwaitTermination();


            Tool.info(logger, "Found " + NumUtils.groupDigits(ans.size() - ansFirst) + " good components " +
                    "extracted from big components");
            Tool.info(logger, "All the following iterations were finished in " + t);
            Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                    "after it = " + Misc.usedMemoryAsString());
        }


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



    class Task implements Runnable {
        ConnectedComponent component;   // task is to split this component to good ones

        Task(ConnectedComponent component) { this.component = component; }
        @Override
        public void run() {
            int curFreqThreshold = component.usedFreqThreshold + 1;

            List<ConnectedComponent> newComps =
                    findAllComponents(component.nextHM, k, b2, curFreqThreshold);

            for (ConnectedComponent comp : newComps) {
                if (comp.size < b1) {
                    // skipping
                } else if (comp.size <= b2) {
                    synchronized (ans) {
                        ans.add(comp);
                    }
                } else {
                    executor.addTask(new Task(comp));
                }
            }
        }
    }

    Comparator<Runnable> comparator = new Comparator<Runnable>() {
        @Override
        public int compare(Runnable o1, Runnable o2) {
            if (!(o1 instanceof Task) || !(o2 instanceof Task)) {
                return 0;
            }
            Task t1 = (Task) o1;
            Task t2 = (Task) o2;
            return -Long.compare(t1.component.size, t2.component.size);
        }
    };


    /**
     * Assuming running in one thread for current hm!
     */
    private static List<ConnectedComponent> findAllComponents(Long2ShortHashMapInterface hm,
                                                   int k, int b2, int curFreqThreshold) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size()/2));

        Iterator<MutableLongShortEntry> iterator = hm.entryIterator();
        while (iterator.hasNext()) {
            MutableLongShortEntry startKmer = iterator.next();
            if (startKmer.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponent comp = bfs(hm, startKmer.getKey(), queue, k, b2, curFreqThreshold);
                ans.add(comp);
            }
        }

        return ans;
    }

    /**
     * Breadth-first search to make the traversal of the component.
     * If the component is small (less than b2 vertices), all its kmers is saved to
     * ConnectedComponent.kmers, else a subset of hm is stored to ConnectedComponent.nextHM structure.
     */
    private static ConnectedComponent bfs(Long2ShortHashMapInterface hm, long startKmer,
                                          LongArrayFIFOQueue queue,
                                          int k, int b2, int curFreqThreshold) {

        ConnectedComponent comp = new ConnectedComponent();
        comp.usedFreqThreshold = curFreqThreshold;

        queue.clear();

        queue.enqueue(startKmer);
        short value = hm.get(startKmer);
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        comp.add(startKmer, value);
        boolean alreadyBigComp = false;

        while (queue.size() > 0) {
            long kmer = queue.dequeue();

            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {    // i.e. if not precessed
                    queue.enqueue(neighbour);
                    hm.put(neighbour, (short) -value);

                    if (!alreadyBigComp) {
                        comp.add(neighbour, value);
                        if (comp.size > b2) {
                            alreadyBigComp = true;
                            comp.nextHM = new BigLong2ShortHashMap(4, 13);
                            for (long kk : comp.kmers) {
                                value = (short) -hm.get(kk);
                                assert value > 0;
                                if (value >= curFreqThreshold+1) {
                                    comp.nextHM.put(kk, value);
                                }
                            }
                            comp.kmers = null;
                        }
                    } else {
                        if (value >= curFreqThreshold+1) {  // for next HM
                            comp.nextHM.put(neighbour, value);
                        }
                        comp.size++;
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
