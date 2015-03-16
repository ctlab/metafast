package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.executors.NonBlockingQueueExecutor;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.Long2ShortHashMap;
import ru.ifmo.genetics.structures.map.Long2ShortHashMapInterface;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;
import structures.ConnectedComponent;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;

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
        logger.info("First iteration...");
        Timer t = new Timer();

        long hmSize = hm.size();
        int freqThreshold = 1;  // i.e. consider only kmers presented in reads
        List<ConnectedComponent> newComps = findAllComponents(hm, k, b2, freqThreshold+1);

        int small = 0, ok = 0, big = 0;
        long smallK = 0, okK = 0;


        List<ConnectedComponent> toProcess = new ArrayList<ConnectedComponent>();
        for (ConnectedComponent comp : newComps) {
            comp.usedFreqThreshold = freqThreshold;
            if (comp.size < b1) {
                small++;
                smallK += comp.size;
            } else if (comp.size < b2) {
                ok++;
                okK += comp.size;
                ans.add(comp);
            } else {
                big++;
                toProcess.add(comp);
            }
        }

        int ansFirst = ans.size();
        logger.info("Found " + NumUtils.groupDigits(ok) + " good components, " +
                "and " + NumUtils.groupDigits(big) + " big ones");
        logger.info("First iteration was finished in " + t);

        logger.debug("Total components found = " + NumUtils.groupDigits(newComps.size()) + ", " +
                "kmers = " + NumUtils.groupDigits(hmSize));
        logger.debug("Components count: small = " + withP(small, newComps.size()) + ", " +
                "ok = " + withP(ok, newComps.size()) + ", " +
                "big = " + withP(big, newComps.size()));
        logger.debug("Components kmers: small = " + withP(smallK, hmSize) + ", " +
                "ok = " + withP(okK, hmSize) + ", " +
                "big = " + withP(hmSize - smallK - okK, hmSize));
        logger.debug("FreqThreshold = " + freqThreshold + ", " +
                "components added = " + ok + ", total components added = " + ans.size());

        logger.debug("Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());

        hm = null;  // for cleaning
        newComps = null;

        logger.debug("Memory used after cleaning = " + Misc.usedMemoryAsString() + ", final time = " + t);


        if (big != 0) {
            logger.info("Following iterations...");
            t.start();

            for (ConnectedComponent comp : toProcess) {
                executor.addTask(new Task(comp));
            }
            toProcess = null;

            ConnectedComponent biggest = ((Task) executor.tasks.peek()).component;
            logger.debug("Biggest component has " +
                    withP(biggest.size, hmSize, "kmers", "of initial hm size"));
            logger.debug("Saved to new hm from it = " +
                    withP(biggest.hm.size(), biggest.size, "kmers", "of its size"));
            biggest = null;

            executor.startWorkers();

            executor.waitForTasksToFinish();
            executor.shutdownAndAwaitTermination();


            logger.info("Found " + NumUtils.groupDigits(ans.size() - ansFirst) + " good components " +
                    "extracted from big components");
            logger.info("All the following iterations were finished in " + t);
            logger.debug("Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                    "after it = " + Misc.usedMemoryAsString());
        }


        // post processing...
        Collections.sort(ans);
        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.size\tcomponent.weight\tusedFreqThreshold");
        for (ConnectedComponent comp : ans) {
            statPW.println(comp.size + "\t" + comp.weight + "\t" + comp.usedFreqThreshold);
        }
        statPW.close();
    }



    class Task implements Runnable {
        ConnectedComponent component;   // task is to split this component to good ones

        Task(ConnectedComponent component) { this.component = component; }
        @Override
        public void run() {
            int freqThreshold = component.usedFreqThreshold + 1;

            List<ConnectedComponent> newComps =
                    findAllComponents(component.hm, k, b2, freqThreshold+1);

            for (ConnectedComponent comp : newComps) {
                comp.usedFreqThreshold = freqThreshold;
                if (comp.size < b1) {
                    // skipping
                } else if (comp.size < b2) {
                    ans.add(comp);
                } else {
                    executor.addTask(new Task(comp));
                }
            }
        }
    }

    Comparator<Runnable> comparator = new Comparator<Runnable>() {
        @Override
        public int compare(Runnable o1, Runnable o2) {
            if (!(o1 instanceof Task) && !(o1 instanceof Task)) {
                return 0;
            }
            Task t1 = (Task) o1;
            Task t2 = (Task) o2;
            return -Long.compare(t1.component.size, t2.component.size);
        }
    };





    private static List<ConnectedComponent> findAllComponents(Long2ShortHashMapInterface hm,
                                                   int k, int b2, int newFreqThreshold) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size()/2));

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            if (entry.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponent comp = findComponent(hm, entry.getKey(), queue, k, b2, newFreqThreshold);
                ans.add(comp);
            }
        }

        return ans;
    }


    /**
     * Breadth-first search to make the traversal of the component.
     * If the component is small (less than b2 vertices), all its kmers is saved to
     * ConnectedComponent.kmers, else a subset of hm is stored to ConnectedComponent.hm structure.
     */
    private static ConnectedComponent findComponent(Long2ShortHashMapInterface hm, long startKmer,
                                             LongArrayFIFOQueue queue,
                                                   int k, int b2, int newFreqThreshold) {

        ConnectedComponent comp = new ConnectedComponent();

        queue.clear();
        boolean bigComp = false;

        queue.enqueue(startKmer);
        short value = hm.get(startKmer);
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        comp.add(startKmer, value);

        while (queue.size() > 0) {
            long kmer = queue.dequeue();

            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {    // i.e. if not precessed
                    queue.enqueue(neighbour);
                    hm.put(neighbour, (short) -value);
                    if (bigComp) {
                        if (value >= newFreqThreshold) {
                            comp.hm.put(neighbour, value);
                        }
                        comp.size++;
                    } else {
                        comp.add(neighbour, value);
                        if (comp.size == b2) {
                            // component is big one
                            bigComp = true;
                            comp.hm = new BigLong2ShortHashMap(4, 13);
                            for (long kk : comp.kmers) {
                                value = (short) -hm.get(kk);
                                assert value > 0;
                                if (value >= newFreqThreshold) {
                                    comp.hm.put(kk, value);
                                }
                            }
                            comp.kmers = null;
                        }
                    }
                }
            }
        }

        return comp;
    }

}
