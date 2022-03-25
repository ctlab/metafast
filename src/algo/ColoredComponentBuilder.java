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
import ru.ifmo.genetics.utils.pairs.MutablePair;
import ru.ifmo.genetics.utils.pairs.Pair;
import ru.ifmo.genetics.utils.tool.Tool;
import structures.ColoredKmers;
import structures.ConnectedComponent;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

import static io.IOUtils.withP;

public class ColoredComponentBuilder{
    public static List<ConnectedComponent> splitStrategy(BigLong2ShortHashMap hm, ColoredKmers coloredKmers,
                                                         int k, int b1, int b2,
                                                         String statFP, Logger logger,
                                                         int availableProcessors, boolean fakeIsCommon) throws FileNotFoundException {
        System.out.println(" " + k + " " +b1 +" " +  b2 +" " +  availableProcessors);
        ColoredComponentBuilder builder = new ColoredComponentBuilder(k, b1, b2, availableProcessors, statFP, logger);
        SPLIT_MODE mode =  (fakeIsCommon)  ? SPLIT_MODE.COMMON : SPLIT_MODE.SEPARATE;
        builder.run(hm, coloredKmers, mode);
        return builder.ans;
    }
    final private List<ConnectedComponent> ans;
    final private NonBlockingQueueExecutor executor;
    final int k;
    final int b1, b2;
    final String statFP;
    final private Logger logger;

    Comparator<Runnable> comparator = new Comparator<Runnable>() {
        @Override
        public int compare(Runnable o1, Runnable o2) {
            if (!(o1 instanceof ComponentsBuilder.Task) || !(o2 instanceof ComponentsBuilder.Task)) {
                return 0;
            }
            ComponentsBuilder.Task t1 = (ComponentsBuilder.Task) o1;
            ComponentsBuilder.Task t2 = (ComponentsBuilder.Task) o2;
            return -Long.compare(t1.component.size, t2.component.size);
        }
    };

    private ColoredComponentBuilder(int k, int b1, int b2, int availableProcessors, String statFP, Logger logger) {
        this.ans = new ArrayList<ConnectedComponent>();
        this.executor = new NonBlockingQueueExecutor(availableProcessors, comparator);
        this.k = k;
        this.b1 = b1;
        this.b2 = b2;
        this.statFP = statFP;
        this.logger = logger;
    }
    private static ConnectedComponent bfs(Long2ShortHashMapInterface hm, ColoredKmers coloredKmers,  long startKmer,
                                          LongArrayFIFOQueue queue,
                                          int k, int b2, int curFreqThreshold, SPLIT_MODE mode) {
        ConnectedComponent comp = new ConnectedComponent();
        comp.usedFreqThreshold = curFreqThreshold;

        queue.clear();

        queue.enqueue(startKmer);
        short value = hm.get(startKmer);
        int startcolor = coloredKmers.getColor(startKmer);
        int fakecolor = coloredKmers.colorsCNT;
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        comp.add(startKmer, value);
        while (queue.size() > 0) {
            long kmer = queue.dequeue();
            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if ((value > 0) && (coloredKmers.getColor(neighbour) == startcolor)) {    // i.e. if not precessed
                    queue.enqueue(neighbour);
                    hm.put(neighbour, (short) -value);
                    comp.add(neighbour, value);
                } else if (mode == SPLIT_MODE.COMMON && coloredKmers.getColor(neighbour)==fakecolor) {
                    queue.enqueue(neighbour);
                    comp.add(neighbour, value);
                }
            }
        }

        return comp;
    }

    public enum SPLIT_MODE {
        SEPARATE, COMMON
    }

    private static List<Pair<ConnectedComponent, Integer>> findAllComponents(BigLong2ShortHashMap hm, ColoredKmers coloredKmers, int k, int b2, int curFreqThreshold, SPLIT_MODE mode) {
        //todo add parallel

        System.out.println("find all components start: ");
        List<Pair<ConnectedComponent, Integer>> ans = new ArrayList<Pair<ConnectedComponent, Integer>>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size()/2));
        int colorsCNT = coloredKmers.colorsCNT;

        Iterator<MutableLongShortEntry> iterator = hm.entryIterator();
        int cnt = 0;
        while (iterator.hasNext()) {
            cnt+=1;
            MutableLongShortEntry startKmer = iterator.next();

            int component_color = coloredKmers.getColor(startKmer.getKey());
            if (mode==SPLIT_MODE.COMMON && component_color == colorsCNT) {
                continue;
            }
            if (startKmer.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponent comp = bfs(hm, coloredKmers, startKmer.getKey(), queue, k, b2, curFreqThreshold, mode);
                ans.add(new MutablePair<ConnectedComponent, Integer>(comp, component_color));
            }
        }

        return ans;
    }

    private void run(BigLong2ShortHashMap hm, ColoredKmers coloredKmers, SPLIT_MODE mode) throws FileNotFoundException {
        Tool.info(logger, "First iteration...");
        Timer t = new Timer();

        long hmSize = hm.size();
        int curFreqThreshold = 1;  // current component is formed of k-mers with frequency >= 1
        List<Pair<ConnectedComponent, Integer>> newComps = findAllComponents(hm, coloredKmers, k, b2, curFreqThreshold, mode);

        int colorsCNT = coloredKmers.colorsCNT + 1;
        int[] smallBYC = new int[colorsCNT], okBYC = new int[colorsCNT], bigBYC = new int[colorsCNT];
        long[] smallKBYC = new long[colorsCNT], okKBYC = new long[colorsCNT];
        for (int i = 0; i<colorsCNT; i++) {
            smallBYC[i] = 0;
            okBYC[i] = 0;
            bigBYC[i] = 0;
            smallKBYC[i] = 0;
            okKBYC[i] = 0;
        }

        int small = 0, ok = 0, big = 0;
        long smallK = 0, okK = 0;

        for (Pair<ConnectedComponent, Integer> compWithColor : newComps) {
            ConnectedComponent comp = compWithColor.first();
            Integer color = compWithColor.second();

            if (comp.size < b1) {
                small++;
                smallBYC[color]++;
                smallK += comp.size;
                smallKBYC[color] += comp.size;
            } else if (comp.size <= b2) {
                ok++;
                ans.add(comp);
                okBYC[color]++;
                okK += comp.size;
                okKBYC[color]+=1;
            } else {
                big++;
                bigBYC[color]++;
            }
        }
        int ansFirst = ans.size();
        System.out.println("total size: " + ans.size());
        Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " good components, " +
                "and " + NumUtils.groupDigits(big) + " big ones");
        Tool.info(logger, "Found " + NumUtils.groupDigits(small) + " small components, ");
        Tool.info(logger, "By colors found: ");
        for (int i = 0; i<colorsCNT; i++) {
            Tool.info(logger, "#" + i + ": " + NumUtils.groupDigits(okBYC[i]) + " good components, " +
                    "and " + NumUtils.groupDigits(bigBYC[i]) + " big ones");
            Tool.info(logger, "Found " + NumUtils.groupDigits(smallBYC[i]) + " small components, ");
        }

        Tool.info(logger, "First iteration was finished in " + t);

        Tool.debug(logger, "Total components found = " + NumUtils.groupDigits(newComps.size()) + ", " +
                "kmers = " + NumUtils.groupDigits(hmSize));
        Tool.debug(logger, "Components count: small = " + withP(small, newComps.size()) + ", " +
                "ok = " + withP(ok, newComps.size()) + ", " +
                "big = " + withP(big, newComps.size()));
        Tool.debug(logger, "Components kmers: small = " + withP(smallK, hmSize) + ", " +
                "ok = " + withP(okK, hmSize) + ", " +
                "big = " + withP(hmSize - smallK - okK, hmSize));

        Tool.debug(logger, "Total components found  by colors: ");
        for (int i = 0; i<colorsCNT; i++) {
            Tool.debug(logger, "Components count: small = " + withP(smallBYC[i], newComps.size()) + ", " +
                    "ok = " + withP(okBYC[i], newComps.size()) + ", " +
                    "big = " + withP(bigBYC[i], newComps.size()));
            Tool.debug(logger, "Components kmers: small = " + withP(smallKBYC[i], hmSize) + ", " +
                    "ok = " + withP(okKBYC[i], hmSize) + ", " +
                    //неправильно в случае по цветам
                    "big = " + withP(hmSize - smallKBYC[i] - okKBYC[i], hmSize));
        }

        Tool.debug(logger, "FreqThreshold = " + curFreqThreshold + ", " +
                "components added = " + ok + ", total components added = " + ans.size());

        Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());

        hm = null;  // for cleaning
        newComps = null;

        Tool.debug(logger, "Memory used after cleaning = " + Misc.usedMemoryAsString() + ", final time = " + t);


        //todo add split big and remove small

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
}
