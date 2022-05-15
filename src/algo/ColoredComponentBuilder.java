package algo;

import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import org.apache.log4j.Logger;
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
import structures.ConnectedSetComponent;
import tools.ComponentColoredCutter;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

import static io.IOUtils.withP;
import static tools.ComponentColoredCutter.*;

public class ColoredComponentBuilder {
    public static List<ConnectedComponent> splitStrategy(BigLong2ShortHashMap hm, ColoredKmers coloredKmers,
                                                         int k, int b1, int b2,
                                                         String statFP, Logger logger,
                                                         int availableProcessors, SPLIT_MODE mode, START_KMER_MODE startMode, BFS_MODE bfsMode, int forEachColorCNT, double minForGreedStart, COMPONENT_SIZES_MODE res_mode) throws FileNotFoundException {
        System.out.println(" " + k + " " + b1 + " " + b2 + " " + availableProcessors);
        ColoredComponentBuilder builder = new ColoredComponentBuilder(k, b1, b2, statFP, logger, forEachColorCNT);
        builder.run(hm, coloredKmers, mode, startMode, bfsMode, minForGreedStart, res_mode);
        return builder.ans;
    }

    final private List<ConnectedComponent> ans;
    final int k;
    final int b1, b2;
    final String statFP;
    final private Logger logger;
    final int forEachColorCNT;


    private ColoredComponentBuilder(int k, int b1, int b2, String statFP, Logger logger, Integer forEachColorCNT) {
        this.ans = new ArrayList<>();
        this.k = k;
        this.b1 = b1;
        this.b2 = b2;
        this.statFP = statFP;
        this.logger = logger;
        this.forEachColorCNT = forEachColorCNT;
    }

    private static boolean kmerInHM(long kmer, Long2ShortHashMapInterface hm) {
        long value = hm.get(kmer);
        return (value != 0 && value != -1);
    }

    private static Pair<Integer, Long> numOfPaths(long kmer, Long2ShortHashMapInterface hm, int k) {
        int cnt = 0;
        long nn = 0;
        for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
            if (kmerInHM(neighbour, hm)) {
                cnt += 1;
                nn = neighbour;
            }
        }
        return new MutablePair<>(cnt, nn);
    }

    private static List<Long> realNeighbours(long kmer, Long2ShortHashMapInterface hm, int k) {
        List<Long> res = new ArrayList<>();
        for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
            if (kmerInHM(neighbour, hm)) res.add(neighbour);
        }
        return res;
    }

    private static void updateComp(Long2ShortHashMapInterface hm, LongArrayFIFOQueue queue, ConnectedComponent comp, long neighbour, short value) {
        queue.enqueue(neighbour);
        hm.put(neighbour, (short) -value);
        comp.add(neighbour, value);
    }

    private static List<Long> getKmersOnPath(long neighbour, int startColour, Long2ShortHashMapInterface hm, int k, ColoredKmers coloredKmers, int fakeColor, SPLIT_MODE split_mode, ConnectedSetComponent comp) {

        List<Long> res = new ArrayList<>();
        long curv = neighbour;
        while (true) {
            int v = hm.get(curv);
            int col = coloredKmers.getColor(curv);
            if (v < 0) {
                if (col==fakeColor && split_mode == SPLIT_MODE.COMMON) {
                    if (comp.contains(neighbour)) {
                        break;
                    }
                } else {
                    break;
                }
            }

            res.add(curv);
            Pair<Integer, Long> cntnn = numOfPaths(curv, hm, k);
            int cntn = cntnn.first();
            if (cntn == 1) {
                curv = cntnn.second();
            } else {
                break;
            }
        }
        return res;
    }

    private static int getColorCNTBeforePathSplit(long neighbour, int startColour, Long2ShortHashMapInterface hm, int k, ColoredKmers coloredKmers, int fakeColor, SPLIT_MODE split_mode, ConnectedSetComponent comp) {
        int res = 0;
        long curv = neighbour;
        while (true) {
            int v = hm.get(curv);
            int col = coloredKmers.getColor(curv);
            if (v < 0) {
                if (col==fakeColor && split_mode == SPLIT_MODE.COMMON) {
                    if (comp.contains(neighbour)) {
                        break;
                    }
                } else {
                    break;
                }
            }
            if (col == startColour) {
                res += 1;
            }
            Pair<Integer, Long> cntnn = numOfPaths(curv, hm, k);
            int cntn = cntnn.first();
            if (cntn == 1) {
                curv = cntnn.second();
            } else {
                break;
            }
        }
        return res;
    }

    public static ConnectedComponent bfsDeep(Long2ShortHashMapInterface hm, ColoredKmers coloredKmers, long startKmer,
                                             LongArrayFIFOQueue queue,
                                             int k, int curFreqThreshold, SPLIT_MODE splitMode) {
        ConnectedSetComponent comp = new ConnectedSetComponent();
        comp.usedFreqThreshold = curFreqThreshold;

        queue.clear();

        queue.enqueue(startKmer);
        short value = hm.get(startKmer);
        int startColour = coloredKmers.getColor(startKmer);
        int fakeColor = coloredKmers.colorsCNT;
        assert value > 0;

        hm.put(startKmer, (short) -value);  // removing
        comp.add(startKmer, value);
        while (queue.size() > 0) {
            long kmer = queue.dequeue();

            List<Long> neighbours = realNeighbours(kmer, hm, k);


            int cntn = neighbours.size();
            if (cntn > 1) {
                long bestNeighbour = 0;
                int maxvalue = -1;
                for (long neighbour : neighbours) {
                    if (hm.get(neighbour) > 0) {
                        int curcnt = getColorCNTBeforePathSplit(neighbour, startColour, hm, k, coloredKmers, fakeColor, splitMode, comp);
                        if (curcnt > maxvalue) {
                            maxvalue = curcnt;
                            bestNeighbour = neighbour;
                        }
                    }
                }
                if (maxvalue > 0) {
                    List<Long> onPath = getKmersOnPath(bestNeighbour, startColour, hm, k, coloredKmers, fakeColor, splitMode, comp);
                    for (long nn : onPath) {
                        short v = hm.get(nn);
                        hm.put(nn, (short) -v);
                        comp.add(nn, v);
                    }
                    queue.enqueue(onPath.get(onPath.size() - 1));
                }
            } else if (cntn == 1) {
                long nn = neighbours.get(0);
                int col = coloredKmers.getColor(nn);
                if (col == startColour) {
                    int v = hm.get(nn);
                    if (v > 0) {
                        updateComp(hm, queue, comp, nn, (short) v);
                    }
                }
            }
        }
        return comp;
    }

    private static ConnectedComponent bfsBest(Long2ShortHashMapInterface hm, ColoredKmers coloredKmers, long startKmer,
                                              LongArrayFIFOQueue queue,
                                              int k, int curFreqThreshold, SPLIT_MODE mode) {
        ConnectedComponent comp = new ConnectedComponent();
        comp.usedFreqThreshold = curFreqThreshold;

        queue.clear();

        queue.enqueue(startKmer);
        short value = hm.get(startKmer);
        int startColor = coloredKmers.getColor(startKmer);
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
        comp.add(startKmer, value);

        while (queue.size() > 0) {
            long kmer = queue.dequeue();
            long bestNeighbour = 0;
            double maxval = -1;
            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if (value > 0) {
                    double curp = coloredKmers.getColorDouble(kmer);
                    int curcol = coloredKmers.getColorInt(kmer); //without cut-off
                    if (curp > maxval) {
                        if (curcol == startColor) {
                            if (curp > 0.75) {
                                updateComp(hm, queue, comp, neighbour, value);
                            } else {
                                maxval = curp;
                                bestNeighbour = neighbour;
                            }
                        }
                    }
                }
            }
            if ((maxval != -1) && (mode == SPLIT_MODE.COMMON)) {
                updateComp(hm, queue, comp, bestNeighbour, hm.get(bestNeighbour));
            }

        }
        return comp;
    }

    private static ConnectedComponent bfs(Long2ShortHashMapInterface hm, ColoredKmers coloredKmers, long startKmer,
                                          LongArrayFIFOQueue queue,
                                          int k, int curFreqThreshold, SPLIT_MODE mode) {
        ConnectedSetComponent comp = new ConnectedSetComponent();
        comp.usedFreqThreshold = curFreqThreshold;

        queue.clear();

        queue.enqueue(startKmer);
        short value = hm.get(startKmer);
        int startColor = coloredKmers.getColor(startKmer);
        int fakeColor = coloredKmers.colorsCNT;
        assert value > 0;
        hm.put(startKmer, (short) -value);  // removing
//        BigLongHashSet used = new BigLongHashSet(1000);
        comp.add(startKmer, value);
        while (queue.size() > 0) {
            long kmer = queue.dequeue();
            //check neighbours if some of them are fake chose only one with best value (maybe need to see deeper)
            for (long neighbour : KmerOperations.possibleNeighbours(kmer, k)) {
                value = hm.get(neighbour);
                if ((value > 0) && (coloredKmers.getColor(neighbour) == startColor)) {    // i.e. if not precessed
                    updateComp(hm, queue, comp, neighbour, value);
                } else if ((mode == SPLIT_MODE.COMMON) && (coloredKmers.getColor(neighbour) == fakeColor) && (value > 0)) {
                    if (!comp.contains(neighbour)) {
                        queue.enqueue(neighbour);
                        comp.add(neighbour, value);
                    }
                }
            }
        }
        return comp;
    }


    private static int getStartKmer(BigLong2ShortHashMap hm, int prevind, List<Long> colorsQueue, int maxvalue) {
        int res = prevind;

        while (true) {
            res += 1;
            if (res == maxvalue) {
                return -1;
            }
            if (hm.get(colorsQueue.get(res)) > 0) {
                return res;
            }
        }
    }


    private static List<Pair<ConnectedComponent, Integer>> findComponentsGreed(BigLong2ShortHashMap hm, ColoredKmers coloredKmers, int k, int curFreqThreshold, SPLIT_MODE mode, int forEachColorMax, BFS_MODE bfsMode, double minForStart) {
        System.out.println("find components greed start: ");
        List<Pair<ConnectedComponent, Integer>> ans = new ArrayList<>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));
        int colorsCNT = coloredKmers.colorsCNT;
        int fakeColor = coloredKmers.colorsCNT;

        Iterator<MutableLongShortEntry> iterator = hm.entryIterator();

        ArrayList<Integer> forEachColorAdded = new ArrayList<>();
        for (int col = 0; col < colorsCNT; col++) {
            forEachColorAdded.add(0);
        }
        //todo make iterator random???
        while (iterator.hasNext()) {
            MutableLongShortEntry startKmer = iterator.next();

            int component_color = coloredKmers.getColor(startKmer.getKey());
            double kmerP = coloredKmers.getColorDouble(startKmer.getKey());
            if (component_color == fakeColor) {
                continue;
            }
            if ((forEachColorAdded.get(component_color) < forEachColorMax) && (kmerP > minForStart)) {
                forEachColorAdded.set(component_color, forEachColorAdded.get(component_color) + 1);
            } else {
                continue;
            }
            ConnectedComponent comp = null;
            switch (bfsMode) {
                case BEST:
                    comp = bfsBest(hm, coloredKmers, startKmer.getKey(), queue, k, curFreqThreshold, mode);
                    break;
                case ALL:
                    comp = bfs(hm, coloredKmers, startKmer.getKey(), queue, k, curFreqThreshold, mode);
                    break;
                case DEEP:
                    comp = bfsDeep(hm, coloredKmers, startKmer.getKey(), queue, k, curFreqThreshold, mode);
                    break;
            }
            if (comp != null) {
                ans.add(new MutablePair<>(comp, component_color));
            }
        }

        return ans;
    }

    //todo  проверить на других данных, чтобы вероятности были не 0, 0.5, 1
    private static List<Pair<ConnectedComponent, Integer>> findBestComponents(BigLong2ShortHashMap hm, ColoredKmers coloredKmers, int k, int curFreqThreshold, SPLIT_MODE mode, int forEachColor, BFS_MODE bfsMode) {
        System.out.println("find best components start: ");
        List<Pair<ConnectedComponent, Integer>> ans = new ArrayList<>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));
        int colorsCNT = coloredKmers.colorsCNT;
        ArrayList<List<Long>> startQueues = new ArrayList<>();
        int mincnt = forEachColor;
        ArrayList<Integer> iterators = new ArrayList<>();
        for (int i = 0; i < colorsCNT; i++) {
            iterators.add(0);
            List<Long> cc = coloredKmers.getValuesForColor(i, forEachColor);

            if (cc.size() < mincnt) {
                mincnt = cc.size();
            }
            startQueues.add(cc);
        }
        boolean finished = false;
        while (!finished) {
            for (int col = 0; col < colorsCNT; col++) {
                int ind = getStartKmer(hm, iterators.get(col), startQueues.get(col), mincnt);
                if (ind == -1) {
                    finished = true;
                    continue;
                }
                long startKmer = startQueues.get(col).get(ind);
                // i.e. if not precessed
                ConnectedComponent comp = null;
                switch (bfsMode) {
                    case BEST:
                        comp = bfsBest(hm, coloredKmers, startKmer, queue, k, curFreqThreshold, mode);
                        break;
                    case ALL:
                        comp = bfs(hm, coloredKmers, startKmer, queue, k, curFreqThreshold, mode);
                        break;
                    case DEEP:
                        comp = bfsDeep(hm, coloredKmers, startKmer, queue, k, curFreqThreshold, mode);
                        break;
                }
                if (comp != null) {
                    ans.add(new MutablePair<>(comp, col));
                }
            }
        }
        return ans;
    }


    private static List<Pair<ConnectedComponent, Integer>> findAllComponents(BigLong2ShortHashMap hm, ColoredKmers coloredKmers, int k, int curFreqThreshold, SPLIT_MODE mode, BFS_MODE bfsMode) {
        System.out.println("find all components start: ");
        List<Pair<ConnectedComponent, Integer>> ans = new ArrayList<>();
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue((int) Math.min(1 << 16, hm.size() / 2));
        int colorsCNT = coloredKmers.colorsCNT;

        Iterator<MutableLongShortEntry> iterator = hm.entryIterator();
        while (iterator.hasNext()) {
            MutableLongShortEntry startKmer = iterator.next();

            int component_color = coloredKmers.getColor(startKmer.getKey());
            if (mode == SPLIT_MODE.COMMON && component_color == colorsCNT) {
                continue;
            }
            if (startKmer.getValue() > 0) {    // i.e. if not precessed
                ConnectedComponent comp = null;
                switch (bfsMode) {
                    case BEST:
                        comp = bfsBest(hm, coloredKmers, startKmer.getKey(), queue, k, curFreqThreshold, mode);
                        break;
                    case ALL:
                        comp = bfs(hm, coloredKmers, startKmer.getKey(), queue, k, curFreqThreshold, mode);
                        break;
                    case DEEP:
                        comp = bfsDeep(hm, coloredKmers, startKmer.getKey(), queue, k, curFreqThreshold, mode);
                        break;
                }
                if (comp != null) {
                    ans.add(new MutablePair<>(comp, component_color));
                }
            }
        }

        return ans;
    }

    private void run(BigLong2ShortHashMap hm, ColoredKmers coloredKmers, SPLIT_MODE mode, ComponentColoredCutter.START_KMER_MODE startMode, ComponentColoredCutter.BFS_MODE bfsMode, double minForGreedStart, COMPONENT_SIZES_MODE res_mode) throws FileNotFoundException {
        Tool.info(logger, "First iteration...");
        Timer t = new Timer();

        long hmSize = hm.size();
        int curFreqThreshold = 1;  // current component is formed of k-mers with frequency >= 1

        List<Pair<ConnectedComponent, Integer>> newComps = null;
        switch (startMode) {
            case RANDOM:
                newComps = findAllComponents(hm, coloredKmers, k, curFreqThreshold, mode, bfsMode);
                break;
            case BEST:
                newComps = findBestComponents(hm, coloredKmers, k, curFreqThreshold, mode, this.forEachColorCNT, bfsMode);
                break;
            case GREED:
                newComps = findComponentsGreed(hm, coloredKmers, k, curFreqThreshold, mode, this.forEachColorCNT, bfsMode, minForGreedStart);
                break;
        }
        if (newComps == null) {
            newComps = new ArrayList<>();
        }

        int colorsCNT = coloredKmers.colorsCNT + 1;
        int[] smallBYC = new int[colorsCNT], okBYC = new int[colorsCNT], bigBYC = new int[colorsCNT];
        long[] smallKBYC = new long[colorsCNT], okKBYC = new long[colorsCNT], bigKbYC = new long[colorsCNT];
        for (int i = 0; i < colorsCNT; i++) {
            smallBYC[i] = 0;
            okBYC[i] = 0;
            bigBYC[i] = 0;
            smallKBYC[i] = 0;
            okKBYC[i] = 0;
            bigKbYC[i] = 0;
        }

        int small = 0, ok = 0, big = 0;
        long smallK = 0, okK = 0;

        for (Pair<ConnectedComponent, Integer> compWithColor : newComps) {
            ConnectedComponent comp = compWithColor.first();
            Integer color = compWithColor.second();

            if (comp.size < b1) {
                if (res_mode==COMPONENT_SIZES_MODE.ALL || res_mode == COMPONENT_SIZES_MODE.SMALL) {
                    ans.add(comp);
                }
//
                small++;
                smallBYC[color]++;
                smallK += comp.size;
                smallKBYC[color] += comp.size;
            } else if (comp.size <= b2) {
                if (res_mode==COMPONENT_SIZES_MODE.ALL || res_mode == COMPONENT_SIZES_MODE.GOOD) {
                    ans.add(comp);
                }
                ok++;
                okBYC[color]++;
                okK += comp.size;
                okKBYC[color] += 1;
            } else {
                if (res_mode==COMPONENT_SIZES_MODE.ALL || res_mode == COMPONENT_SIZES_MODE.BIG) {
                    ans.add(comp);
                }
                big++;
                bigBYC[color]++;
                bigKbYC[color] += comp.size;
            }
        }
        System.out.println("total size: " + ans.size());
        Tool.info(logger, "Found " + NumUtils.groupDigits(ok) + " good components, " +
                "and " + NumUtils.groupDigits(big) + " big ones");
        Tool.info(logger, "Found " + NumUtils.groupDigits(small) + " small components, ");
        Tool.info(logger, "By colors found: ");
        for (int i = 0; i < colorsCNT; i++) {
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
        for (int i = 0; i < colorsCNT; i++) {
            Tool.debug(logger, "#" + i + ": " + "Components count: small = " + withP(smallBYC[i], newComps.size()) + ", " +
                    "ok = " + withP(okBYC[i], newComps.size()) + ", " +
                    "big = " + withP(bigBYC[i], newComps.size()));
            Tool.debug(logger, "#" + i + ": " + "Components kmers: small = " + withP(smallKBYC[i], hmSize) + ", " +
                    "ok = " + withP(okKBYC[i], hmSize) + ", " +
                    //неправильно в случае по цветам (возможно, если переполнение)
                    "big = " + withP(bigKbYC[i], hmSize));
        }

        Tool.debug(logger, "FreqThreshold = " + curFreqThreshold + ", " +
                "components added = " + ok + ", total components added = " + ans.size());

        Tool.debug(logger, "Memory used: without GC = " + Misc.usedMemoryWithoutRunningGCAsString() + ", " +
                "after it = " + Misc.usedMemoryAsString());


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
