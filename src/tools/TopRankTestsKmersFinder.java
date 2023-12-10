package tools;

import algo.KruskalWallisTest;
import io.IOUtils;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.pairs.Pair;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.map.BigLong2BitShortaHashMap;
import structures.map.MutableLongBitShortaEntry;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.stream.IntStream;
import java.util.stream.Stream;


public class TopRankTestsKmersFinder extends Tool {

    public static final String NAME = "top-rank-test-kmers";
    public static final String DESCRIPTION = "Output specified number of best (by p-value) k-mers that statistically" +
            " significant to one of two or three groups of samples based on Mann-Whitney test (for two groups) and " +
            "and Kruskal-Wallis (for three groups), save files for faster subsequent extractions";

    public final Parameter<File[]> Afiles = addParameter(new FileMVParameterBuilder("a-kmers")
            .mandatory()
            .withShortOpt("A")
            .withDescription("list of input files with k-mers in binary format for group A")
            .create());

    public final Parameter<File[]> Bfiles = addParameter(new FileMVParameterBuilder("b-kmers")
            .mandatory()
            .withShortOpt("B")
            .withDescription("list of input files with k-mers in binary format for group B")
            .create());

    public final Parameter<File[]> Cfiles = addParameter(new FileMVParameterBuilder("c-kmers")
            .optional()
            .withShortOpt("C")
            .withDescription("optional list of input files with k-mers in binary format for group C")
            .create());

    public final Parameter<Integer> NBest = addParameter(new IntParameterBuilder("num-kmers")
            .mandatory()
            .withShortOpt("nk")
            .withDescription("number of most specific k-mers to be extracted")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous (while detecting presence in samples)")
            .withDefaultValue(0)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    private final InMemoryValue<File> filteredKmersFilePr = new InMemoryValue<File>();
    public final InValue<File> filteredKmersFile =
            addOutput("filtered-kmers-file", filteredKmersFilePr, File.class);


    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        info("Loading k-mers occurrences...");
        Timer t = new Timer();
        int Alength = Afiles.get().length;
        int Blength = Bfiles.get().length;
        int Clength = 0;
        File[] all_files;
        boolean is3classes = Cfiles.get() != null;
        if (is3classes) {
            Clength = Cfiles.get().length;
            all_files = Stream.of(Afiles.get(), Bfiles.get(), Cfiles.get()).flatMap(Stream::of).toArray(File[]::new);
        } else {
            all_files = Stream.of(Afiles.get(), Bfiles.get()).flatMap(Stream::of).toArray(File[]::new);
        }

        int totalLength = Alength + Blength + Clength;
        BigLong2BitShortaHashMap allKmers = IOUtils.loadBitShortaKmers(all_files, maximalBadFrequency.get(), availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }

        int n = 0;
        long nUnique = 0;
        long nInAll = 0;
        long nScarce = 0;

        org.apache.commons.math3.distribution.ChiSquaredDistribution dist = new org.apache.commons.math3.distribution.ChiSquaredDistribution(2);
        BigLong2ShortHashMap hm_all = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        Iterator<MutableLongBitShortaEntry> it_all = allKmers.entryIterator();
        while (it_all.hasNext()) {
            MutableLongBitShortaEntry entry = it_all.next();
            long key = entry.getKey();
            int n_1_A = allKmers.getCardinality(key, 0, Alength);
            int n_1_B = allKmers.getCardinality(key, Alength, Alength + Blength);
            int n_1_C = 0;
            if (is3classes) {
                n_1_C = allKmers.getCardinality(key, Alength + Blength, totalLength);
            }

            if (n_1_A + n_1_B + n_1_C <= Math.ceil(totalLength * 0.05)) { //skip scarce k-mers
                nScarce++;
                continue;
            }
            if (n_1_A + n_1_B + n_1_C == totalLength) {  // if present in all files then useless
                nInAll++;
                continue;
            }

            boolean isUniq = (n_1_A + n_1_C) == 0 || (n_1_B + n_1_A) == 0 || (n_1_B + n_1_C) == 0; // if present only in one group then its unique
            if (isUniq) {
                nUnique++;
            }

            hm_all.put(key, (short) 1);
            n++;
        }

        info("Loading k-mers frequencies...");
        ArrayList<BigLong2ShortHashMap> AkmersHMs = new ArrayList<>(Alength);
        ArrayList<BigLong2ShortHashMap> BkmersHMs = new ArrayList<>(Blength);
        ArrayList<BigLong2ShortHashMap> CkmersHMs = new ArrayList<>(Clength);

        double[] AkmersFreq = new double[Alength];
        double[] BkmersFreq = new double[Blength];
        double[] CkmersFreq = new double[Clength];

        int i = 0;
        for (File file : Afiles.get()) {
            BigLong2ShortHashMap add_hm = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
            Pair<BigLong2ShortHashMap, Long> tmp = IOUtils.loadKmersFreq(new File[]{file}, 0, availableProcessors.get(), logger);
            BigLong2ShortHashMap tmp_hm = tmp.first();
            Iterator<MutableLongShortEntry> it = hm_all.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                add_hm.put(key, tmp_hm.getWithZero(key));
            }
            AkmersHMs.add(add_hm);
            AkmersFreq[i] = tmp.second();
            i++;
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }

        i = 0;
        for (File file : Bfiles.get()) {
            BigLong2ShortHashMap add_hm = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
            Pair<BigLong2ShortHashMap, Long> tmp = IOUtils.loadKmersFreq(new File[]{file}, 0, availableProcessors.get(), logger);
            BigLong2ShortHashMap tmp_hm = tmp.first();
            Iterator<MutableLongShortEntry> it = hm_all.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                add_hm.put(key, tmp_hm.getWithZero(key));
            }
            BkmersHMs.add(add_hm);
            BkmersFreq[i] = tmp.second();
            i++;
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }
        if (is3classes) {
            i = 0;
            for (File file : Cfiles.get()) {
                BigLong2ShortHashMap add_hm = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
                Pair<BigLong2ShortHashMap, Long> tmp = IOUtils.loadKmersFreq(new File[]{file}, 0, availableProcessors.get(), logger);
                BigLong2ShortHashMap tmp_hm = tmp.first();
                Iterator<MutableLongShortEntry> it = hm_all.entryIterator();
                while (it.hasNext()) {
                    MutableLongShortEntry entry = it.next();
                    long key = entry.getKey();
                    add_hm.put(key, tmp_hm.getWithZero(key));
                }
                CkmersHMs.add(add_hm);
                CkmersFreq[i] = tmp.second();
                i++;
                debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
            }
        }

        if (is3classes) {
            info("Applying Kruskal-Wallis test...");
        } else {
            info("Applying Mann-Whitney test...");
        }
        MannWhitneyUTest mw = new MannWhitneyUTest();
        double meanSumKmers = Stream.of(AkmersFreq, BkmersFreq, CkmersFreq).flatMapToDouble(Arrays::stream).sum() / totalLength;
        Iterator<MutableLongShortEntry> it = hm_all.entryIterator();
        Double[] pvalues = new Double[(int) hm_all.size()];
        i = 0;
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            double[] groupA = new double[Alength];
            double[] groupB = new double[Blength];
            double[] groupC = new double[Clength];

            final long key = entry.getKey();
            double[] tmp_freqA = AkmersHMs.stream().mapToDouble(tmp -> tmp.getWithZero(key)).map(v -> v * meanSumKmers).toArray();
            Arrays.setAll(groupA, j -> tmp_freqA[j] / AkmersFreq[j]);
            double[] tmp_freqB = BkmersHMs.stream().mapToDouble(tmp -> tmp.getWithZero(key)).map(v -> v * meanSumKmers).toArray();
            Arrays.setAll(groupB, j -> tmp_freqB[j] / BkmersFreq[j]);
            if (is3classes) {
                double[] tmp_freqC = CkmersHMs.stream().mapToDouble(tmp -> tmp.getWithZero(key)).map(v -> v * meanSumKmers).toArray();
                Arrays.setAll(groupC, j -> tmp_freqC[j] / CkmersFreq[j]);
            }

            if (is3classes) {
                algo.KruskalWallisTest KWTest = new KruskalWallisTest(3);
                KWTest.addAll(groupA, 0);
                KWTest.addAll(groupB, 1);
                KWTest.addAll(groupC, 2);
                pvalues[i] = 1.0D - dist.cumulativeProbability(KWTest.test());
            } else {
                pvalues[i] = mw.mannWhitneyUTest(groupA, groupB);
            }

            entry.setValue((short) Math.max(Math.max(mean(groupA), mean(groupB)), mean(groupC)));
            i++;
        }

        String rankTestName = "";
        if (is3classes) {
            rankTestName = "kw";
        } else {
            rankTestName = "mw";
        }

        File allKmersFile = new File(outDir, "all.kmers.bin");
        File allRanksFile = new File(outDir, "all_" + rankTestName + "_ranks.bin");

        DataOutputStream streamKmers = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(allKmersFile), 1 << 24));
        DataOutputStream streamRanks = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(allRanksFile), 1 << 24));

        File filteredKmersFile = new File(outDir, "top_" + NBest.get().toString() + "_" + rankTestName + "_specific.kmers.bin");
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(filteredKmersFile), 1 << 24));

        int[] sortedIndices = IntStream.range(0, n)
                .boxed().sorted((a, b) -> pvalues[a].compareTo(pvalues[b])) // sort in descending order
                .mapToInt(ele -> ele).toArray();

        int[] kmersRanks = new int[n];
        for (i = 0; i < n; i++) {
            kmersRanks[sortedIndices[i]] = i;
        }

        i = 0;
        int n_left = 0;
        it = hm_all.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            int value = (int) entry.getValue();
            streamKmers.writeLong(key);
            streamKmers.writeShort((short) value);
            streamRanks.writeInt(kmersRanks[i]);
            if (kmersRanks[i] < NBest.get()) {
                stream.writeLong(key);
                stream.writeShort((short) value);
                n_left += 1;
            }
            i++;
        }

        streamRanks.close();
        streamKmers.close();
        stream.close();

        debug("All non-unique and not scarce k-mers printed to " + allKmersFile.getPath());
        debug("Ranks for non-unique and not scarce k-mers printed to " + allRanksFile.getPath());
        debug("Filtered k-mers printed to " + filteredKmersFile.getPath());
        debug("Total k-mers count = " + n);
        debug("Total unique k-mers = " + nUnique);
        debug("Total k-mers present in all files = " + nInAll);
        debug("Total k-mers left = " + n_left);
        debug("Total scarce k-mers = " + nScarce);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        info("top-stats-kmers has finished! Time elapsed = " + t);
    }

    private static double mean(double[] m) {
        double sum = 0;
        for (double v : m) {
            sum += v;
        }
        return sum / m.length;
    }

    public TopRankTestsKmersFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new TopRankTestsKmersFinder().mainImpl(args);
    }
}