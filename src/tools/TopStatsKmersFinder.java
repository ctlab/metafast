package tools;

import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.MutableLongLongEntry;
import ru.ifmo.genetics.utils.Misc;
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
import java.util.Iterator;
import java.util.stream.IntStream;
import java.util.stream.Stream;


public class TopStatsKmersFinder extends Tool {

    public static final String NAME = "top-stats-kmers";
    public static final String DESCRIPTION = "Output specified number of best (by p-value) k-mers that statistically" +
            " significant to one of two or three groups of samples based on chi-squared test and save files for faster " +
            "subsequent extractions";

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

        info("Applying chi-squared test...");
        BigLong2LongHashMap hm_chisq = new BigLong2LongHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        Iterator<MutableLongBitShortaEntry> it_all = allKmers.entryIterator();
        Double[] pvalues = new Double[(int) allKmers.size()];
        while (it_all.hasNext()) {
            MutableLongBitShortaEntry entry = it_all.next();
            long key = entry.getKey();

            int n_1_A = allKmers.getCardinality(key, 0, Alength);
            int n_1_B = allKmers.getCardinality(key, Alength, Alength + Blength);
            int n_0_A = Alength - n_1_A;
            int n_0_B = Blength - n_1_B;

            int n_0_C = 0;
            int n_1_C = 0;
            if (is3classes) {
                n_1_C = allKmers.getCardinality(key, Alength + Blength, totalLength);
                n_0_C = Clength - n_1_C;
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

            double chisq;
            if (is3classes) {
                chisq = chisq_3gr(n_0_A, n_1_A, n_0_B, n_1_B, n_0_C, n_1_C);
            } else {
                chisq = chisq_2gr(n_0_A, n_1_A, n_0_B, n_1_B);
            }

            pvalues[n] = chisq;
            hm_chisq.put(key, n);
            n++;
        }

        File allKmersFile = new File(outDir, "all.kmers.bin");
        File allRanksFile = new File(outDir, "all_chi_squared_ranks.bin");
        DataOutputStream streamKmers = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(allKmersFile), 1 << 24));
        DataOutputStream streamRanks = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(allRanksFile), 1 << 24));

        File filteredKmersFile = new File(outDir, "top_" + NBest.get().toString() + "_chi_squared_specific.kmers.bin");
        DataOutputStream stream = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(filteredKmersFile), 1 << 24));

        int[] sortedIndices = IntStream.range(0, n)
                .boxed().sorted((i, j) -> pvalues[j].compareTo(pvalues[i])) // sort in descending order
                .mapToInt(ele -> ele).toArray();

        int[] kmersRanks = new int[n];
        for (int i = 0; i < n; i++) {
            kmersRanks[sortedIndices[i]] = i;
        }

        int i = 0;
        int n_left = 0;
        Iterator<MutableLongLongEntry> it = hm_chisq.entryIterator();
        while (it.hasNext()) {
            MutableLongLongEntry entry = it.next();
            long key = entry.getKey();
            int value = (int) entry.getValue();
            streamKmers.writeLong(key);
            streamKmers.writeShort((short) 1);
            streamRanks.writeInt(kmersRanks[value]);
            if (kmersRanks[i] < NBest.get()) {
                stream.writeLong(key);
                stream.writeShort((short) 1);
                n_left+=1;
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


    public static double chisq_3gr(float c0, float c1, float p0, float p1, float q0, float q1) {
        float tmp = c0;
        c0 = 100 * c0 / (c0 + c1);
        c1 = 100 * c1 / (tmp + c1);
        tmp = p0;
        p0 = 100 * p0 / (p0 + p1);
        p1 = 100 * p1 / (tmp + p1);
        tmp = q0;
        q0 = 100 * q0 / (q0 + q1);
        q1 = 100 * q1 / (tmp + q1);

        float gr_1 = c0 + c1;
        float gr_2 = p0 + p1;
        float gr_3 = q0 + q1;
        float all = gr_1 + gr_2 + gr_3;
        float x1 = gr_1 / all * (p1 + c1 + q1);
        float x2 = gr_1 / all * (p0 + c0 + q0);
        float x3 = gr_2 / all * (p1 + c1 + q1);
        float x4 = gr_2 / all * (p0 + c0 + q0);
        float x5 = gr_3 / all * (p1 + c1 + q1);
        float x6 = gr_3 / all * (p0 + c0 + q0);
        return Math.pow(Math.abs(p1 - x1) - 0.5, 2) / x1 + Math.pow(Math.abs(p0 - x2) - 0.5, 2) / x2 + Math.pow(Math.abs(c1 - x3) - 0.5, 2) / x3 + Math.pow(Math.abs(c0 - x4) - 0.5, 2) / x4 + Math.pow(Math.abs(q1 - x5) - 0.5, 2) / x5 + Math.pow(Math.abs(q0 - x6) - 0.5, 2) / x6;
    }

    public static double chisq_2gr(float c0, float c1, float p0, float p1) {
        float tmp = c0;
        c0 = 100 * c0 / (c0 + c1);
        c1 = 100 * c1 / (tmp + c1);
        tmp = p0;
        p0 = 100 * p0 / (p0 + p1);
        p1 = 100 * p1 / (tmp + p1);

        float gr_1 = c0 + c1;
        float gr_2 = p0 + p1;
        float all = gr_1 + gr_2;
        float x1 = gr_1 / all * (p1 + c1);
        float x2 = gr_1 / all * (p0 + c0);
        float x3 = gr_2 / all * (p1 + c1);
        float x4 = gr_2 / all * (p0 + c0);
        return Math.pow(Math.abs(p1 - x1) - 0.5, 2) / x1 + Math.pow(Math.abs(p0 - x2) - 0.5, 2) / x2 + Math.pow(Math.abs(c1 - x3) - 0.5, 2) / x3 + Math.pow(Math.abs(c0 - x4) - 0.5, 2) / x4;
    }

    public TopStatsKmersFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new TopStatsKmersFinder().mainImpl(args);
    }
}