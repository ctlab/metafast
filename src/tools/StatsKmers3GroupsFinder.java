package tools;

import io.IOUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.pairs.Pair;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.DoubleParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.map.BigLong2BitLongaHashMap;
import structures.map.MutableLongBitLongaEntry;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.stream.Stream;

/**
 Implementation of statistical tests with custom bitset based on long[].
 */
public class StatsKmers3GroupsFinder extends Tool {

    public static final String NAME = "stats-kmers-3";
    public static final String DESCRIPTION = "Output k-mers statistically significant to each of three groups of samples based on chi-squared & Mann-Whitney test";

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
            .mandatory()
            .withShortOpt("C")
            .withDescription("list of input files with k-mers in binary format for group C")
            .create());

    public final Parameter<Double> PValueChi2 = addParameter(new DoubleParameterBuilder("p-value-chi2")
            .withShortOpt("pchi2")
            .withDescription("p-value for chi-squared test")
            .withDefaultValue(0.05)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    public final Parameter<Double> PValueMW = addParameter(new DoubleParameterBuilder("p-value-mw")
            .withShortOpt("pmw")
            .withDescription("p-value for Mann-Whitney test")
            .withDefaultValue(0.05)
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
        int Clength = Cfiles.get().length;
        int totalLength = Alength + Blength + Clength;

        File[] all_files = Stream.of(Afiles.get(), Bfiles.get(), Cfiles.get()).flatMap(Stream::of).toArray(File[]::new);
        BigLong2BitLongaHashMap allKmers = IOUtils.loadBitLongaKmers(all_files, 0, availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

        ChiSquaredDistribution xi = new ChiSquaredDistributionImpl(1, 1e-15D);
        double qvalue;
        try {
            qvalue = xi.inverseCumulativeProbability(1 - PValueChi2.get());
        } catch (MathException e) {
            throw new ExecutionFailedException("Error calculating chi-squared value!", e);
        }


        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }
        File filteredKmers = new File(outDir, "filtered_chisquared.kmers.bin");
        File stFile = new File(outDir, "filtered_chisquared.stat.txt");

        long n = 0;
        long nUnique = 0;
        long nUniqueLeft = 0;
        long nFilteredChi = 0;
        long nFilteredMW = 0;
        long nA = 0;
        long nB = 0;
        long nC = 0;
        long nInAll = 0;
        long nScarce = 0;

        info("Applying chi-squared test...");
        BigLong2ShortHashMap hm_chisq = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        Iterator<MutableLongBitLongaEntry> it_all = allKmers.entryIterator();
        while (it_all.hasNext()) {
            n++;
            MutableLongBitLongaEntry entry = it_all.next();
            long key = entry.getKey();

            int n_1_A = allKmers.getCardinality(key, 0, Alength);
            int n_1_B = allKmers.getCardinality(key, Alength, Alength+Blength);
            int n_1_C = allKmers.getCardinality(key, Alength+Blength, totalLength);

            int n_0_A = Alength - n_1_A;
            int n_0_B = Blength - n_1_B;
            int n_0_C = Clength - n_1_C;

            if (n_1_A+n_1_B+n_1_C <= Math.ceil(totalLength * 0.05)) { //skip scarce k-mers
                nScarce++;
                continue;
            }
            if (n_1_A+n_1_B+n_1_C == totalLength) {  // if present in all files then useless
                nInAll++;
                continue;
            }

            boolean isUniq = (n_1_A + n_1_C) == 0 || (n_1_B + n_1_A) == 0 || (n_1_B + n_1_C) == 0; // if present only in one group then its unique
            if (isUniq) {
                nUnique++;
            }

            boolean chisq_bool = chisq(n_0_A, n_1_A, n_0_B, n_1_B, n_0_C, n_1_C, qvalue);

            if (chisq_bool) {
                hm_chisq.put(key, (short)1);
            } else {
                nFilteredChi++;
            }
        }

        long c = IOUtils.printKmers(hm_chisq, 0, filteredKmers, stFile);
        if (c != n - nInAll - nScarce - nFilteredChi) {
            throw new ExecutionFailedException("Wrong count of k-mers left after chi-squared: " + c + "(expected:" + (n - nInAll - nScarce - nFilteredChi) + ")");
        }
        filteredKmersFilePr.set(filteredKmers);
        info("Survived after chi-squared test k-mers printed to: " + filteredKmers.getPath());

        allKmers.reset();
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


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
            Iterator<MutableLongShortEntry> it = hm_chisq.entryIterator();
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
            Iterator<MutableLongShortEntry> it = hm_chisq.entryIterator();
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

        i = 0;
        for (File file : Cfiles.get()) {
            BigLong2ShortHashMap add_hm = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
            Pair<BigLong2ShortHashMap, Long> tmp = IOUtils.loadKmersFreq(new File[]{file}, 0, availableProcessors.get(), logger);
            BigLong2ShortHashMap tmp_hm = tmp.first();
            Iterator<MutableLongShortEntry> it = hm_chisq.entryIterator();
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

        info("Applying Mann-Whitney test...");
        MannWhitneyUTest mw = new MannWhitneyUTest();
        File fileA = new File(outDir, "filtered_groupA.kmers.bin");
        File fileB = new File(outDir, "filtered_groupB.kmers.bin");
        File fileC = new File(outDir, "filtered_groupC.kmers.bin");

        DataOutputStream streamA = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileA), 1 << 24));
        DataOutputStream streamB = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileB), 1 << 24));
        DataOutputStream streamC = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileC), 1 << 24));

        double meanSumKmers = Stream.of(AkmersFreq, BkmersFreq, CkmersFreq).flatMapToDouble(Arrays::stream).sum() / totalLength;
        Iterator<MutableLongShortEntry> it = hm_chisq.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            double[] groupA = new double[Alength];
            double[] groupB = new double[Blength];
            double[] groupC = new double[Clength];

            final long key = entry.getKey();
            double[] tmp_freqA = AkmersHMs.stream().mapToDouble(tmp->tmp.getWithZero(key)).map(v->v*meanSumKmers).toArray();
            Arrays.setAll(groupA, j -> tmp_freqA[j] / AkmersFreq[j]);
            double[] tmp_freqB = BkmersHMs.stream().mapToDouble(tmp->tmp.getWithZero(key)).map(v->v*meanSumKmers).toArray();
            Arrays.setAll(groupB, j -> tmp_freqB[j] / BkmersFreq[j]);
            double[] tmp_freqC = CkmersHMs.stream().mapToDouble(tmp->tmp.getWithZero(key)).map(v->v*meanSumKmers).toArray();
            Arrays.setAll(groupC, j -> tmp_freqC[j] / CkmersFreq[j]);


            int pass = 0;
            if (PValueMW.get() > 0) {   // if less than zero, then skip this phase of filtration
                double A_B = mw.mannWhitneyUTest(groupA, groupB);
                double B_C = mw.mannWhitneyUTest(groupB, groupC);
                double A_C = mw.mannWhitneyUTest(groupA, groupC);
                if (A_B < PValueMW.get() || B_C < PValueMW.get() || A_C < PValueMW.get()) {
                    pass = 1;
                }
            } else {
                pass = 1;
            }

            if (pass == 1) {
                double meanA = mean(groupA);
                double meanB = mean(groupB);
                double meanC = mean(groupC);

                if (meanA > meanB && meanA > meanC) {
                    streamA.writeLong(entry.getKey());
                    streamA.writeShort((int) meanA);
                    nA++;
                } else {
                    if (meanB > meanA && meanB > meanC) {
                        streamB.writeLong(entry.getKey());
                        streamB.writeShort((int) meanB);
                        nB++;
                    } else {
                        streamC.writeLong(entry.getKey());
                        streamC.writeShort((int) meanC);
                        nC++;
                    }
                }

                if ((meanA + meanB) == 0 || (meanA + meanC) == 0 || (meanB + meanC) == 0) {
                    nUniqueLeft++;
                }
            } else {
                nFilteredMW++;
            }
        }

        streamA.close();
        streamB.close();
        streamC.close();

        info("Group A k-mers printed to " + fileA.getPath());
        info("Group B k-mers printed to " + fileB.getPath());
        info("Group C k-mers printed to " + fileC.getPath());


        long n_left = nA + nB + nC;
        if (n_left != n - nInAll - nScarce - nFilteredChi - nFilteredMW) {
            throw new ExecutionFailedException("Wrong count of k-mers left after chi-squared: " + n_left + "(expected:" + (n - nInAll - nScarce - nFilteredChi - nFilteredMW) + ")");
        }


        debug("Total k-mers count = " + n);
        debug("Total unique k-mers = " + nUnique);
        debug("Total k-mers present in all files = " + nInAll);
        debug("Total k-mers left = " + n_left);
        debug("Total unique left = " + nUniqueLeft);
        info("Total group A k-mers = " + nA);
        info("Total group B k-mers = " + nB);
        info("Total group C k-mers = " + nC);
        debug("Total scarce k-mers = " + nScarce);
        debug("Total skipped by Chi-squared test = " + nFilteredChi);
        debug("Total skipped by Mann-Whitney test = " + nFilteredMW);
        debug("Parameters used : P-value for Chi-squared test = " + PValueChi2.get() + " that gives " + qvalue + " threshold, p-value for Kruskal-Wallis test = " + PValueMW.get());
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        info("chi-squared-kmers-3 has finished! Time elapsed = " + t);
    }


    public static boolean chisq(float c0, float c1, float p0, float p1, float q0, float q1, double value) {
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
        double stat = Math.pow(Math.abs(p1 - x1) - 0.5, 2) / x1 + Math.pow(Math.abs(p0 - x2) - 0.5, 2) / x2 + Math.pow(Math.abs(c1 - x3) - 0.5, 2) / x3 + Math.pow(Math.abs(c0 - x4) - 0.5, 2) / x4 + Math.pow(Math.abs(q1 - x5) - 0.5, 2) / x5 + Math.pow(Math.abs(q0 - x6) - 0.5, 2) / x6;
        return value < stat;
    }

    private static double mean(double[] m) {
        double sum = 0;
        for (double v : m) {
            sum += v;
        }
        return sum / m.length;
    }

    public StatsKmers3GroupsFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new StatsKmers3GroupsFinder().mainImpl(args);
    }
}