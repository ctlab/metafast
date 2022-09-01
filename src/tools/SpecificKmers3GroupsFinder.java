package tools;

import io.IOUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.moeaframework.util.statistics.KruskalWallisTest;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;


public class SpecificKmers3GroupsFinder extends Tool {

    public static final String NAME = "specific-kmers-3";
    public static final String DESCRIPTION = "Output k-mers specific to each of three groups of samples based on chi-squared & Kruskal-Wallis tests";

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

    public final Parameter<Double> PValueKW = addParameter(new DoubleParameterBuilder("p-value-kw")
            .withShortOpt("pkw")
            .withDescription("p-value for Kruskal-Wallis test")
            .withDefaultValue(0.05)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        ArrayList<BigLong2ShortHashMap> kmersHMs = new ArrayList<>();
        int A_length = Afiles.get().length;
        int B_length = Bfiles.get().length;
        int C_length = Cfiles.get().length;

        for (File file : Afiles.get()) {
            kmersHMs.add(IOUtils.loadKmers(new File[]{file}, 0, availableProcessors.get(), logger));
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }
        for (File file : Bfiles.get()) {
            kmersHMs.add(IOUtils.loadKmers(new File[]{file}, 0, availableProcessors.get(), logger));
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }
        for (File file : Cfiles.get()) {
            kmersHMs.add(IOUtils.loadKmers(new File[]{file}, 0, availableProcessors.get(), logger));
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }

        ChiSquaredDistribution xi = new ChiSquaredDistributionImpl(1, 0.00001);
        double qvalue;
        try {
            qvalue = xi.inverseCumulativeProbability(1 - PValueChi2.get());
        } catch (MathException e) {
            throw new ExecutionFailedException("Error calculating chi-squared value!", e);
        }

        BigLong2ShortHashMap hm_all =  new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 8);
        hm_all.resetValues();
        HashMap<String,Long> files_kmers_hm = new HashMap<String,Long>();
        for (BigLong2ShortHashMap tmp_hm : kmersHMs) {
            long n_kmers = 0;
            Iterator<MutableLongShortEntry> it = tmp_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                short value = entry.getValue();
                if (value > maximalBadFrequency.get()) {
                    hm_all.put(key, (short) (hm_all.getWithZero(key) + 1));
                }
                n_kmers = n_kmers + value;
            }
            files_kmers_hm.put(tmp_hm.toString(), n_kmers);
        }

        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }
        File stDir = workDir.append("kmers").get();
        if (!stDir.exists()) {
            stDir.mkdirs();
        }
        File outFile = new File(outDir, "n_samples.kmers.bin");
        File stFile = new File(stDir, "n_samples.stat.txt");

        debug("Starting to print k-mers to " + outFile.getPath());
        long c = 0;
        try {
            c = IOUtils.printKmers(hm_all, 0, outFile, stFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        info(NumUtils.groupDigits(hm_all.size()) + " k-mers found, "
                + NumUtils.groupDigits(c) + " (" + String.format("%.1f", c * 100.0 / hm_all.size()) + "%) of them is good (not erroneous)");

        long n = 0;
        long n_unique = 0;
        long n_unique_left = 0;
        long n_filtered_chi = 0;
        long n_filtered_2 = 0;
        long n_A = 0;
        long n_B = 0;
        long n_C = 0;
        long n_error = 0;
        long n_in_all = 0;

        // File outDir = outputDir.get();
        if (!outDir.exists())
            outDir.mkdirs();
        File fileA = new File(outDir, "filtered_groupA.kmers.bin");
        File fileB = new File(outDir, "filtered_groupB.kmers.bin");
        File fileC = new File(outDir, "filtered_groupC.kmers.bin");

        DataOutputStream streamA = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileA), 1 << 24));
        DataOutputStream streamB = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileB), 1 << 24));
        DataOutputStream streamC = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileC), 1 << 24));

        info("Splitting k-mers...");
        int cur_file = 0;
        Iterator<MutableLongShortEntry> it_all = hm_all.entryIterator();
        while (it_all.hasNext()) {
            n++;
            MutableLongShortEntry entry = it_all.next();
            cur_file = 0;

            int n_0_A = 0;
            int n_0_B = 0;
            int n_0_C = 0;
            double[] groupA = new double[A_length];
            double[] groupB = new double[B_length];
            double[] groupC = new double[C_length];

            for (BigLong2ShortHashMap tmp_hm : kmersHMs) {
                double d = tmp_hm.get(entry.getKey());
                if (d != -1) {
                    d = (short) (d * files_kmers_hm.get(kmersHMs.get(0).toString()) /  files_kmers_hm.get(tmp_hm.toString()));
                    if (cur_file < A_length) {
                        groupA[cur_file] = d;
                    } else {
                        if (cur_file < A_length + B_length) {
                            groupB[cur_file - A_length] = d;
                        } else {
                            groupC[cur_file - A_length - B_length] = d;
                        }
                    }
                } else {
                    if (cur_file < A_length) {
                        n_0_A += 1;
                        groupA[cur_file] = 0;
                    } else {
                        if (cur_file < A_length + B_length) {
                            n_0_B += 1;
                            groupB[cur_file - A_length] = 0;
                        } else {
                            n_0_C += 1;
                            groupC[cur_file - A_length - B_length] = 0;
                        }
                    }
                }
                cur_file = cur_file + 1;
            }

            int n_1_A = A_length - n_0_A;
            int n_1_B = B_length - n_0_B;
            int n_1_C = C_length - n_0_C;

            if ((n_1_A + n_1_C + n_1_B) == 0) { // if absent in all groups then it`s definitely erroneous
                n_error++;
                continue;
            }

            if ((n_0_A + n_0_B + n_0_C) == 0) {  // if present in all files then useless
                n_in_all++;
                continue;
            }

            boolean isUniq = (n_1_A + n_1_C) == 0 || (n_1_B + n_1_A) == 0 || (n_1_B + n_1_C) == 0; // if present only in one group then its unique
            if (isUniq)
                n_unique++;

            boolean chisq_bool = chisq(n_0_A, n_1_A, n_0_B, n_1_B, n_0_C, n_1_C, qvalue);

            if (chisq_bool || isUniq) { //if k-mer is unique it`s also passes
                int pass = 0;
                if (PValueKW.get() > 0) {   // if less then zero then skip this phase of filtration
                    KruskalWallisTest kwtest = new KruskalWallisTest(3);
                    kwtest.addAll(groupA, 0);
                    kwtest.addAll(groupB, 1);
                    kwtest.addAll(groupC, 2);
                    if (!kwtest.test(PValueKW.get()))
                        pass = 1;
                }

                if (pass == 1)
                    n_filtered_2++;
                else {
                    double mean_A = mean(groupA);
                    double mean_B = mean(groupB);
                    double mean_C = mean(groupC);

                    if (mean_A > mean_B && mean_A > mean_C) {
                        streamA.writeLong(entry.getKey());
                        streamA.writeShort((int) mean_A);
                        n_A++;
                    } else {
                        if (mean_B > mean_A && mean_B > mean_C) {
                            streamB.writeLong(entry.getKey());
                            streamB.writeShort((int) mean_B);
                            n_B++;
                        } else {
                            streamC.writeLong(entry.getKey());
                            streamC.writeShort((int) mean_C);
                            n_C++;
                        }
                    }
                    if (isUniq)
                        n_unique_left++;
                }
            } else n_filtered_chi++;
        }

        streamA.close();
        streamB.close();
        streamC.close();
        long n_left = n_A + n_B + n_C;

        info("Group A k-mers printed to " + fileA.getPath());
        info("Group B k-mers printed to " + fileB.getPath());
        info("Group C k-mers printed to " + fileC.getPath());

        debug("Total k-mers count = " + n);
        debug("Total unique k-mers = " + n_unique);
        debug("Total erroneous k-mers = " + n_error);
        debug("Total k-mers present in all files = " + n_in_all);
        debug("Total k-mers left = " + n_left);
        debug("Total group A k-mers = " + n_A);
        debug("Total group B k-mers = " + n_B);
        debug("Total group C k-mers = " + n_C);
        debug("Total skipped by Chi-squared test = " + n_filtered_chi);
        debug("Total skipped by Kruskal-Wallis test = " + n_filtered_2);
        debug("Total unique left = " + n_unique_left);
        debug("Parameters used : P-value for Chi-squared test = " + PValueChi2.get() + " that gives " + qvalue + " threshold, p-value for Kruskal-Wallis test = " + PValueKW.get());
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        info("Specific-kmers has finished! Time = " + t);
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

    private double mean(double[] m) {
        double sum = 0;
        for (double v : m) {
            sum += v;
        }
        return sum / m.length;
    }

    public SpecificKmers3GroupsFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new SpecificKmers3GroupsFinder().mainImpl(args);
    }
}
