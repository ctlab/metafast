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
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.DoubleParameterBuilder;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;


public class SpecificKmers3GroupsFinder extends Tool {

    public static final String NAME = "specific-kmers-3";
    public static final String DESCRIPTION = "Output k-mers specific to each of three groups of samples based on chi-squared & Mann-Whitney tests";

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

    public final Parameter<Double> PValueMW = addParameter(new DoubleParameterBuilder("p-value-mw")
            .withShortOpt("pmw")
            .withDescription("p-value for Mann-Whitney test")
            .withDefaultValue(0.05)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        ArrayList<BigLong2ShortHashMap> kmersHMs = new ArrayList<>();
        int aLength = Afiles.get().length;
        int bLength = Bfiles.get().length;
        int cLength = Cfiles.get().length;
        int totalLength = aLength + bLength + cLength;

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
                hm_all.put(key, (short) (hm_all.getWithZero(key) + 1));
                n_kmers = n_kmers + value;
            }
            files_kmers_hm.put(tmp_hm.toString(), n_kmers);
        }

        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }

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
        int curFile = 0;
        Iterator<MutableLongShortEntry> it_all = hm_all.entryIterator();
        while (it_all.hasNext()) {
            n++;
            MutableLongShortEntry entry = it_all.next();
            curFile = 0;
            if (entry.getValue() <= Math.ceil(totalLength * 0.05)) { //skip scarce k-mers
                nScarce++;
                continue;
            }
            if (entry.getValue() == totalLength) {  // if present in all files then useless
                nInAll++;
                continue;
            }

            int n_0_A = 0;
            int n_0_B = 0;
            int n_0_C = 0;
            double[] groupA = new double[aLength];
            double[] groupB = new double[bLength];
            double[] groupC = new double[cLength];

            for (BigLong2ShortHashMap tmp_hm : kmersHMs) {
                double d = tmp_hm.get(entry.getKey());
                if (d != -1) {
                    d = (d * files_kmers_hm.get(kmersHMs.get(0).toString()) /  files_kmers_hm.get(tmp_hm.toString()));
                    if (curFile < aLength) {
                        groupA[curFile] = d;
                    } else {
                        if (curFile < aLength + bLength) {
                            groupB[curFile - aLength] = d;
                        } else {
                            groupC[curFile - aLength - bLength] = d;
                        }
                    }
                } else {
                    if (curFile < aLength) {
                        n_0_A += 1;
                        groupA[curFile] = 0;
                    } else {
                        if (curFile < aLength + bLength) {
                            n_0_B += 1;
                            groupB[curFile - aLength] = 0;
                        } else {
                            n_0_C += 1;
                            groupC[curFile - aLength - bLength] = 0;
                        }
                    }
                }
                curFile = curFile + 1;
            }

            int n_1_A = aLength - n_0_A;
            int n_1_B = bLength - n_0_B;
            int n_1_C = cLength - n_0_C;

            boolean isUniq = (n_1_A + n_1_C) == 0 || (n_1_B + n_1_A) == 0 || (n_1_B + n_1_C) == 0; // if present only in one group then its unique
            if (isUniq)
                nUnique++;

            boolean chisq_bool = chisq(n_0_A, n_1_A, n_0_B, n_1_B, n_0_C, n_1_C, qvalue);

            if (chisq_bool) {
                int pass = 0;
                if (PValueMW.get() > 0) {   // if less then zero then skip this phase of filtration
                    MannWhitneyUTest mw = new MannWhitneyUTest();
                    double A_B = mw.mannWhitneyUTest(groupA, groupB);
                    double B_C = mw.mannWhitneyUTest(groupB, groupC);
                    double A_C = mw.mannWhitneyUTest(groupA, groupC);
                    if (!(A_B < PValueMW.get() || B_C < PValueMW.get() || A_C < PValueMW.get()))
                        pass = 1;
                }

                if (pass == 1){
                    nFilteredMW++;
                }
                else {
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
                    if (isUniq)
                        nUniqueLeft++;
                }
            } else nFilteredChi++;
        }

        streamA.close();
        streamB.close();
        streamC.close();
        long n_left = nA + nB + nC;

        info("Group A k-mers printed to " + fileA.getPath());
        info("Group B k-mers printed to " + fileB.getPath());
        info("Group C k-mers printed to " + fileC.getPath());

        debug("Total k-mers count = " + n);
        debug("Total unique k-mers = " + nUnique);
        debug("Total k-mers present in all files = " + nInAll);
        debug("Total k-mers left = " + n_left);
        debug("Total unique left = " + nUniqueLeft);
        debug("Total group A k-mers = " + nA);
        debug("Total group B k-mers = " + nB);
        debug("Total group C k-mers = " + nC);
        info("Total scarce k-mers = " + nScarce);
        debug("Total skipped by Chi-squared test = " + nFilteredChi);
        debug("Total skipped by Mann-Whitney test = " + nFilteredMW);
        debug("Parameters used : P-value for Chi-squared test = " + PValueChi2.get() + " that gives " + qvalue + " threshold, p-value for Kruskal-Wallis test = " + PValueMW.get());
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        info("Specific-k-mers has finished! Time = " + t);
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