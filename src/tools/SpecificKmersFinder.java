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

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;

public class SpecificKmersFinder extends Tool {

    public static final String NAME = "specific-kmers";

    public static final String DESCRIPTION = "Output k-mers specific to groups of samples based on chi-squared & Mann-Whitney tests";


    public final Parameter<File[]> Afiles = addParameter(new FileMVParameterBuilder("a-kmers") //ibd
            .mandatory()
            .withShortOpt("A")
            .withDescription("list of input files with k-mers in binary format for group A")
            .create());

    public final Parameter<File[]> Bfiles = addParameter(new FileMVParameterBuilder("b-kmers") //nonibd
            .mandatory()
            .withShortOpt("B")
            .withDescription("list of input files with k-mers in binary format for group B")
            .create());

    public final Parameter<Double> PValueChi2 = addParameter(new DoubleParameterBuilder("p-value-chi2")
            .withShortOpt("pchi2")
            .withDescription("p-value for chi-squared test")
            .create());

    public final Parameter<Double> PValueMW = addParameter(new DoubleParameterBuilder("p-value-mw")
            .withShortOpt("pmw")
            .withDescription("p-value for Mann-Whitney test")
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
        int A_length = Afiles.get().length;
        int B_length = Bfiles.get().length;

        for (File file : Afiles.get()) {
            kmersHMs.add(IOUtils.loadKmers(new File[]{file}, 0, availableProcessors.get(), logger));
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }
        for (File file : Bfiles.get()) {
            kmersHMs.add(IOUtils.loadKmers(new File[]{file}, 0, availableProcessors.get(), logger));
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        }

        ChiSquaredDistribution xi = new ChiSquaredDistributionImpl(1, 0.0001);
        double qvalue;
        try {
            qvalue = xi.inverseCumulativeProbability(1 - PValueChi2.get());
        } catch (MathException e) {
            throw new ExecutionFailedException("Error calculating chi-squared value!", e);
        }

        long total_kmers = 0;
        long unique = 0;
        long unique_left = 0;
        long skip_chi2 = 0;
        long skip_MW = 0;
        long a_kmers = 0;
        long b_kmers = 0;

        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }
        File fileA = new File(outDir, "filtered_groupA.kmers.bin");
        File fileB = new File(outDir, "filtered_groupB.kmers.bin");
        DataOutputStream streamA = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileA), 1 << 24));
        DataOutputStream streamB = new DataOutputStream(new BufferedOutputStream(
                new FileOutputStream(fileB), 1 << 24));


        info("Splitting k-mers...");
        int cur_file = 0;
        for (BigLong2ShortHashMap hm : kmersHMs) {
            debug("Started processing file #" + cur_file);
            Iterator<MutableLongShortEntry> it = hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                if (entry.getValue() > 0) {
                    total_kmers++;
                    int num_files = 0;
                    int n_0_A = 0;
                    int n_0_B = 0;
                    double[] groupA = new double[A_length];
                    double[] groupB = new double[B_length];

                    for (BigLong2ShortHashMap tmp_hm : kmersHMs) {
                        short d = tmp_hm.get(key);
                        if (d > 0) {
                            if (num_files < A_length) {
                                groupA[num_files] = d;
                            } else {
                                groupB[num_files - A_length] = d;
                            }
                            tmp_hm.put(key, (short) 0);
                        } else {
                            if (num_files < A_length) {
                                n_0_A += 1;
                                groupA[num_files] = 0;
                            } else {
                                n_0_B += 1;
                                groupB[num_files - A_length] = 0;
                            }
                        }
                        num_files++;
                    }

                    int n_1_A = A_length - n_0_A;
                    int n_1_B = B_length - n_0_B;


                    boolean chisq_bool = chisq(n_0_A, n_1_A, n_0_B, n_1_B, qvalue);
                    boolean isUniq = (n_1_A == 0) || (n_1_B == 0);
                    if (isUniq) {
                        chisq_bool = true;
                        unique++;
                    }
                    if((n_0_A+n_0_B)==0)    //in all files
                        chisq_bool = true;

                    if (chisq_bool) {
                        int pass = 0;
                        if (PValueMW.get() > 0) {
                            MannWhitneyUTest mw = new MannWhitneyUTest();
                            double pvalue = mw.mannWhitneyUTest(groupA, groupB);
                            if (pvalue > PValueMW.get())
                                pass = 1;
                        }
                        if (pass == 1) {
                            skip_MW++;
                        } else {
                            if (isUniq) {
                                unique_left++;
                            }
                            double mean_A = mean(groupA);
                            double mean_B = mean(groupB);
                            if (mean_A > mean_B) {
                                a_kmers++;
                                streamA.writeLong(key);
                                streamA.writeShort((int) mean_A);
                            } else {
                                b_kmers++;
                                streamB.writeLong(key);
                                streamB.writeShort((int) mean_B);
                            }
                        }
                    } else {
                        skip_chi2++;
                    }
                }
            }
            debug("Finished processing file #" + cur_file);
            debug("Processed " + total_kmers + " k-mers");
            cur_file++;
        }

        streamA.close();
        streamB.close();
        long totalkmersleft = a_kmers+b_kmers;
        info("Group A k-mers printed to " + fileA.getPath());
        info("Group B k-mers printed to " + fileB.getPath());
        info("Total specific k-mers in Group A = " + a_kmers);
        info("Total specific k-mers in Group B = " + b_kmers);

        debug("Total unique kmers = " + unique);
        debug("Total skipped by chi-squared test = " + skip_chi2);
        debug("Total skipped by Mann-Whitney test = " + skip_MW);
        debug("Total unique left = " + unique_left);
        debug("Total kmers left = " + totalkmersleft);
        debug("Parameters used : PValueChi = " + PValueChi2.get() + " that gives " + qvalue + " threshold" );

        debug("Specific-kmers has finished! Time = " + t);
        debug("Processed " + total_kmers + " k-mers");

    }

    private boolean chisq(int c0, int c1, int p0, int p1, double value) {
        int gr_1 = c0 + c1;
        int gr_2 = p0 + p1;
        float all = gr_1 + gr_2;
        float x1 = gr_1 / all * (p1 + c1);
        float x2 = gr_1 / all * (p0 + c0);
        float x3 = gr_2 / all * (p1 + c1);
        float x4 = gr_2 / all * (p0 + c0);
        double kk = Math.pow(Math.abs(p1 - x1) - 0.5, 2) / x1 + Math.pow(Math.abs(p0 - x2)- 0.5, 2) / x2 + Math.pow(Math.abs(c1 - x3)- 0.5, 2) / x3 + Math.pow(Math.abs(c0 - x4)- 0.5, 2) / x4;
        return value < kk;
        }

    private double mean(double[] m) {
        double sum = 0;
        for (double v : m) {
            sum += v;
        }
        return sum / m.length;
    }

    public SpecificKmersFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new SpecificKmersFinder().mainImpl(args);
    }

}
