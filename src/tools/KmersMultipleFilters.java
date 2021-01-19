package tools;


import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;
import java.io.IOException;

/**
 * Created by -- on 19.01.2021.
 */
public class KmersMultipleFilters extends Tool {

    public static final String NAME = "kmers-multiple-filters";

    public static final String DESCRIPTION = "Filter k-mers from test set according to three specified sets";


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size (maximum 31 due to realization details)")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("i")
            .withDescription("list of input files with k-mers in binary format")
            .create());

    public final Parameter<File[]> CDfilterFiles = addParameter(new FileMVParameterBuilder("cd-filter-kmers")
            .mandatory()
            .withShortOpt("cd")
            .withDescription("list of input files with k-mers in binary format used for CD filtering")
            .create());

    public final Parameter<File[]> UCfilterFiles = addParameter(new FileMVParameterBuilder("uc-filter-kmers")
            .mandatory()
            .withShortOpt("uc")
            .withDescription("list of input files with k-mers in binary format used for UC filtering")
            .create());

    public final Parameter<File[]> nonIBDfilterFiles = addParameter(new FileMVParameterBuilder("nonibd-filter-kmers")
            .mandatory()
            .withShortOpt("nonibd")
            .withDescription("list of input files with k-mers in binary format used for nonIBD filtering")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous (while saving kmers from sequences)")
            .withDefaultValue(1)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    public final Parameter<File> statsDir = addParameter(new FileParameterBuilder("stats-dir")
            .withDescription("Directory with statistics")
            .withDefaultValue(workDir.append("stats"))
            .create());


    @Override
    protected void runImpl() throws ExecutionFailedException {
        if (k.get() <= 0) {
            error("The size of k-mer must be at least 1.");
            System.exit(1);
        }
        if (k.get() > 31) {
            error("The size of k-mer must be no more than 31.");
            System.exit(1);
        }


        Timer t = new Timer();
        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }
        File stDir = statsDir.get();
        if (!stDir.exists()) {
            stDir.mkdirs();
        }

        BigLong2ShortHashMap cd_filter_hm = IOUtils.loadKmers(CDfilterFiles.get(), 0,
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        BigLong2ShortHashMap uc_filter_hm = IOUtils.loadKmers(UCfilterFiles.get(), 0,
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        BigLong2ShortHashMap nonibd_filter_hm = IOUtils.loadKmers(nonIBDfilterFiles.get(), 0,
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        for (File file : inputFiles.get()) {
            BigLong2ShortHashMap hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            String name = file.getName().replaceAll(".kmers.bin", "");
            File outFile = new File(outDir, name + ".kmers.bin");
            File stFile = new File(stDir, name + ".stat.txt");


            debug("Starting to print k-mers to " + outFile.getPath());
            long c = 0;
            try {
                c = IOUtils.MultipleFiltersAndPrintKmers(hm, cd_filter_hm, uc_filter_hm, nonibd_filter_hm,
                        maximalBadFrequency.get(), outFile, stFile);
            } catch (IOException e) {
                e.printStackTrace();
            }
            info(NumUtils.groupDigits(hm.size()) + " k-mers found, " + NumUtils.groupDigits(c) +
                    " (" + String.format("%.1f", c * 100.0 / hm.size()) + "%) of them survived after filtering");

            info("Filtered k-mers printed to " + outFile.getPath());
        }
    }

    @Override
    protected void cleanImpl() {
    }

    public KmersMultipleFilters() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new KmersMultipleFilters().mainImpl(args);
    }


}
