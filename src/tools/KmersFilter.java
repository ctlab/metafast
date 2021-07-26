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
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;

import java.io.File;
import java.io.IOException;

/**
 * Created by -- on 14.02.2020.
 */
public class KmersFilter extends Tool {

    public static final String NAME = "kmers-filter";

    public static final String DESCRIPTION = "Filter k-mers from test set according to known samples";


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

    public final Parameter<File[]> filterFiles = addParameter(new FileMVParameterBuilder("filter-kmers")
            .mandatory()
            .withDescription("list of input files with k-mers in binary format used for filtering")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous (while saving kmers from sequences)")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> maximalThreshold = addParameter(new IntParameterBuilder("max-thresh")
            .optional()
            .withDescription("maximal frequency for a k-mer in filtering file to be assumed not found")
            .withDefaultValue(0)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    public final Parameter<File> statsDir = addParameter(new FileParameterBuilder("stats-dir")
            .withDescription("Directory with statistics")
            .withDefaultValue(workDir.append("stats"))
            .create());


    private final InMemoryValue<File> resultingKmerFilesPr = new InMemoryValue<File>();
    public final InValue<File> resultingKmerFiles =
            addOutput("resulting-kmers-file", resultingKmerFilesPr, File.class);


    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
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

        BigLong2ShortHashMap filter_hm = IOUtils.loadKmers(filterFiles.get(), maximalBadFrequency.get(),
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        for (File file : inputFiles.get()) {
            BigLong2ShortHashMap hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            String name = file.getName().replaceAll(".kmers.bin", "");
            File outFile = new File(outputDir.get(), name + ".kmers.bin");

            debug("Starting to print k-mers to " + outFile.getPath());
            long c = 0;
            try {
                c = IOUtils.filterAndPrintKmers(hm, filter_hm, maximalBadFrequency.get(),
                        maximalThreshold.get() * filterFiles.get().length, outFile);
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

    public KmersFilter() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new KmersFilter().mainImpl(args);
    }


}
