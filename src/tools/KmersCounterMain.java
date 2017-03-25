package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.*;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.File;
import java.io.IOException;

public class KmersCounterMain extends Tool {

    public static final String NAME = "kmer-counter";

    public static final String DESCRIPTION = "Count k-mers in given reads with ArrayLong2IntHashMap";


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size (maximum 31 due to realization details)")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of reads files from single environment. FASTQ, BINQ, FASTA (ignored reads with 'N')")
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
        BigLong2ShortHashMap hm = IOUtils.loadReads(inputFiles.get(), k.get(), 0,
                availableProcessors.get(),  logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }
        File stDir = statsDir.get();
        if (!stDir.exists()) {
            stDir.mkdirs();
        }

        String name = ReadersUtils.readDnaLazy(inputFiles.get()[0]).name()
                            + (inputFiles.get().length > 1 ? "+" : "");
        File outFile = new File(outDir, name + ".kmers.bin");
        File stFile = new File(stDir, name + ".stat.txt");


        debug("Starting to print k-mers to " + outFile.getPath());
        long c = 0;
        try {
            c = IOUtils.printKmers(hm, maximalBadFrequency.get(), outFile, stFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        info(NumUtils.groupDigits(hm.size()) + " k-mers found, "
                + NumUtils.groupDigits(c) + " (" + String.format("%.1f", c * 100.0 / hm.size()) + "%) of them is good (not erroneous)");

        if (hm.size() == 0) {
            warn("No k-mers found in reads! Perhaps you reads file is empty or k-mer size is too big");
        } else if (c == 0 || c < (long) (hm.size() * 0.03)) {
            warn("Too few good k-mers were found! Perhaps you should decrease k-mer size or --maximal-bad-frequency value");
        }
        long allKmersNumber = (1L << (2*k.get())) / 2;  // (4^k)/2
        if (hm.size() == allKmersNumber) {
            warn("All possible k-mers were found in reads! Perhaps you should increase k-mer size");
        } else if (hm.size() >= (long) (allKmersNumber * 0.99)) {
            warn("Almost all possible k-mers were found in reads! Perhaps you should increase k-mer size");
        }

        info("Good k-mers printed to " + outFile.getPath());
        resultingKmerFilesPr.set(outFile);
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new KmersCounterMain().mainImpl(args);
    }

    public KmersCounterMain() {
        super(NAME, DESCRIPTION);
    }

}
