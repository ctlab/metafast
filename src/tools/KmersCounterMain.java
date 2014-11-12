package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.statistics.QuantitativeStatistics;
import ru.ifmo.genetics.statistics.QuickQuantitativeStatistics;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.structures.map.BigLong2IntHashMap;
import ru.ifmo.genetics.structures.map.MutableLongIntEntry;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.*;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class KmersCounterMain extends Tool {

    public static final String NAME = "kmer-counter";

    public static final String DESCRIPTION = "Count k-mers in given reads with ArrayLong2IntHashMap";
//            "\nBinary output format: 64 bits to k-mer itself + 32 bits to frequency";

    static final int LOAD_TASK_SIZE = 1 << 15;  // 32 K reads

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


    /*
    To think about that params

    public final Parameter<Boolean> plainText = addParameter(new BoolParameterBuilder("plainText")
            .withShortOpt("pt")
            .optional()
            .withDescription("print k-mers in plain text, not in binary (64 + 32) format")
            .withDefaultValue(false)
            .create());

    public final Parameter<Boolean> notPrintNames = addParameter(new BoolParameterBuilder("not-print-names")
            .withShortOpt("npn")
            .optional()
            .withDescription("print only k-mers counts for all possible k-mers (k < 15)")
            .withDefaultValue(false)
            .create());
    */

    private final InMemoryValue<File> resultingKmerFilesPr = new InMemoryValue<File>();
    public final InValue<File> resultingKmerFiles =
            addOutput("resulting-kmers-file", resultingKmerFilesPr, File.class);


    @Override
    protected void runImpl() throws ExecutionFailedException {
        int LEN = k.get();

        Timer t = new Timer();
        BigLong2IntHashMap hm;
        try {
            hm = IOUtils.loadReads(inputFiles.get(), LEN, LOAD_TASK_SIZE,
                    new ShortKmerIteratorFactory(), availableProcessors.get(), this.logger);
        } catch (IOException e) {
            throw new ExecutionFailedException("Couldn't load k-mers", e);
        }
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
