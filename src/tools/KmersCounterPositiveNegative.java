package tools;

import ru.ifmo.genetics.Runner;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;

import java.io.File;

public class KmersCounterPositiveNegative extends Tool {
    public static final String NAME = "kmer-counter-posneg";
    public static final String DESCRIPTION = "Count k-mers for files from two groups independently";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size (in nucleotides, maximum 31 due to realization details)")
            .withDescriptionShort("k-mer size")
            .withDescriptionRu("Длина k-мера при построении графа де Брейна")
            .withDescriptionRuShort("Длина k-мера")
            .create());

    public final Parameter<File[]> positiveFiles = addParameter(new FileMVParameterBuilder("positiveReads")
            .withShortOpt("pos")
            .mandatory()
            .withDescription("list of reads files from positive group. FASTQ, FASTA")
            .create());

    public final Parameter<File[]> negativeFiles = addParameter(new FileMVParameterBuilder("negativeReads")
            .withShortOpt("neg")
            .mandatory()
            .withDescription("list of reads files from negative group. FASTQ, FASTA")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .important()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .withDescriptionShort("Maximal bad frequency")
            .withDescriptionRu("Максимальная частота ошибочного k-мера")
            .withDescriptionRuShort("Максимальная частота ошибочного k-мера")
            .create());


    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers_posneg"))
            .create());



    // adding sub tools
    public KmersCounterForManyFilesMain kmersCounterPositive = new KmersCounterForManyFilesMain();
    {
        setFix(kmersCounterPositive.k, k);
        setFix(kmersCounterPositive.inputFiles, positiveFiles);
        setFix(kmersCounterPositive.maximalBadFrequency, maximalBadFrequency);
        setFix(kmersCounterPositive.outputDir, workDir.append("pos").append("kmers"));
        kmersCounterPositive.workDir = workDir.append("pos");
    }

    public KmersCounterForManyFilesMain kmersCounterNegative = new KmersCounterForManyFilesMain();
    {
        setFix(kmersCounterNegative.k, k);
        setFix(kmersCounterNegative.inputFiles, negativeFiles);
        setFix(kmersCounterNegative.maximalBadFrequency, maximalBadFrequency);
        setFix(kmersCounterNegative.outputDir, workDir.append("neg").append("kmers"));
        kmersCounterNegative.workDir = workDir.append("neg");
    }

    public File[] outputDescFiles = null;

    private final InMemoryValue<File[]> resultingPositiveKmerFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> resultingPositiveKmerFiles =
            addOutput("resulting-pos-kmers-files", resultingPositiveKmerFilesPr, File[].class);

    private final InMemoryValue<File[]> resultingNegativeKmerFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> resultingNegativeKmerFiles =
            addOutput("resulting-neg-kmers-files", resultingNegativeKmerFilesPr, File[].class);


    private Timer t;

    @Override
    protected void runImpl() throws ExecutionFailedException {
        // preparing
        t = new Timer();

        info("Found " + positiveFiles.get().length + " samples in positive class and "
                + negativeFiles.get().length + " samples in negative class to process");
        if (positiveFiles.get().length == 0 || negativeFiles.get().length == 0) {
            throw new ExecutionFailedException("No libraries to process!!! Can't continue the calculations.");
        }

        // running steps
        addStep(kmersCounterPositive);
        addStep(kmersCounterNegative);
    }


    @Override
    protected void cleanImpl() {
        resultingPositiveKmerFilesPr.set(kmersCounterPositive.resultingKmerFiles.get());
        resultingNegativeKmerFilesPr.set(kmersCounterNegative.resultingKmerFiles.get());
        debug("K-mers counter has finished! Time = " + t);
    }


    @Override
    public void mainImpl(String[] args) {
        if (Runner.containsOption(args, Runner.getOptKeys(continueParameter)) ||
                Runner.containsOption(args, Runner.getOptKeys(Tool.startParameter))) {
            forceParameter.set(false);
        }
        super.mainImpl(args);
    }

    public static void main(String[] args) {
        new KmersCounterPositiveNegative().mainImpl(args);
    }

    public KmersCounterPositiveNegative() {
        super(NAME, DESCRIPTION);
    }
}