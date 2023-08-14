package tools;

import ru.ifmo.genetics.Runner;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.DoubleParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;

public class StatsFeaturesBuilderMain extends Tool {
    public static final String NAME = "stats-features";
    public static final String DESCRIPTION = "Builds statistically significant features for group of metagenomic samples";


    static {
        forceParameter.set(true);
        launchOptions.remove(forceParameter);
    }

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

    public final Parameter<Boolean> splitComponents = addParameter(new BoolParameterBuilder("split")
            .withDescription("Save each component in separate file?")
            .withDefaultValue(false)
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


    public final File[] outputDescFiles = new File[]{
            new File("output_description.txt"),
            workDir.append("output_description.txt").get()};



    // adding sub tools
    public KmersCounterPositiveNegative kmersCounterPosNeg = new KmersCounterPositiveNegative();
    {
        setFix(kmersCounterPosNeg.k, k);
        setFix(kmersCounterPosNeg.positiveFiles, positiveFiles);
        setFix(kmersCounterPosNeg.negativeFiles, negativeFiles);
        setFix(kmersCounterPosNeg.maximalBadFrequency, maximalBadFrequency);
        setFixDefault(kmersCounterPosNeg.outputDir);
        kmersCounterPosNeg.outputDescFiles = outputDescFiles;
        addSubTool(kmersCounterPosNeg);
    }

    public StatsKmersFinder statsKmers = new StatsKmersFinder();
    {
        setFix(statsKmers.Afiles, kmersCounterPosNeg.resultingPositiveKmerFiles);
        setFix(statsKmers.Bfiles, kmersCounterPosNeg.resultingNegativeKmerFiles);
        setFix(statsKmers.PValueChi2, PValueChi2);
        setFix(statsKmers.PValueMW, PValueMW);
        setFix(statsKmers.maximalBadFrequency, maximalBadFrequency);
        setFixDefault(statsKmers.outputDir);
        addSubTool(statsKmers);
    }


    public ComponentExtractorMain compExtractor = new ComponentExtractorMain();
    {
        setFix(compExtractor.k, k);
        setFix(compExtractor.inputFiles, kmersCounterPosNeg.resultingPositiveKmerFiles);
        setFix(compExtractor.pivotFiles, statsKmers.resultingKmerFiles);
        addSubTool(compExtractor);
    }

    public FeaturesCalculatorMain featuresCalculator = new FeaturesCalculatorMain();
    {
        setFix(featuresCalculator.k, k);
        setFix(featuresCalculator.kmersFiles, kmersCounterPosNeg.resultingPositiveKmerFiles);
        setFix(featuresCalculator.componentsFile, compExtractor.componentsFile);
        setFix(featuresCalculator.selectedKmers, statsKmers.resultingKmerFiles);
        addSubTool(featuresCalculator);
    }

    public ComponentsToSequences comp2seq = new ComponentsToSequences();
    {
        setFix(comp2seq.k, k);
        setFix(comp2seq.componentsFile, compExtractor.componentsFile);
        setFix(comp2seq.splitComponents, splitComponents);
        addSubTool(comp2seq);
    }


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


        outputDescFiles[1] = workDir.append("output_description.txt").get();    // updating workdir
        createOutputDescFiles();

        // running steps
        addStep(kmersCounterPosNeg);
        addStep(statsKmers);
        addStep(compExtractor);
        addStep(featuresCalculator);
        addStep(comp2seq);
    }

    private void createOutputDescFiles() {
        for (File f : outputDescFiles) {
            try {
                PrintWriter out = new PrintWriter(f);
                out.println("# Output files' description for run started at " +
                        new SimpleDateFormat("dd-MMM-yyyy (EEE) HH:mm:ss").format(startDate));
                out.println();
                for (File ff : getLogFiles()) {
                    out.println(ff);
                }
                out.println("   Identical files with run log");

                out.println();
                for (File ff : outputDescFiles) {
                    out.println(ff);
                }
                out.println("   Identical files with output files' description");
                out.close();
            } catch (FileNotFoundException e) {
                warn("Can't create file " + f + ", skipping");
                debug(e.getClass().getName() + ": " + e.getMessage());
            }
        }
    }


    @Override
    protected void cleanImpl() {
        debug("Statistically significant features builder has finished! Time = " + t);
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
        new StatsFeaturesBuilderMain().mainImpl(args);
    }

    public StatsFeaturesBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
