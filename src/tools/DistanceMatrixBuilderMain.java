package tools;

import ru.ifmo.genetics.Runner;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.*;

import java.io.*;
import java.util.Arrays;

public class DistanceMatrixBuilderMain extends Tool {
    public static final String NAME = "matrix-builder";
    public static final String DESCRIPTION = "Builds the distance matrix for input sequences";


    static {
        forceParameter.set(true);
        launchOptions.remove(forceParameter);
    }

    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .important()
            .withShortOpt("k")
            .withDefaultValue(31)
            .withDescription("k-mer size (in nucleotides, maximum 31 due to realization details)")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of reads files from single environment. FASTQ, BINQ, FASTA")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .important()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len")
            .important()
            .withShortOpt("l")
            .withDescription("minimum sequence length to be added to a component (in nucleotides)")
            .withDefaultValue(100)
            .create());

    public final Parameter<File> matrixFile = addParameter(new FileParameterBuilder("matrix-file")
            .withDefaultValue(workDir.append("matrices").append("dist_matrix_$DT.txt"))
            .withDefaultComment("<workDir>/matrices/dist_matrix_<date>_<time>.txt")
            .withDescription("resulting distance matrix file")
            .create());



    // running tools
    public KmersCounterForManyFilesMain kmersCounter = new KmersCounterForManyFilesMain();
    {
        setFix(kmersCounter.k, k);
        setFix(kmersCounter.inputFiles, inputFiles);
        setFix(kmersCounter.maximalBadFrequency, maximalBadFrequency);
        setFixDefault(kmersCounter.outputDir);
        addSubTool(kmersCounter);
    }

    public SeqBuilderForManyFilesMain seqBuilder = new SeqBuilderForManyFilesMain();
    {
        setFix(seqBuilder.k, k);
        setFix(seqBuilder.inputFiles, kmersCounter.resultingKmerFiles);
        setFix(seqBuilder.maximalBadFrequency, maximalBadFrequency);
        setFix(seqBuilder.sequenceLen, minLen);
        setFixDefault(seqBuilder.outputDir);
        addSubTool(seqBuilder);
    }

    public ComponentCutterMain compCutter = new ComponentCutterMain();
    {
        setFix(compCutter.k, k);
        setFix(compCutter.sequencesFiles, seqBuilder.outputFilesOut);
        setFix(compCutter.minLen, minLen);
        setFixDefault(compCutter.componentsFile);
        addSubTool(compCutter);
    }

    public FeaturesCalculatorMain featuresCalculator = new FeaturesCalculatorMain();
    {
        setFix(featuresCalculator.k, k);
        setFix(featuresCalculator.componentsFile, compCutter.componentsFileOut);
        setFix(featuresCalculator.readsFiles, inputFiles);
        setFixDefault(featuresCalculator.kmersFiles);
        setFixDefault(featuresCalculator.threshold);
        addSubTool(featuresCalculator);
    }

    public DistanceMatrixCalculatorMain distMatrixCalculator = new DistanceMatrixCalculatorMain();
    {
        setFix(distMatrixCalculator.featuresFiles, featuresCalculator.featuresFilesOut);
        setFix(distMatrixCalculator.matrixFile, matrixFile);
        addSubTool(distMatrixCalculator);
    }

    public HeatMapMakerMain heatMapMaker = new HeatMapMakerMain();
    {
        setFix(heatMapMaker.matrixFile, distMatrixCalculator.matrixFile);
        addSubTool(heatMapMaker);
    }


    @Override
    protected void runImpl() throws ExecutionFailedException {
        info("Found " + inputFiles.get().length + " libraries to process");
        addStep(kmersCounter);
        addStep(seqBuilder);
        addStep(compCutter);
        addStep(featuresCalculator);
        addStep(distMatrixCalculator);
        addStep(heatMapMaker);
    }


    @Override
    protected void cleanImpl() {
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
        new DistanceMatrixBuilderMain().mainImpl(args);
    }

    public DistanceMatrixBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
