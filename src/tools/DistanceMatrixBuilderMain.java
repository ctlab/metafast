package tools;

import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.FilesFromOneFileYielder;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class DistanceMatrixBuilderMain extends Tool {
    public static final String NAME = "build-matrix";
    public static final String DESCRIPTION = "Builds the distance matrix for input sequences/kmers";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .optional()
            .withShortOpt("k")
            .withDefaultValue(31)
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
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(0)
            .create());

    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len")
            .withShortOpt("l")
            .withDescription("minimum sequence length to be added to a component")
            .withDefaultValue(100)
            .create());



    // running tools
    public KmersCounterMain kmersCounter = new KmersCounterMain();
    {
        setFix(kmersCounter.k, k);
        setFix(kmersCounter.inputFiles, inputFiles);
        setFix(kmersCounter.maximalBadFrequency, maximalBadFrequency);
        addSubTool(kmersCounter);
    }

    public SeqBuilderMain seqBuilder = new SeqBuilderMain();
    {
        setFix(seqBuilder.maximalBadFrequency, maximalBadFrequency);
        setFix(seqBuilder.k, k);
        setFix(seqBuilder.inputFiles, new FilesFromOneFileYielder(kmersCounter.resultingKmerFiles));
        setFix(seqBuilder.sequenceLen, minLen);
        setFix(seqBuilder.outputFile, (File) null); // setFixDefault with null default value
        addSubTool(seqBuilder);
    }

    public SeqMergerMain seqMerger = new SeqMergerMain();
    {
        setFix(seqMerger.k, k);
        setFix(seqMerger.sequencesFiles, new FilesFromOneFileYielder(seqBuilder.outputFileOut));
        setFix(seqMerger.minLen, minLen);
        addSubTool(seqMerger);
    }

    public ReadsFeaturesBuilderMain featuresBuilder = new ReadsFeaturesBuilderMain();
    {
        setFix(featuresBuilder.k, k);
        setFix(featuresBuilder.componentsFiles, new FilesFromOneFileYielder(seqMerger.componentsFileOut));
        setFix(featuresBuilder.readsFiles, inputFiles);
        setFix(featuresBuilder.threshold, kmersCounter.maximalBadFrequency);
        addSubTool(featuresBuilder);
    }

    public CalculateDistanceMatrixMain calcDist = new CalculateDistanceMatrixMain();
    {
        setFix(calcDist.featuresFiles, featuresBuilder.featuresFilesOut);
        setFixDefault(calcDist.separator);
        addSubTool(calcDist);
    }


    @Override
    protected void runImpl() throws ExecutionFailedException {
        addStep(kmersCounter);
        addStep(seqBuilder);
        addStep(seqMerger);
        addStep(featuresBuilder);
        addStep(calcDist);
    }


    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new DistanceMatrixBuilderMain().mainImpl(args);
    }

    public DistanceMatrixBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
