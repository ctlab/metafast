package tools;

import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.FilesFromOneFileYielder;
import ru.ifmo.genetics.utils.tool.values.IfYielder;
import ru.ifmo.genetics.utils.tool.values.InValue;
import ru.ifmo.genetics.utils.tool.values.SimpleFixingInValue;

import java.io.*;

public class DistanceMatrixBuilderMain extends Tool {
    public static final String NAME = "matrix-builder";
    public static final String DESCRIPTION = "Builds the distance matrix for input sequences/kmers";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .optional()
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
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len")
            .withShortOpt("l")
            .withDescription("minimum sequence length to be added to a component (in nucleotides)")
            .withDefaultValue(100)
            .create());



    // running tools
    public KmersCounterForManyFilesMain kmersCounter = new KmersCounterForManyFilesMain();
    {
        setFix(kmersCounter.k, k);
        setFix(kmersCounter.inputFiles, inputFiles);
        setFix(kmersCounter.maximalBadFrequency,
                new IfYielder<Integer>(new InValue<Boolean>() {
                        @Override
                        public Boolean get() {
                            return maximalBadFrequency.get() >= 1;
                        }
                    },
                    new SimpleFixingInValue<Integer>(1),
                    maximalBadFrequency
                )
        );
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
        setFix(featuresCalculator.componentsFiles, new FilesFromOneFileYielder(compCutter.componentsFileOut));
        setFix(featuresCalculator.readsFiles, inputFiles);
        setFix(featuresCalculator.threshold, kmersCounter.maximalBadFrequency);
        setFixDefault(featuresCalculator.kmersFiles);
        addSubTool(featuresCalculator);
    }

    public DistanceMatrixCalculatorMain distMatrixCalculator = new DistanceMatrixCalculatorMain();
    {
        setFix(distMatrixCalculator.featuresFiles, featuresCalculator.featuresFilesOut);
        setFixDefault(distMatrixCalculator.separator);
        addSubTool(distMatrixCalculator);
    }


    @Override
    protected void runImpl() throws ExecutionFailedException {
        addStep(kmersCounter);
        addStep(seqBuilder);
        addStep(compCutter);
        addStep(featuresCalculator);
        addStep(distMatrixCalculator);
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
