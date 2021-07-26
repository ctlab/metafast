package tools;

import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;

/**
 * Created by -- on 20.07.2021.
 */
public class ComponentsToSequences extends Tool {
    public static final String NAME = "comp2seq";
    public static final String DESCRIPTION = "Transforms components in binary format to FASTA sequences (contigs)";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .important()
            .withShortOpt("k")
            .withDefaultValue(31)
            .withDescription("k-mer size, used while saving object (NEED with ANY operation)")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .important()
            .withShortOpt("cf")
            .withDescription("binary components file")
            .create());

    public final Parameter<Boolean> splitComponents = addParameter(new BoolParameterBuilder("split")
            .withDescription("Save each component in separate file? Works only for components converter")
            .withDefaultValue(false)
            .create());


    //adding sub tools
    public BinaryToFasta bin2fasta = new BinaryToFasta();
    {
        setFix(bin2fasta.k, k);
        setFix(bin2fasta.componentsFile, componentsFile);
        setFix(bin2fasta.splitComponents, splitComponents);
        setFix(bin2fasta.outputFile, workDir.append("kmers_fasta").append("component"));
        addSubTool(bin2fasta);
    }

    public KmersCounterForManyFilesMain kmersCounter = new KmersCounterForManyFilesMain();
    {
        setFix(kmersCounter.k, k);
        setFix(kmersCounter.inputFiles, bin2fasta.resultingKmerFiles);
        setFix(kmersCounter.maximalBadFrequency, 0);
        setFixDefault(kmersCounter.outputDir);
        addSubTool(kmersCounter);
    }

    public SeqBuilderForManyFilesMain seqBuilder = new SeqBuilderForManyFilesMain();
    {
        setFix(seqBuilder.k, k);
        setFix(seqBuilder.inputFiles, kmersCounter.resultingKmerFiles);
        setFix(seqBuilder.maximalBadFrequency, 0);
        setFix(seqBuilder.sequenceLen, k);
        setFixDefault(seqBuilder.outputDir);
        addSubTool(seqBuilder);
    }


    @Override
    protected void runImpl() throws ExecutionFailedException {
        // running steps
        addStep(bin2fasta);
        addStep(kmersCounter);
        addStep(seqBuilder);
    }

    @Override
    protected void cleanImpl() { }

    public static void main(String[] args) {
        new ComponentsToSequences().mainImpl(args);
    }

    public ComponentsToSequences() {
        super(NAME, DESCRIPTION);
    }
}