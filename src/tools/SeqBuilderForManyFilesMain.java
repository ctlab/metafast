package tools;

import algo.SequencesFinders;
import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.Sequence;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SeqBuilderForManyFilesMain extends Tool {
    public static final String NAME = "seq-builder-many";
    public static final String DESCRIPTION = "Metagenome De Bruijn graph analysis and sequences building for many files independently";


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files with k-mers in binary format")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .create());

    public final Parameter<Integer> bottomCutPercent = addParameter(new IntParameterBuilder("bottom-cut-percent")
            .optional()
            .withShortOpt("bp")
            .withDescription("k-mers percent to be assumed erroneous while building sequences in seq-builder (if specified, --maximal-bad-frequency wouldn't be used in seq-builder)")
            .create());

    public final Parameter<Integer> sequenceLen = addParameter(new IntParameterBuilder("sequence-len")
            .mandatory()
            .withShortOpt("l")
            .withDescription("minimal sequence length to be written to sequences.fasta")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("sequences"))
            .withDescription("Destination of resulting FASTA sequences")
            .create());

    public File[] outputDescFiles = null;


    // output values
    private final InMemoryValue<File[]> outputFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> outputFilesOut = addOutput("output-files", outputFilesPr, File[].class);


    private List<SeqBuilderMain> builders = new ArrayList<SeqBuilderMain>();

    private Timer t;

    @Override
    protected void runImpl() throws ExecutionFailedException {
        if (maximalBadFrequency.get() != null && bottomCutPercent.get() != null) {
            throw new IllegalArgumentException("-b and -bp can not be set both");
        }
        t = new Timer();

        for (File f : inputFiles.get()) {
            SeqBuilderMain builder = new SeqBuilderMain();
            builder.workDir.set(workDir.append("sub-builder"));
            builder.k.set(k);
            builder.inputFiles.set(new File[]{f});
            builder.maximalBadFrequency.set(maximalBadFrequency);
            builder.bottomCutPercent.set(bottomCutPercent);
            builder.sequenceLen.set(sequenceLen);
            builder.outputDir.set(outputDir);

            addStep(builder);
            builders.add(builder);
        }
    }

    @Override
    protected void cleanImpl() {
        File[] outFiles = new File[builders.size()];
        for (int i = 0; i < builders.size(); i++) {
            outFiles[i] = builders.get(i).outputFileOut.get();
        }
        outputFilesPr.set(outFiles);
        debug("Seq-builder-many has finished! Time = " + t);
    }

    @Override
    protected void postprocessing() {
        IOUtils.tryToAppendDescription(outputDescFiles,
                outputDir.get(),
                "Directory with FASTA files - paths from reads for every library"
        );
    }


    public static void main(String[] args) {
        new SeqBuilderForManyFilesMain().mainImpl(args);
    }

    public SeqBuilderForManyFilesMain() {
        super(NAME, DESCRIPTION);
    }
}
