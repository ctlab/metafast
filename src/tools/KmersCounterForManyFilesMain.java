package tools;

import io.IOUtils;
import ru.ifmo.genetics.io.ReadersUtils;
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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class KmersCounterForManyFilesMain extends Tool {

    public static final String NAME = "kmer-counter-many";
    public static final String DESCRIPTION = "Count k-mers for many files independently";

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

    public File[] outputDescFiles = null;



    private final InMemoryValue<File[]> resultingKmerFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> resultingKmerFiles =
            addOutput("resulting-kmers-files", resultingKmerFilesPr, File[].class);


    private List<KmersCounterMain> counters = new ArrayList<KmersCounterMain>();

    private Timer t;

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        t = new Timer();
        counters.clear();

        File[] files = inputFiles.get();
        Arrays.sort(files);
        List<String> names = new ArrayList<String>();
        for (File f : files) {
            names.add(ReadersUtils.readDnaLazy(f).name());
        }

        int i=0;
        while (i < files.length) {
            if ((names.get(i).endsWith("_r1") && i+1<files.length && names.get(i+1).endsWith("_r2")) ||
                    (names.get(i).endsWith("_R1") && i+1<files.length && names.get(i+1).endsWith("_R2"))) {
                KmersCounterMain counter = new KmersCounterMain();
                counter.workDir.set(workDir.append("sub-counter"));
                counter.k.set(k);
                counter.inputFiles.set(new File[]{files[i], files[i+1]});
                counter.maximalBadFrequency.set(maximalBadFrequency);
                counter.outputDir.set(outputDir);
                counter.statsDir.set(statsDir);

                addStep(counter);
                counters.add(counter);
                i += 2;
            } else {
                KmersCounterMain counter = new KmersCounterMain();
                counter.workDir.set(workDir.append("sub-counter"));
                counter.k.set(k);
                counter.inputFiles.set(new File[]{files[i]});
                counter.maximalBadFrequency.set(maximalBadFrequency);
                counter.outputDir.set(outputDir);
                counter.statsDir.set(statsDir);

                addStep(counter);
                counters.add(counter);
                i += 1;
            }
        }

    }

    @Override
    protected void cleanImpl() {
        File[] outFiles = new File[counters.size()];
        for (int i = 0; i < counters.size(); i++) {
            outFiles[i] = counters.get(i).resultingKmerFiles.get();
        }
        resultingKmerFilesPr.set(outFiles);
        debug("Kmer-counter-many has finished! Time = " + t);
    }

    @Override
    protected void postprocessing() {
        IOUtils.tryToAppendDescription(outputDescFiles,
                statsDir.get(),
                "Directory with kmer frequency statistics " +
                        "(statistics files is in text format for every input reads file)"
        );
    }


    public static void main(String[] args) {
        new KmersCounterForManyFilesMain().mainImpl(args);
    }

    public KmersCounterForManyFilesMain() {
        super(NAME, DESCRIPTION);
    }

}
