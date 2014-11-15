package tools;

import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;
import io.IOUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.FileUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class FeaturesCalculatorMain extends Tool {
    public static final String NAME = "features-calculator";
    public static final String DESCRIPTION = "Calculates features values for input reads files";


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .mandatory()
            .withShortOpt("cm")
            .withDescription("file with connected components (one component is considered as one feature)")
            .create());

    public final Parameter<File[]> readsFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .withDescription("FASTQ, BINQ, FASTA reads")
            .create());

    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("kmers")
            .withShortOpt("ka")
            .withDescription("additional k-mers files in binary (long+int) format")
            .withDefaultValue(new File[]{})
            .create());

    public final Parameter<Integer> threshold = addParameter(new IntParameterBuilder("threshold")
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(0)
            .create());


    // output values
    private final InMemoryValue<File[]> featuresFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> featuresFilesOut = addOutput("features-files", featuresFilesPr, File[].class);


    @Override
    protected void runImpl() throws ExecutionFailedException {
        Timer t = new Timer();

        debug("Loading components...");
        List<ConnectedComponent> components;

        try {
            components = ConnectedComponent.loadComponents(componentsFile.get());
            info(NumUtils.groupDigits(components.size()) + " components loaded from " + componentsFile.get());
        } catch (IOException e) {
            throw new ExecutionFailedException("Couldn't load components", e);
        }

        String compName = FileUtils.removeExtension(componentsFile.get().getName(), ".bin");
        File outDir = new File(workDir.get(), "vectors");
        outDir.mkdirs();
        debug(outDir + " directory was created for components file " + componentsFile.get().getName());


        // preparing
        long kmers = 0;
        for (ConnectedComponent component : components) {
            kmers += component.size();
        }
        BigLong2LongHashMap hm = new BigLong2LongHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        for (ConnectedComponent component : components) {
            for (long kmer : component.kmers) {
                hm.put(kmer, 0);
            }
        }
        debug("Memory used (before processing files) = " + Misc.usedMemoryAsString() + ", time = " + t);


        int featuresFilesCount = (readsFiles.get() == null ? 0 : readsFiles.get().length)
                                  + (kmersFiles.get() == null ? 0 : kmersFiles.get().length);
        File[] featuresFiles = new File[featuresFilesCount];
        int curFiles = 0;

        if (readsFiles.get() != null) {
            for (File readsFile : readsFiles.get()) {
                hm.resetValues();
                IOUtils.calculatePresenceForReads(new File[]{readsFile}, k.get(), hm,
                        availableProcessors.get(), logger);

                File outFile = new File(outDir, ReadersUtils.readDnaLazy(readsFile).name() + ".vec");
                buildAndPrintVector(components, hm, threshold.get(), outFile);
                info("Features for file " + readsFile.getName() + " printed to " + outFile);
                featuresFiles[curFiles] = outFile;
                curFiles++;
            }
        }

        if (kmersFiles.get() != null) {
            for (File kmersFile : kmersFiles.get()) {
                hm.resetValues();
                IOUtils.calculatePresenceForKmers(new File[]{kmersFile}, hm,
                        availableProcessors.get(), logger);

                File outFile = new File(outDir, FileUtils.removeExtension(kmersFile.getName()) + ".vec");
                buildAndPrintVector(components, hm, threshold.get(), outFile);
                info("Features for file " + kmersFile + " printed to " + outFile);
                featuresFiles[curFiles] = outFile;
                curFiles++;
            }
        }

        featuresFilesPr.set(featuresFiles);
        debug("Features-calculator has finished! Time = " + t);
    }

    private static void buildAndPrintVector(List<ConnectedComponent> components,
                                            BigLong2LongHashMap hm, int threshold, File outFile)
            throws ExecutionFailedException {
        List<Long> vector = new ArrayList<Long>();

        for (ConnectedComponent component : components) {
            long kmersInComponent = 0;
            for (long kmer : component.kmers) {
                long value = hm.getWithZero(kmer);
                if (value > threshold) {
                    kmersInComponent += value;
                }
            }
            vector.add(kmersInComponent);
        }

        PrintWriter vectorPW = null;
        try {
            vectorPW = new PrintWriter(outFile);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Can't create file " + outFile, e);
        }
        for (long x : vector) {
            vectorPW.println(x);
        }
        vectorPW.close();
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new FeaturesCalculatorMain().mainImpl(args);
    }

    public FeaturesCalculatorMain() {
        super(NAME, DESCRIPTION);
    }
}
