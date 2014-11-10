package tools;

import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.structures.map.BigLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;
import io.IOUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.FileUtils;

import java.io.*;
import java.util.*;

public class FeaturesCalculatorMain extends Tool {
    public static final String NAME = "features-calculator";
    public static final String DESCRIPTION = "Calculates features values for input reads files";

    static final int LOAD_TASK_SIZE = 1 << 15;

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> componentsFiles = addParameter(new FileMVParameterBuilder("components-files")
            .mandatory()
            .withShortOpt("cm")
            .withDescription("files with connected components (one component is considered as one feature)")
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
            .withDescription("threshold for k-mers from <reads>")
            .withDefaultValue(0)
            .create());


    // output values
    private final InMemoryValue<File[]> featuresFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> featuresFilesOut = addOutput("features-files", featuresFilesPr, File[].class);


    @Override
    protected void runImpl() throws ExecutionFailedException {
        Timer t = new Timer();
        List<List<ConnectedComponent>> models = new ArrayList<List<ConnectedComponent>>();
        List<String> modelsDirs = new ArrayList<String>();

        for (File componentsFile : componentsFiles.get()) {
            try {
                List<ConnectedComponent> components =
                        ConnectedComponent.loadComponents(componentsFile.getAbsolutePath());
                models.add(components);
                info(components.size() + " components loaded from " + componentsFile);
            } catch (IOException e) {
                throw new ExecutionFailedException("Couldn't load components", e);
            }

            String compName = FileUtils.removeExtension(componentsFile.getName(), ".bin");
            String modelDir = String.format(workDir + File.separator + "%03d-" + compName, models.size());
            modelsDirs.add(modelDir);
            (new File(modelDir)).mkdir();
            debug("Dir created for vectors: " + modelDir);
        }

        int featuresFilesCount = (readsFiles.get() == null ? 0 : readsFiles.get().length)
                                  + (kmersFiles.get() == null ? 0 : kmersFiles.get().length);
        File[] featuresFiles = new File[featuresFilesCount];
        int curFiles = 0;

        if (readsFiles.get() != null) {
            for (File readsFile : readsFiles.get()) {
                BigLong2IntHashMap readsHM;
                try {
                    readsHM = IOUtils.loadReads(new File[]{readsFile}, k.get(), LOAD_TASK_SIZE,
                            new ShortKmerIteratorFactory(), availableProcessors.get(), this.logger);
                } catch (IOException e) {
                    throw new ExecutionFailedException("Couldn't load k-mers from " + readsFile, e);
                }

                for (int i = 0; i < models.size(); i++) {
                    File outFile = new File(modelsDirs.get(i),
                            ReadersUtils.readDnaLazy(readsFile).name() + ".vec");
                    try {
                        buildAndPrintVector(readsHM, models.get(i), outFile);
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                        return;
                    }
                    info("Features for " + readsFile + " printed to " + outFile);
                    featuresFiles[curFiles] = outFile;
                    curFiles++;
                }
            }
        }

        if (kmersFiles.get() != null && kmersFiles.get().length != 0) {
            for (File kmersFile : kmersFiles.get()) {
                BigLong2IntHashMap kmersHM;
                try {
                    kmersHM = IOUtils.loadKmers(new File[]{kmersFile},
                            threshold.get(), availableProcessors.get(), this.logger);
                } catch (IOException e) {
                    throw new ExecutionFailedException("Couldn't load k-mers from " + kmersFile, e);
                }

                for (int i = 0; i < models.size(); i++) {
                    File outFile = new File(modelsDirs.get(i),
                            ReadersUtils.readDnaLazy(kmersFile).name() + ".vec");
                    try {
                        buildAndPrintVector(kmersHM, models.get(i), outFile);
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }
                    info("Features for " + kmersFile + " printed to " + outFile);
                    featuresFiles[curFiles] = outFile;
                    curFiles++;
                }
            }
        }

        featuresFilesPr.set(featuresFiles);
        debug("Features-calculator has finished! Time = " + t);
    }

    private static void buildAndPrintVector(BigLong2IntHashMap readsHM,
                                            List<ConnectedComponent> components,
                                            File outFile) throws FileNotFoundException {
        List<Long> vector = new ArrayList<Long>();

        for (ConnectedComponent component : components) {
            long kmersInComponent = 0;
            for (long kmer : component.kmers) {
                int value = readsHM.getWithZero(kmer);
                kmersInComponent += value;
            }
            vector.add(kmersInComponent);
        }

        PrintWriter vectorPW = new PrintWriter(outFile);
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
