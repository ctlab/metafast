package tools;

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

import java.io.*;
import java.util.*;

public class ReadsFeaturesBuilderMain extends Tool {
    public static final String NAME = "features-builder";
    public static final String DESCRIPTION = "Features builder";

    static final int LOAD_TASK_SIZE = 1 << 15;

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> componentsFiles = addParameter(new FileMVParameterBuilder("components-files")
            .mandatory()
            .withShortOpt("m")
            .withDescription("files with connected components")
            .create());

    public final Parameter<File[]> readsFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .withDescription("FASTQ, BINQ, FASTA (ignored reads with 'N') reads")
            .create());

    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("kmers")
            .withShortOpt("ka")
            .withDescription("additional k-mers files in binary (long+int) format")
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

            String modelDir =
                    String.format(workDir + File.separator + "%03d-" + componentsFile.getName(), models.size());
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
                ArrayLong2IntHashMap readsHM;
                try {
                    readsHM = IOUtils.loadReads(new File[]{readsFile}, k.get(), LOAD_TASK_SIZE,
                            new ShortKmerIteratorFactory(), availableProcessors.get(), this.logger);
                } catch (IOException e) {
                    throw new ExecutionFailedException("Couldn't load k-mers from " + readsFile, e);
                }

                for (int i = 0; i < models.size(); i++) {
                    String vectorFP = modelsDirs.get(i) + File.separator + readsFile.getName() + ".vec";
                    try {
                        buildAndPrintVector(readsHM, models.get(i), vectorFP);
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                        return;
                    }
                    info("Features for " + readsFile + " printed to " + vectorFP);
                    featuresFiles[curFiles] = new File(vectorFP);
                    curFiles++;
                }
            }
        }

        if (kmersFiles.get() != null) {
            for (File kmersFile : kmersFiles.get()) {
                ArrayLong2IntHashMap kmersHM;
                try {
                    kmersHM = IOUtils.loadKmers(new File[]{kmersFile},
                            threshold.get(), availableProcessors.get(), this.logger);
                } catch (IOException e) {
                    throw new ExecutionFailedException("Couldn't load k-mers from " + kmersFile, e);
                }

                for (int i = 0; i < models.size(); i++) {
                    String vectorFP = modelsDirs.get(i) + File.separator + kmersFile.getName() + ".vec";
                    try {
                        buildAndPrintVector(kmersHM, models.get(i), vectorFP);
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }
                    info("Features for " + kmersFile + " printed to " + vectorFP);
                    featuresFiles[curFiles] = new File(vectorFP);
                    curFiles++;
                }
            }
        }

        featuresFilesPr.set(featuresFiles);
    }

    private static void buildAndPrintVector(ArrayLong2IntHashMap readsHM,
                                            List<ConnectedComponent> components,
                                            String vectorFP) throws FileNotFoundException {
        List<Long> vector = new ArrayList<Long>();

        for (ConnectedComponent component : components) {
            long kmersInComponent = 0;
            for (long kmer : component.kmers) {
                int value = readsHM.get(kmer);
                kmersInComponent += value;
            }
            vector.add(kmersInComponent);
        }

        PrintWriter vectorPW = new PrintWriter(vectorFP);
        for (long x : vector) {
            vectorPW.println(x);
        }
        vectorPW.close();
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new ReadsFeaturesBuilderMain().mainImpl(args);
    }

    public ReadsFeaturesBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
