package tools;

import io.IOUtils;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;

import java.io.*;
import java.util.Arrays;
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
            .withDefaultValue(new File[]{})
            .create());

    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("kmers")
            .withShortOpt("ka")
            .withDescription("k-mers files in binary (long+int) format")
            .withDefaultValue(new File[]{})
            .create());


    public final Parameter<File[]> selectedKmers = addParameter(new FileMVParameterBuilder("selected")
            .optional()
            .withDescription("Use only selected k-mers to calculate features")
            .create());

    public final Parameter<Integer> threshold = addParameter(new IntParameterBuilder("threshold")
            //.withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(0)
            .create());

    public File[] outputDescFiles = null;


    // output values
    private final InMemoryValue<File[]> featuresFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> featuresFilesOut = addOutput("features-files", featuresFilesPr, File[].class);
    private final InMemoryValue<File> featuresDirPr = new InMemoryValue<File>();
    public final InValue<File> featuresDirOut = addOutput("features-dir", featuresDirPr, File.class);


    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();

        debug("Loading components...");
        List<ConnectedComponent> components =
                ConnectedComponent.loadComponents(componentsFile.get());
        info(NumUtils.groupDigits(components.size()) + " components loaded from " + componentsFile.get());

        if (components.size() == 0) {
            throw new ExecutionFailedException("No components were found in input files! Can't continue the calculations.");
        }

        String compName = FileUtils.removeExtension(componentsFile.get().getName(), ".bin");
        File outDir = new File(workDir.get(), "vectors");
        outDir.mkdirs();
        debug(outDir + " directory was created for components file " + componentsFile.get().getName());
        featuresDirPr.set(outDir);


        // preparing
        BigLong2LongHashMap hm = new BigLong2LongHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        for (ConnectedComponent component : components) {
            for (long kmer : component.kmers) {
                hm.put(kmer, 0);
            }
        }
        debug("Kmers in components = " + NumUtils.groupDigits(hm.size()));
        final long[] vector = new long[components.size()];
        final double[] breadth = new double[components.size()];
        debug("Memory used (before processing files) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);


        int featuresFilesCount = (readsFiles.get() == null ? 0 : readsFiles.get().length)
                                  + (kmersFiles.get() == null ? 0 : kmersFiles.get().length);
        File[] featuresFiles = new File[featuresFilesCount];
        int curFiles = 0;

        BigLong2ShortHashMap selected = null;
        if (selectedKmers.get() != null) {
            selected = IOUtils.loadKmers(selectedKmers.get(), 0, availableProcessors.get(), logger);
        }

        if (readsFiles.get() != null) {
            for (File readsFile : readsFiles.get()) {
                hm.resetValues();
                IOUtils.calculatePresenceForReads(new File[]{readsFile}, k.get(), hm,
                        availableProcessors.get(), logger);

                File outFile = new File(outDir, ReadersUtils.readDnaLazy(readsFile).name() + ".vec");
                File outBreadthFile = new File(outDir, ReadersUtils.readDnaLazy(readsFile).name() + ".breadth");
                buildAndPrintVector(components, hm, threshold.get(), selected, vector, breadth, outFile, outBreadthFile);
                info("Features for file " + readsFile.getName() + " printed to " + outFile);
                info("Components breadth coverage for file " + readsFile.getName() + " printed to " + outBreadthFile);
                featuresFiles[curFiles] = outFile;
                curFiles++;
            }
        }

        if (kmersFiles.get() != null) {
            for (File kmersFile : kmersFiles.get()) {
                /*
                // normalize on total amount of k-mers in kmersFile
                BigLong2ShortHashMap hm_tmp =IOUtils.loadKmers(new File[]{kmersFile}, threshold.get(), availableProcessors.get(), logger);
                debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
                long totalKmers = 0;
                Iterator<MutableLongShortEntry> it = hm_tmp.entryIterator();
                while (it.hasNext()) {
                    MutableLongShortEntry entry = it.next();
                    int value = entry.getValue();
                    totalKmers += value;
                }
                */

                hm.resetValues();
                IOUtils.calculatePresenceForKmers(new File[]{kmersFile}, hm,
                        availableProcessors.get(), logger);

                File outFile = new File(outDir, FileUtils.removeExtension(kmersFile.getName(), ".kmers.bin") + ".vec");
                File outBreadthFile = new File(outDir, FileUtils.removeExtension(kmersFile.getName(), ".kmers.bin") + ".breadth");
                buildAndPrintVector(components, hm, threshold.get(), selected, vector, breadth, outFile, outBreadthFile);
                info("Features for file " + kmersFile.getName() + " printed to " + outFile);
                info("Components breadth coverage for file " + kmersFile.getName() + " printed to " + outBreadthFile);
                featuresFiles[curFiles] = outFile;
                curFiles++;
            }
        }

        featuresFilesPr.set(featuresFiles);
        debug("Features-calculator has finished! Time = " + t);
    }

    private void buildAndPrintVector(final List<ConnectedComponent> components, final BigLong2LongHashMap hm,
                                     final int threshold, final BigLong2ShortHashMap selected, final long[] vector, final double[] breadth, File outFile,
                                     File outBreadthFile) throws ExecutionFailedException {

        try {
            debug("Building vector...");

            // calculating
            Arrays.fill(vector, 0);
            Arrays.fill(breadth, 0);
            Thread[] workers = new Thread[availableProcessors.get()];
            int compPerWorker = (components.size() / workers.length) + ((components.size() % workers.length == 0) ? 0 : 1);
            for (int i = 0; i < workers.length; i++) {
                final int from = i * compPerWorker;
                final int to = Math.min(components.size(), (i + 1) * compPerWorker);
                if (to > from) {
                    workers[i] = new Thread(new Runnable() {
                        @Override
                        public void run() {
                            for (int i = from; i < to; i++) {
                                ConnectedComponent component = components.get(i);
                                long kmers = 0;
                                long kmersCount = 0, kmersFound = 0;
                                for (long kmer : component.kmers) {
                                    if (selected == null || selected.getWithZero(kmer) > 0) {
                                        long value = hm.getWithZero(kmer);
                                        if (value > threshold) {
                                            kmers += value;
                                            kmersFound++;
                                        }
                                        kmersCount++;
                                    }
                                }
                                vector[i] = kmers;
                                breadth[i] = ((double) kmersFound) / kmersCount;
                            }
                        }
                    });
                    workers[i].start();
                }
            }
            for (Thread thread : workers) {
                if (thread != null) {
                    thread.join();
                }
            }

            // writing to file
            PrintWriter out = new PrintWriter(new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(outFile)), 1 << 20));   // 1 Mb buffer
            for (long kmers : vector) {
                // out.println((double) kmers / totalKmers);
                out.println(kmers);
            }
            out.close();

            PrintWriter bout = new PrintWriter(new BufferedWriter(
                    new OutputStreamWriter(new FileOutputStream(outBreadthFile)), 1 << 20));   // 1 Mb buffer
            for (double kmers : breadth) {
                bout.println(kmers);
            }
            bout.close();
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't write vector to file " + outFile, e);
        } catch (InterruptedException e) {
            throw new ExecutionFailedException("Calculating thread was interrupted!", e);
        }
    }

    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void postprocessing() {
        IOUtils.tryToAppendDescription(outputDescFiles,
                featuresDirOut.get(),
                "Directory with features values files for every library (in text format)"
        );
    }


    public static void main(String[] args) {
        new FeaturesCalculatorMain().mainImpl(args);
    }

    public FeaturesCalculatorMain() {
        super(NAME, DESCRIPTION);
    }
}
