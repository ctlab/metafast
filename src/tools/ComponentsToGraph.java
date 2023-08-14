package tools;

import algo.Comp2Graph;
import io.IOUtils;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import structures.ConnectedComponent;


import java.io.*;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class ComponentsToGraph extends Tool {
    public static final String NAME = "comp2graph";
    public static final String DESCRIPTION = "Transforms components in binary format to de Bruijn graph in GFA format";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .mandatory()
            .withShortOpt("cf")
            .withDescription("binary components file")
            .create());

    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .optional()
            .withShortOpt("i")
            .withDescription("list of input files with k-mers in binary format for graph coverage")
            .create());


    public final Parameter<Boolean> coverageMode = addParameter(new BoolParameterBuilder("coverage")
            .optional()
            .withShortOpt("cov")
            .withDefaultValue(false)
            .withDescription("if SET count graph coverage as total k-mer occurrences in samples, otherwise as number of samples with k-mer (works only with `-i` option)")
            .create());


    public final Parameter<File> graphFile = addParameter(new FileParameterBuilder("graph-file")
            .withDescription("file to write found components to")
            .withDefaultValue(workDir.append("components-graph.gfa"))
            .create());


    @Override
    protected void cleanImpl() {

    }

    @Override
    protected void runImpl() throws ExecutionFailedException {
        if (k.get() <= 0) {
            error("The size of k-mer must be at least 1.");
            System.exit(1);
        }
        if (k.get() > 31) {
            error("The size of k-mer must be no more than 31.");
            System.exit(1);
        }

        List<ConnectedComponent> components = ConnectedComponent.loadComponents(componentsFile.get());
        info(components.size() + " components loaded from " + componentsFile.get());

        ExecutorService execService = Executors.newFixedThreadPool(availableProcessors.get());

        BigLong2ShortHashMap hm = null;
        if (kmersFiles.get() != null) {
            hm = IOUtils.loadKmers(kmersFiles.get(), 0, availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString());
            if (!coverageMode.get()) {
                hm.resetValues();
                for (File file : kmersFiles.get()) {
                    BigLong2ShortHashMap filt_hm = IOUtils.loadKmers(new File[]{file}, 0, availableProcessors.get(), logger);
                    debug("Memory used = " + Misc.usedMemoryAsString());
                    Iterator<MutableLongShortEntry> it = filt_hm.entryIterator();
                    while (it.hasNext()) {
                        MutableLongShortEntry entry = it.next();
                        long key = entry.getKey();
                        hm.put(key, (short) (hm.getWithZero(key) + 1));
                    }
                }
            }
        }
        debug("Memory used = " + Misc.usedMemoryAsString());



        ByteArrayOutputStream[] compFiles = new ByteArrayOutputStream[components.size()];
        for (int icomp=0; icomp < components.size(); icomp++) {
            ConnectedComponent cmp = components.get(icomp);
            compFiles[icomp] = new ByteArrayOutputStream();
            execService.execute(new Comp2Graph(k.get(), icomp, cmp, hm, compFiles[icomp], logger));
        }

        execService.shutdown();
        try {
            execService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS);
        } catch (InterruptedException e) {
            throw new ExecutionFailedException("Error while transforming components to graph image: " + e);
        }

        try(OutputStream outputStream = new FileOutputStream(graphFile.get())) {
            for (ByteArrayOutputStream out: compFiles) {
                out.writeTo(outputStream);
            }
        } catch (IOException e) {
            throw new ExecutionFailedException("Error while writing GFA image to file: " + e);
        }

        info("Graph components saved to GFA format!");
    }



    public static void main(String[] args) {
        new ComponentsToGraph().mainImpl(args);
    }

    public ComponentsToGraph() {
        super(NAME, DESCRIPTION);
    }
}

