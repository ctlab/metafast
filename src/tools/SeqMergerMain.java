package tools;

import structures.ConnectedComponent;
import algo.KmerOperations;
import io.IOUtils;
import it.unimi.dsi.fastutil.longs.Long2IntMap.*;
import it.unimi.dsi.fastutil.longs.LongArrayFIFOQueue;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;

import java.io.*;
import java.util.*;

public class SeqMergerMain extends Tool {
    public static final String NAME = "sequences-merger";
    public static final String DESCRIPTION = "Connected components builder";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len")
            .withShortOpt("l")
            .withDescription("minimum sequence length to be added")
            .withDefaultValue(100)
            .create());

    public final Parameter<Integer> minComponentSize = addParameter(new IntParameterBuilder("min-component-size")
            .withShortOpt("b1")
            .withDescription("minimum component size to be added")
            .withDefaultValue(1000)
            .create());

    public final Parameter<Integer> maxComponentSize = addParameter(new IntParameterBuilder("max-component-size")
            .withShortOpt("b2")
            .withDescription("maximum component size to be added")
            .withDefaultValue(10000)
            .create());

    public final Parameter<File[]> sequencesFiles = addParameter(new FileMVParameterBuilder("sequences")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .withDescription("file with found components")
            .withDefaultValue(workDir.append("components"))
            .create());

    @Override
    protected void runImpl() throws ExecutionFailedException {
        ArrayLong2IntHashMap hm =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4);
        try {
            IOUtils.addFASTASequences(sequencesFiles.get(), hm, k.get(), minLen.get(), this.logger);
        } catch (IOException e) {
            e.printStackTrace();
        }

        List<ConnectedComponent> components = null;
        try {
            String statFP = workDir + File.separator + "components-stat-" +
                    minComponentSize.get() + "-" + maxComponentSize.get();
            components = splitStrategy(hm, minComponentSize.get(), maxComponentSize.get(), statFP);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        info(components.size() + " components found");

        try {
            ConnectedComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
            info("Components printed to " + componentsFile.get());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private List<ConnectedComponent> splitStrategy(ArrayLong2IntHashMap hm,
                                                   int b1,
                                                   int b2,
                                                   String statFP) throws FileNotFoundException {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();

        PrintWriter statPW = new PrintWriter(statFP);
        for (int freqThreshold = 0; ; freqThreshold++) {
            List<ConnectedComponent> components = getComponents(hm, freqThreshold);
            if (components.size() == 0) {
                break;
            }
            for (ConnectedComponent comp : components) {
                if (comp.size() < b1) {
                    banComponent(hm, comp);
                } else if (comp.size() < b2) {
                    ans.add(comp);
                    statPW.println(comp.size() + " " + comp.getWeight() + " " + freqThreshold);
                    banComponent(hm, comp);
                }
            }
            debug("Freq = " + freqThreshold + ", components count = " + ans.size());
        }
        statPW.close();

        return ans;
    }

    private static void banComponent(ArrayLong2IntHashMap hm, ConnectedComponent component) {
        for (long kmer : component.kmers) {
            hm.add(kmer, -hm.get(kmer));
        }
    }

    private List<ConnectedComponent> getComponents(ArrayLong2IntHashMap hm, int freqThreshold) {
        List<ConnectedComponent> ans = new ArrayList<ConnectedComponent>();

        ArrayLong2IntHashMap processedKmers =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4);
        for (int i = 0; i < hm.hm.length; ++i) {
            for (Entry entry : hm.hm[i].long2IntEntrySet()) {
                long kmer = entry.getLongKey();
                if (processedKmers.get(kmer) > 0) {
                    continue;
                }

                int value = entry.getIntValue();
                if (value > freqThreshold) {
                    ConnectedComponent comp = getComponent(hm, k.get(), kmer, freqThreshold, processedKmers);
                    ans.add(comp);
                }
            }
        }

        return ans;
    }

    private static ConnectedComponent getComponent(ArrayLong2IntHashMap hm,
                                                   int kValue,
                                                   long kmer,
                                                   int freqThreshold,
                                                   ArrayLong2IntHashMap processedKmers) {
        if (hm.get(kmer) <= freqThreshold || processedKmers.get(kmer) > 0) {
            return null;
        }
        ConnectedComponent ans = new ConnectedComponent();

        long weight = hm.get(kmer);
        LongArrayFIFOQueue queue = new LongArrayFIFOQueue();

        queue.enqueue(kmer);
        processedKmers.add(kmer, 1);

        while (queue.size() > 0) {
            long kmerRepr = queue.dequeue();
            ans.add(kmerRepr);

            for (long neighbour : KmerOperations.possibleNeighbours(kmerRepr, kValue)) {
                int value = hm.get(neighbour);
                if (value <= freqThreshold || processedKmers.get(neighbour) > 0) {
                    continue;
                }
                weight += value;
                processedKmers.add(neighbour, 1);
                queue.enqueue(neighbour);
            }
        }

        ans.setWeight(weight);
        return ans;
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new SeqMergerMain().mainImpl(args);
    }

    public SeqMergerMain() {
        super(NAME, DESCRIPTION);
    }
}
