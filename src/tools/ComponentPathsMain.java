package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.set.LongHashSet;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;
import structures.Sequence;

import java.io.*;
import java.util.*;

public class ComponentPathsMain extends Tool {
    public static final String NAME = "component-paths";
    public static final String DESCRIPTION = "Extracts paths in the components";

    public static final int MAX_PATHS_COUNT = (int) 1e6;


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .mandatory()
            .withShortOpt("cf")
            .withDescription("binary file with connected components " +
                    "(usually saved in workDir/component-cutter/components.bin)")
            .create());

    public final Parameter<File[]> sequenceFiles = addParameter(new FileMVParameterBuilder("seq")
            .mandatory()
            .withDescription("files with paths (sequences) " +
                    "(may use workDir/seq-builder-many/sequences/*.seq.fasta)")
            .create());

    public final Parameter<Integer[]> components = addParameter(new IntMVParameterBuilder("components")
            .mandatory()
            .withShortOpt("cm")
            .withDescription("components' numbers to print paths to")
            .withDefaultValue(new Integer[]{})
            .create());

    public final Parameter<Boolean> allComponents = addParameter(new BoolParameterBuilder("all-components")
            .important()
            .withShortOpt("a")
            .withDescription("print paths for all components in component file (--components option will be ignored)")
            .create());


    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-length")
            .important()
            .withShortOpt("l")
            .withDescription("minimum length of the path to be printed (in nucleotides)")
            .withDefaultValue(50)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .important()
            .withShortOpt("o")
            .withDefaultValue(workDir.append("paths"))
            .withDescription("Destination of resulting FASTA sequences")
            .create());



    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();

        if (!allComponents.get() && (components.get().length == 0)) {
            error("No components to process!!! Do you forget to set --all-components or --components n1,n2,...?");
            return;
        }

        debug("Loading components...");
        List<ConnectedComponent> allComps;

        allComps = ConnectedComponent.loadComponents(componentsFile.get());
        info(NumUtils.groupDigits(allComps.size()) + " components loaded from " + componentsFile.get());


        debug("Preparing...");
        int n = allComponents.get() ? allComps.size() : components.get().length;
        if (n == 0) {
            error("No components to process!!!");
            System.exit(0);
        }
        ConnectedComponent[] usedComps = new ConnectedComponent[n];
        LongHashSet[] compKmers = new LongHashSet[n];
        if (allComponents.get()) {
            allComps.toArray(usedComps);
        } else {
            for (int i = 0; i < n; i++) {
                usedComps[i] = allComps.get(components.get()[i]-1);
            }
        }
        for (int i = 0; i < n; i++) {
            compKmers[i] = new LongHashSet((int) (usedComps[i].size / 0.75f) +10);
            for (long kmer : usedComps[i].kmers) {
                compKmers[i].add(kmer);
            }
        }
        List<Sequence>[] ans = new List[n];
        for (int i = 0; i < n; i++) {
            ans[i] = new ArrayList<Sequence>();
        }
        int k = this.k.get();



        debug("Loading sequences and extracting paths...");
        for (File seqFile : sequenceFiles.get()) {
            info("Loading file " + seqFile.getName() + "...");

            NamedSource<Dna> reader = ReadersUtils.readDnaLazy(seqFile);
            Iterator<Dna> iterator = reader.iterator();
            while (iterator.hasNext()) {
                Dna dna = iterator.next();

                for (int i = 0; i < n; i++) {
                    // checking sequence dna on component i
                    int first = -1;
                    int cur = 0;
                    for (ShortKmer kmer : ShortKmer.kmersOf(dna, k)) {
                        if (compKmers[i].contains(kmer.toLong())) {
                            if (first == -1) {
                                first = cur;
                            } else {
                                // continue...
                            }
                        } else {
                            if (first != -1) {
                                checkAndAddPath(dna, first, cur, ans[i], usedComps[i]);
                                first = -1;
                            }
                        }
                        cur++;
                    }
                    if (first != -1) {
                        checkAndAddPath(dna, first, cur, ans[i], usedComps[i]);
                    }
                }
            }
        }

        for (int i = 0; i < n; i++) {
            if (ans[i].size() == MAX_PATHS_COUNT) {
                warn("Too many paths in component " + usedComps[i].no + ", " +
                        "keeping only first " + MAX_PATHS_COUNT + " of them!");
            }
        }

        info("Sorting...");
        for (List<Sequence> ansI : ans) {
            Collections.sort(ansI, new Comparator<Sequence>() {
                @Override
                public int compare(Sequence o1, Sequence o2) {
                    return -Integer.compare(o1.length(), o2.length());
                }
            });
        }

        info("Saving to files...");
        outputDir.get().mkdir();
        for (int i = 0; i < n; i++) {
            File file = new File(outputDir.get(), "component-" + usedComps[i].no + ".seq.fasta");
            try {
                Sequence.printSequences(ans[i], file);
            } catch (IOException e) {
                throw new RuntimeException("Can't write sequences to file " + file, e);
            }
        }
        info("Paths for " + n + " component(s) were saved in directory " + outputDir.get());
    }

    private void checkAndAddPath(Dna dna, int first, int cur, List<Sequence> ans, ConnectedComponent comp) {
        // possible adding from first to cur
        int len = cur - first - 1 + k.get();
        if (len >= minLen.get()) {
            double avgKmerWeight = comp.weight / (double) comp.size;
            Sequence seq = new Sequence(
                    dna.substring(first, first + len),
                    (int) Math.round(avgKmerWeight),
                    0,0 // minWeight and maxWeight are unknown
            );
            if (ans.size() < MAX_PATHS_COUNT) {
                ans.add(seq);
            }
        }
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new ComponentPathsMain().mainImpl(args);
    }

    public ComponentPathsMain() {
        super(NAME, DESCRIPTION);
    }
}
