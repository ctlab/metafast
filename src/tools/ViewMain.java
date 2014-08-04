package tools;

import io.IOUtils;
import it.unimi.dsi.fastutil.longs.Long2IntMap;
import ru.ifmo.genetics.ToolTemplate;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class ViewMain extends Tool {
    public static final String NAME = "view";
    public static final String DESCRIPTION = "Converts different binary objects to text format";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .withShortOpt("k")
            .withDefaultValue(31)
            .withDescription("k-mer size, used while saving object")
            .create());

    public final Parameter<File> kmersFile = addParameter(new FileParameterBuilder("kmers-file")
            .withShortOpt("kf")
            .withDescription("binary file with kmers")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .withShortOpt("cf")
            .withDescription("binary components file")
            .create());


    public final Parameter<File> outputFile = addParameter(new FileParameterBuilder("output-file")
            .withShortOpt("o")
            .withDescription("file to print to")
            .withDefaultValue((File) null)
            .withDefaultComment("print to the screen")
            .create());



    @Override
    protected void runImpl() throws ExecutionFailedException {
        if (kmersFile.get() == null && componentsFile.get() == null) {
            logger.warn("No input file is selected  --->  no data to display!");
            return;
        }
        if (kmersFile.get() != null && componentsFile.get() != null) {
            logger.warn("Kmers file and components file are selected simultaneously!");
            logger.warn("Be careful, going to print all information in this case...");
        }


        PrintWriter out;
        try {
            out = (outputFile.get() != null) ? new PrintWriter(outputFile.get())
                    : new PrintWriter(System.out);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Couldn't open output file", e);
        }

        if (kmersFile.get() != null) {
            ArrayLong2IntHashMap kmersHM;
            try {
                kmersHM = IOUtils.loadKmers(new File[]{kmersFile.get()},
                        0, availableProcessors.get(), this.logger);
            } catch (IOException e) {
                throw new ExecutionFailedException("Couldn't load k-mers from " + kmersFile, e);
            }
            logger.info("Printing kmers...");

            out.println("Kmer\tCount");
            for (Long2IntMap map : kmersHM.hm) {
                for (Long2IntMap.Entry e : map.long2IntEntrySet()) {
                    out.println(new ShortKmer(e.getLongKey(), k.get()) + "\t" + e.getIntValue());
                }
            }
        }


        if (componentsFile.get() != null) {
            List<ConnectedComponent> components;
            try {
                components = ConnectedComponent.loadComponents(componentsFile.get().getPath());
                info(components.size() + " components loaded from " + componentsFile.get());
            } catch (IOException e) {
                throw new ExecutionFailedException("Couldn't load components", e);
            }
            logger.info("Printing components...");

            out.println(components.size() + " components:");

            for (int i = 0; i < components.size(); i++) {
                ConnectedComponent component = components.get(i);
                out.println("Component " + (i + 1) + ", size = " + component.size() + " kmers, " +
                        "weight = " + component.getWeight() +". Kmers:");

                for (long kmer : component.kmers) {
                    out.println(new ShortKmer(kmer, k.get()).toString());
                }
                out.println();
            }
        }


        out.close();
    }



    @Override
    protected void cleanImpl() {
    }

    public ViewMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new ToolTemplate().mainImpl(args);
    }
}
