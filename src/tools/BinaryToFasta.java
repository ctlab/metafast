package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;


/**
 * Created by -- on 09.02.2021.
 */
public class BinaryToFasta extends Tool {
    public static final String NAME = "bin2fasta";
    public static final String DESCRIPTION = "Converts different binary objects to FASTA format";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .important()
            .withShortOpt("k")
            .withDefaultValue(31)
            .withDescription("k-mer size, used while saving object (NEED with ANY operation)")
            .create());

    public final Parameter<File> kmersFile = addParameter(new FileParameterBuilder("kmers-file")
            .important()
            .withShortOpt("kf")
            .withDescription("binary file with kmers")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .important()
            .withShortOpt("cf")
            .withDescription("binary components file")
            .create());

    public final Parameter<Boolean> splitComponents = addParameter(new BoolParameterBuilder("split")
            .withDescription("Save each component in separate file? Works only for components converter")
            .withDefaultValue(false)
            .create());


    public final Parameter<File> outputFile = addParameter(new FileParameterBuilder("output-file")
            .important()
            .withShortOpt("o")
            .withDescription("file prefix to print to")
            .withDefaultValue((File) null)
            .withDefaultComment("print to the screen")
            .create());


    private final InMemoryValue<File[]> resultingKmerFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> resultingKmerFiles =
            addOutput("resulting-kmers-files", resultingKmerFilesPr, File[].class);


    @Override
    protected void runImpl() throws ExecutionFailedException {
        if (outputFile.get() != null && outputFile.get().getParent() != null) {
            (new File(outputFile.get().getParent())).mkdirs();
        }
        if (kmersFile.get() == null && componentsFile.get() == null) {
            logger.warn("No input file is selected  --->  no data to display!");
            return;
        }
        if (kmersFile.get() != null && componentsFile.get() != null) {
            logger.warn("Kmers file and components file are selected simultaneously!");
            logger.warn("Be careful, going to print all information in this case...");
        }


        if (kmersFile.get() != null) {
            PrintWriter out;
            try {
                out = (outputFile.get() != null) ? new PrintWriter(outputFile.get() + ".fasta")
                        : new PrintWriter(System.out);
            } catch (FileNotFoundException e) {
                throw new ExecutionFailedException("Couldn't open output file", e);
            }

            BigLong2ShortHashMap kmersHM =
                    IOUtils.loadKmers(new File[]{kmersFile.get()}, 0, availableProcessors.get(), logger);

            info("Printing kmers...");
            Iterator<MutableLongShortEntry> it = kmersHM.entryIterator();
            int i = 1;
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                out.println(">" + i);
                out.println(new ShortKmer(entry.getKey(), k.get()));
                i++;
            }
            out.close();
        }


        if (componentsFile.get() != null) {
            List<ConnectedComponent> components =
                    ConnectedComponent.loadComponents(componentsFile.get());
            info(components.size() + " components loaded from " + componentsFile.get());
            info("Printing " + components.size() + " components...");

            File[] outFiles = new File[splitComponents.get() ? components.size() : 1];

            if (splitComponents.get()) {
                for (int i = 0; i < components.size(); i++) {
                    PrintWriter out;
                    try {
                        out = (outputFile.get() != null)
                                ? new PrintWriter(outputFile.get() + "_" + (i + 1) + ".fasta")
                                : new PrintWriter(System.out);
                    } catch (FileNotFoundException e) {
                        throw new ExecutionFailedException("Couldn't open output file", e);
                    }
                    outFiles[i] = new File(outputFile.get() + "_" + (i + 1) + ".fasta");

                    ConnectedComponent component = components.get(i);
                    info("Component " + (i + 1) + ", size = " + component.size + " kmers, " +
                            "weight = " + component.weight);

                    int j = 1;
                    for (long kmer : component.kmers) {
                        out.println(">" + j);
                        out.println(new ShortKmer(kmer, k.get()).toString());
                        j++;
                    }
                    out.close();
                }
            } else {
                PrintWriter out;
                try {
                    out = (outputFile.get() != null)
                            ? new PrintWriter(outputFile.get() + ".fasta")
                            : new PrintWriter(System.out);
                } catch (FileNotFoundException e) {
                    throw new ExecutionFailedException("Couldn't open output file", e);
                }
                outFiles[0] = new File(outputFile.get() + ".fasta");

                for (int i = 0; i < components.size(); i++) {
                    ConnectedComponent component = components.get(i);
                    info("Component " + (i + 1) + ", size = " + component.size + " kmers, " +
                            "weight = " + component.weight);

                    int j = 1;
                    for (long kmer : component.kmers) {
                        out.println(">" + (i + 1) + "_" + j);
                        out.println(new ShortKmer(kmer, k.get()).toString());
                        j++;
                    }
                }
                out.close();
            }

            resultingKmerFilesPr.set(outFiles);
        }

    }



    @Override
    protected void cleanImpl() {
    }

    public BinaryToFasta() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new BinaryToFasta().mainImpl(args);
    }
}