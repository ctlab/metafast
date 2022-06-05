package tools;

import io.IOUtils;
import ru.ifmo.genetics.ToolTemplate;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Iterator;

public class DoubleViewMain extends Tool {
    public static final String NAME = "double-view";
    public static final String DESCRIPTION = "Converts two binary k-mers files to single text file";


    // input parameters
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .important()
            .withShortOpt("k")
            .withDefaultValue(31)
            .withDescription("k-mer size, used while saving object")
            .create());

    public final Parameter<File> mgxFile = addParameter(new FileParameterBuilder("kmers-mgx")
            .important()
            .withShortOpt("mgx")
            .withDescription("first binary file with k-mers")
            .create());

    public final Parameter<File> mtxFile = addParameter(new FileParameterBuilder("kmers-mtx")
            .important()
            .withShortOpt("mtx")
            .withDescription("second binary file with k-mers")
            .create());


    public final Parameter<File> outputFile = addParameter(new FileParameterBuilder("output-file")
            .important()
            .withShortOpt("o")
            .withDescription("file to print to")
            .withDefaultValue((File) null)
            .withDefaultComment("print to the screen")
            .create());



    @Override
    protected void runImpl() throws ExecutionFailedException {
        int _k = k.get();

        PrintWriter out;
        try {
            out = (outputFile.get() != null) ? new PrintWriter(outputFile.get())
                    : new PrintWriter(System.out);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Couldn't open output file", e);
        }

        BigLong2ShortHashMap mtxHM = IOUtils.loadKmers(new File[]{mtxFile.get()}, 0, availableProcessors.get(), logger);
        BigLong2ShortHashMap mgxHM = IOUtils.loadKmers(new File[]{mgxFile.get()}, 0, availableProcessors.get(), logger);

        logger.info("Printing kmers...");
        out.println("Kmer\tmtx_count\tmgx_count");
        Iterator<MutableLongShortEntry> it = mtxHM.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            short value = entry.getValue();

            out.println(new ShortKmer(key, _k) + "\t" + value + "\t" + mgxHM.getWithZero(key));
        }

        out.close();
    }



    @Override
    protected void cleanImpl() {
    }

    public DoubleViewMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new ToolTemplate().mainImpl(args);
    }
}
