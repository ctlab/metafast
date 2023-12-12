package tools;

import algo.BestKmersExtractor;
import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.*;


public class SubsetSpecificExtractor extends Tool {

    public static final String NAME = "subset-specific";
    public static final String DESCRIPTION = "Output subset of top most specific k-mers based on given statistical ranking";

    public final Parameter<File> inputFile = addParameter(new FileParameterBuilder("input-kmers")
            .mandatory()
            .withShortOpt("i")
            .withDescription("file with filtered k-mers in binary format")
            .create());

    public final Parameter<File> ranksFile = addParameter(new FileParameterBuilder("ranks-kmers")
            .mandatory()
            .withShortOpt("rk")
            .withDescription("file with k-mers ranks in binary format")
            .create());

    public final Parameter<Integer> NBest = addParameter(new IntParameterBuilder("num-kmers")
            .mandatory()
            .withShortOpt("n")
            .withDescription("number of most specific k-mers to be extracted")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());


    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        DataInputStream ranksStream = new DataInputStream(new BufferedInputStream(new FileInputStream(ranksFile.get()), 1 << 24));
        BigLong2ShortHashMap kmersFile = IOUtils.loadKmers(new File[]{inputFile.get()}, 0, availableProcessors.get(), logger);

        if (kmersFile.size() < NBest.get()){
            throw new ExecutionFailedException("Trying to extract more k-mers then present in input file!");
        }

        BestKmersExtractor extractor = new BestKmersExtractor(kmersFile, ranksStream);
        File file = new File(outputDir.get(), ranksFile.get().getName().split("\\.")[0].split("_ranks")[0] + "_top_" + NBest.get().toString() + ".kmers.bin");
        DataOutputStream outputStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file), 1 << 24));
        extractor.outTopNKmers(NBest.get(), outputStream);
        outputStream.close();

        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        info("subset-specific has finished! Time elapsed = " + t);
    }

    public SubsetSpecificExtractor() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new SubsetSpecificExtractor().mainImpl(args);
    }
}

