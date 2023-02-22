package tools;

import algo.ColoredKmerOperations;
import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class ColorKmersMain extends Tool {

    public static final String NAME = "kmers-color";
    public static final String DESCRIPTION = "Color k-mers based on occurrences in classes";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("kf")
            .withDescription("list of input files with k-mers in binary format")
            .create());

    public final Parameter<File> classFile = addParameter(new FileParameterBuilder("class")
            .mandatory()
            .withDescription("tab-separated file with two columns: sample_name [string] and class [0|1|2]")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    public final Parameter<Boolean> countValues = addParameter(new BoolParameterBuilder("val")
            .optional()
            .withDefaultValue(false)
            .withDescription("if SET count k-mer occurrence as total coverage in samples, otherwise as number of samples")
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("colored-kmers"))
            .withDescription("Output directory")
            .create());


    public static Map<String, Integer> readFileToColor(File file) throws ExecutionFailedException {
        Map<String, Integer> file2Color = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] splitLine = line.split("\t");
                file2Color.put(splitLine[0], Integer.parseInt(splitLine[1]));
            }
        } catch (IOException e) {
            throw new ExecutionFailedException("Cannot read file: " + file);
        }
        return file2Color;
    }


    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        BigLong2LongHashMap coloredKmers = new BigLong2LongHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);

        info("Loading class file...");
        Map<String, Integer> file2Color = readFileToColor(classFile.get());
        debug("Memory used (before processing files) = " + Misc.usedMemoryAsString() + ", Time = " + t);
        info("Loading kmers files...");
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        for (File kmersFile : kmersFiles.get()) {
            hm.resetValues();
            hm = IOUtils.loadKmers(new File[]{kmersFile}, maximalBadFrequency.get(), availableProcessors.get(), logger);
            int color = file2Color.get(FileUtils.removeExtension(kmersFile.getName(), ".kmers.bin"));

            Iterator<MutableLongShortEntry> it = hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                long newValue;
                if (countValues.get()) {
                    short value = entry.getValue();
                    newValue = ColoredKmerOperations.addValue(coloredKmers.getWithZero(key), color, value);
                } else {
                    newValue = ColoredKmerOperations.addValue(coloredKmers.getWithZero(key), color);
                }
                coloredKmers.put(key, newValue);
            }
        }
        debug("Memory used (after cycle) = " + Misc.usedMemoryAsString() + ", Time = " + t);
        hm.resetValues();

        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }

        File outFile = new File(outDir, "colored_kmers.kmers.bin");
        File stFile = new File(outDir, "colored_kmers.stat.txt");
        debug("Starting to print k-mers to " + outFile.getPath());
        long c;
        try {
            c = IOUtils.printKmers(coloredKmers, 0, outFile, stFile);
        } catch (IOException e) {
            throw new ExecutionFailedException("Cannot write results to: " + outFile);
        }
        info(NumUtils.groupDigits(c) + " colored k-mers printed to " + outFile.getPath());
    }

    public ColorKmersMain() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new ColorKmersMain().mainImpl(args);
    }

}
