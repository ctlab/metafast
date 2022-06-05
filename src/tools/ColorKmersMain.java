package tools;

import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import structures.ColoredKmers;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class ColorKmersMain extends Tool {

    public static final String NAME = "kmers-color";
    public static final String DESCRIPTION = "Color kmers graph";
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());
    public final Parameter<Integer> colorsN = addParameter(new IntParameterBuilder("color-size")
            .mandatory()
            .withShortOpt("cs")
            .withDescription("num of colors")
            .create());
    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("kf")
            .withDescription("list of input files with k-mers in binary format")
            .create());
    public final Parameter<File> classFile = addParameter(new FileParameterBuilder("class-file")
            .mandatory()
            .withShortOpt("cf")
            .withDescription("file with classes for each input data")
            .create());
    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());
    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("colored-kmers"))
            .withDescription("Destination of resulting Color Kmers")
            .create());

    public ColorKmersMain() {
        super(NAME, DESCRIPTION);
    }

    protected ColorKmersMain(String name, String description) {
        super(name, description);
    }

    public static void main(String[] args) {
        new ColorKmersMain().mainImpl(args);
    }

    private static int groupToColorMap(String group) {
        if (group.equals("nonIBD")) {
            return 0;
        } else if (group.equals("UC")) {
            return 1;
        } else {
            return 2;
        }
    }

    public static Map<String, Integer> readFileToColor(File file) throws ExecutionFailedException {
        Map<String, Integer> res = new HashMap<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            br.readLine(); //first line is header: "row sampleID group"
            String line = br.readLine();
            while (line != null) {
                String[] splitLine = line.split("\t");

                res.put(splitLine[1], groupToColorMap(splitLine[2]));
                line = br.readLine();
            }
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Can't load components: file not found", e);
        } catch (IOException e) {
            throw new ExecutionFailedException("Can't load components: unknown IOException", e);
        }

        return res;
    }

    @Override
    protected void cleanImpl() {

    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        ColoredKmers coloredKmers = new ColoredKmers(colorsN.get(), availableProcessors.get());

        debug("Loading colors...");
        Map<String, Integer> fileToColorMap = readFileToColor(classFile.get());
        debug("Loading components...");
        debug("Memory used (before processing files) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        int cnt_loaded = 0;
        for (File kmersFile : kmersFiles.get()) {
            cnt_loaded += 1;
            hm.resetValues();
            System.out.println("kmers file #" + cnt_loaded + ": " + kmersFile.getName());
            hm = IOUtils.loadKmers(new File[]{kmersFile}, 0, availableProcessors.get(), logger);
            String compName = FileUtils.removeExtension(kmersFile.getName(), ".bin");
            String compName2 = FileUtils.removeExtension(compName, ".kmers");
            int color = fileToColorMap.get(compName2);
            System.out.println(kmersFile.getName() + " " + compName + " " + compName2 + " " + color);
            Iterator<MutableLongShortEntry> iterator = hm.entryIterator();
            while (iterator.hasNext()) {
                MutableLongShortEntry curKmer = iterator.next();
                coloredKmers.addColor(curKmer.getKey(), color);
            }
        }
        debug("Memory used (after cycle) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);
        hm.resetValues();
        System.out.println(coloredKmers.kmersColors.size());
        System.out.println(coloredKmers.size);
        System.out.println(coloredKmers.colorsCNT);
        try {
            coloredKmers.saveColorInt(outputDir.get().getAbsolutePath());
            info("Components saved to " + outputDir.get());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
