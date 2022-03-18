package tools;

import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongLongEntry;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.BoolParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import structures.ColoredKmers;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class ColorKmersMain  extends Tool {
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

    //input -- kmers <kmer, cnt>
    //param -- countonereadasone
    //param -- ncolor
    //results --  kmer color1p, colo2p, color3p ??? float or int???

    public static final String NAME = "kmers-color";

    public static final String DESCRIPTION = "Color kmers graph ";

    public static void main(String[] args) {
        new ColorKmersMain().mainImpl(args);
    }
    public ColorKmersMain() {
        super(NAME, DESCRIPTION);
    }
    protected ColorKmersMain(String name, String description) {
        super(name, description);
    }

    @Override
    protected void cleanImpl() {

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

    //todo нормирование на кол-во файлов одного кмера
    public static Map<String, Integer> readFileToColor(File file) throws ExecutionFailedException  {
        Map<String, Integer> res = new HashMap<String, Integer>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            br.readLine(); //first line is header: "row sampleID group"
            String line = br.readLine();
            while (line!=null) {
                String[] sline = line.split("\t");

                res.put(sline[1], groupToColorMap(sline[2]));
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
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        ColoredKmers coloredKmers = new ColoredKmers(colorsN.get(), availableProcessors.get());

        debug("Loading colors...");
        Map<String, Integer> fileToColorMap = readFileToColor(classFile.get());
        debug("Loading components...");
//        BigLong2ShortHashMap hm;
        debug("Memory used (before processing files) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);
        int cnt_loaded = 0;
        for (File kmersFile : kmersFiles.get()) {
            cnt_loaded+=1;
            hm.resetValues();
            System.out.println("kmers file #" + cnt_loaded + ": " + kmersFile.getName());
            hm = IOUtils.loadKmers(new File[]{kmersFile}, 0, availableProcessors.get(), logger);
            String compName = FileUtils.removeExtension(kmersFile.getName(), ".bin");
            String compName2 = FileUtils.removeExtension(compName, ".kmers");
            int color = fileToColorMap.get(compName2);
            Iterator<MutableLongShortEntry> iterator = hm.entryIterator();
            while (iterator.hasNext()) {
                MutableLongShortEntry curkmer = iterator.next();
                coloredKmers.addColor(curkmer.getKey(), color);
            }
            debug("Memory used (in cycle) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);

        }
        debug("Memory used (after cycle) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);
        hm.resetValues();
        debug("Memory used (after reset) = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);
        System.out.println(coloredKmers.kmersColors.size());
        System.out.println(coloredKmers.size);
        System.out.println(coloredKmers.colorsCNT);
        try {
            coloredKmers.saveColorInt(outputDir.get().getAbsolutePath());
            info("Components saved to " + outputDir.get());
        } catch (IOException e) {
            e.printStackTrace();
        }

//       try {
//            System.out.println("test");
//            ColoredKmers comp2 = new ColoredKmers(outputDir.get(), availableProcessors.get());
//            System.out.println(comp2.size);
//            System.out.println(comp2.colorsCNT);
//            System.out.println(comp2.kmersColors.size());
//            Integer[] colorsstat = new Integer[comp2.colorsCNT + 1];
//            Arrays.fill(colorsstat, 0);
//
//            Map<Long, Integer> col = comp2.getColors();
//            int cnt = 0;
//           Iterator<MutableLongLongEntry> iterator = comp2.kmersColors.entryIterator();
//           while (iterator.hasNext()) {
//               long kmer = iterator.next().getKey();
//               if (cnt%100==1) {
//                   System.out.println(kmer + " " + col.get(kmer) + " " + comp2.kmersColors.get(kmer) + " " + Arrays.toString(comp2.get_color_from_int(kmer)));
//               }
//               cnt+=1;
//               int c =  col.get(kmer);
//               colorsstat[c]+=1;
//           }
//           System.out.println(cnt);
//           for (int i = 0; i<comp2.colorsCNT + 1; i++) {
//               System.out.println(i + " = " + colorsstat[i]);
//           }
//            info("Loaded saved to " + outputDir.get());
//
//        } catch (ExecutionFailedException e) {
//            e.printStackTrace();
//        }
        debug("In the end used = " + Misc.usedMemoryAsString() + ", Time for preparing = " + t);
    }
}
