package tools;

import algo.ColoredComponentBuilder;
import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongLongEntry;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.tool.parameters.ParameterDescription;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import structures.ColoredKmers;
import structures.ConnectedComponent;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class ComponentColoredCutter extends Tool {
    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File> coloredKmersFile = addParameter(new FileParameterBuilder("colored k-mers")
            .mandatory()
            .withShortOpt("ckf")
            .withDescription("list of input files with k-mers in binary format")
            .create());

    public final Parameter<File[]> kmersFiles = addParameter(new FileMVParameterBuilder("kmers")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files")
            .create());


    /*doesn't matter now, just count*/
    public final Parameter<Integer> minComponentSize = addParameter(new IntParameterBuilder("min-component-size")
            .important()
            .withShortOpt("b1")
            .withDescription("minimum component size in component-cutter (in k-mers)")
            .withDefaultValue(10)
            .withDescriptionShort("Minimal component size")
            .withDescriptionRu("Минимальный размер компоненты для использования")
            .withDescriptionRuShort("Минимальный размер компоненты")
            .create());

    public final Parameter<Double> minForGreedStart = addParameter(new ParameterBuilder<Double>(Double.class, "min-for-greed-start")
            .withShortOpt("mfg")
            .withDescription("max count of components for each color")
            .withDefaultValue(0.90)
            .create());

    /*doesn't matter now, just count*/
    public final Parameter<Integer> maxComponentSize = addParameter(new IntParameterBuilder("max-component-size")
            .important()
            .withShortOpt("b2")
            .withDescription("maximum component size in component-cutter (in k-mers)")
            .withDefaultValue(1000)
            .withDescriptionShort("Maximal component size")
            .withDescriptionRu("Максимальный размер компоненты для использования")
            .withDescriptionRuShort("Максимальный размер компоненты")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .withDescription("file to write found components to")
            .withDefaultValue(workDir.append("components.bin"))
            .create());

    public final Parameter<String> splitMode  = addParameter(new StringParameterBuilder("mode-for-split")
            .withDescription("split mode, one of:  SEPARATE, COMMON, CHOOSE_THE_BEST")
            .withShortOpt("sp_mode")
            .withDefaultValue("SEPARATE")
            .create());
    public final Parameter<String> startMode  = addParameter(new StringParameterBuilder("mode-for-start")
            .withDescription("split mode, one of:  RANDOM, BEST")
            .withShortOpt("st_mode")
            .withDefaultValue("RANDOM")
            .create());
    public final Parameter<String> bfsMode  = addParameter(new StringParameterBuilder("mode-for-bfs")
            .withDescription("split mode, one of:  ALL, BEST")
            .withShortOpt("bfs_mode")
            .withDefaultValue("ALL")
            .create());

    public final Parameter<Integer> componentsForColor  = addParameter(new IntParameterBuilder("components-for-color")
            .withDescription("max count of components for each color")
            .withShortOpt("cfc")
            .withDefaultValue(1000)
            .create());


    public static final String NAME = "component-cutter-color";

    public static final String DESCRIPTION = "Build graph components from colored kmers ";
    //input --  kmer color1p, colo2p, color3p, colornp
    //output --ConnectedCompomemt produced by new rule
    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("colored-cutted-compomemts"))
            .withDescription("Destination of resulting Color Kmers compomemts")
            .create());

    public static void main(String[] args) {
        new ComponentColoredCutter().mainImpl(args);
    }

    public ComponentColoredCutter() {
        super(NAME, DESCRIPTION);
    }

    protected ComponentColoredCutter(String name, String description) {
        super(name, description);
    }

    @Override
    protected void cleanImpl() {

    }

    public enum SPLIT_MODE {
        SEPARATE, COMMON
    }

    //DEEP от развилки по всем соседям, до следующей, подсчет, сколько нашего цвета на пути
    //Для deep только common автоматом
    public enum BFS_MODE {
        ALL, BEST, DEEP
    }


    //greed --  берем случайный и смотрим, насколько он хорош, например, >0.99 или 0.9 и ограничиваю также по компонентам каждого цвета
    public enum START_KMER_MODE {
        RANDOM, BEST, GREED
    }

    private final InMemoryValue<File> componentsStatPr = new InMemoryValue<File>();


    private static final boolean IS_TEST = false;

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        ColoredKmers coloredKmers = new ColoredKmers(coloredKmersFile.get(), availableProcessors.get());
        if (IS_TEST) {
            System.out.println("test");
            System.out.println(coloredKmers.size);
            System.out.println(coloredKmers.colorsCNT);
            System.out.println(coloredKmers.kmersColors.size());
            Integer[] colorsstat = new Integer[coloredKmers.colorsCNT + 1];
            Arrays.fill(colorsstat, 0);
            int cnt = 0;
            Iterator<MutableLongLongEntry> iterator = coloredKmers.kmersColors.entryIterator();
            while (iterator.hasNext()) {
                long kmer = iterator.next().getKey();
                cnt += 1;
                int c = coloredKmers.getColor(kmer);
                colorsstat[c] += 1;
            }
            System.out.println(cnt);
            for (int i = 0; i < coloredKmers.colorsCNT + 1; i++) {
                System.out.println(i + " = " + colorsstat[i]);
            }
        } else {
            BigLong2ShortHashMap hm = IOUtils.loadKmers(kmersFiles.get(), 0, availableProcessors.get(), logger);
            info("Searching for components...");
            List<ConnectedComponent> components;
            String statFP = workDir + File.separator + "components-stat-" +
                    minComponentSize.get() + "-" + maxComponentSize.get() + ".txt";
            SPLIT_MODE splitmode = SPLIT_MODE.valueOf(splitMode.get());
            START_KMER_MODE startmode = START_KMER_MODE.valueOf(startMode.get());
            BFS_MODE bfsmode = BFS_MODE.valueOf(bfsMode.get());
            System.out.println("CFC is: " +componentsForColor.get());
            components = ColoredComponentBuilder.splitStrategy(hm, coloredKmers, (int) k.get(), (int) minComponentSize.get(),
                    (int) maxComponentSize.get(), statFP, logger, availableProcessors.get(), splitmode, startmode, bfsmode, componentsForColor.get(), minForGreedStart.get());
            componentsStatPr.set(new File(statFP));
            System.out.println("in main, components size is: " + components.size());
            info("Total " + NumUtils.groupDigits(components.size()) + " components were found");
            if (components.size() == 0) {
                warn("No components were extracted! Perhaps you should decrease --min-component-size value");
            }
            try {
                ConnectedComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
                info("Components saved to " + componentsFile.get());
            } catch (IOException e) {
                e.printStackTrace();
            }
            debug("Components-cutter has finished! Time = " + t);
        }
//        components.an


    }
}
