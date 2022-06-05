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
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import structures.ColoredKmers;
import structures.ConnectedSetComponent;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class ComponentColoredCutter extends Tool {
    public static final String NAME = "component-cutter-color";
    public static final String DESCRIPTION = "Build graph components from colored kmers";
    private static final boolean IS_TEST = false;

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
    public final Parameter<Integer> minComponentSize = addParameter(new IntParameterBuilder("min-component-size")
            .important()
            .withShortOpt("b1")
            .withDescription("minimum component size in component-cutter (in k-mers)")
            .withDefaultValue(10)
            .withDescriptionShort("Minimal component size")
            .withDescriptionRu("Минимальный размер компоненты для использования")
            .withDescriptionRuShort("Минимальный размер компоненты")
            .create());
    public final Parameter<Double> minForGreedStart = addParameter(new ParameterBuilder<>(Double.class, "min-for-greed-start")
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
    public final Parameter<String> splitMode = addParameter(new StringParameterBuilder("mode-for-split")
            .withDescription("split mode, one of:  SEPARATE, COMMON")
            .withShortOpt("sp_mode")
            .withDefaultValue("SEPARATE")
            .create());
    public final Parameter<String> startMode = addParameter(new StringParameterBuilder("mode-for-start")
            .withDescription("split mode, one of:  RANDOM, BEST, GREED")
            .withShortOpt("st_mode")
            .withDefaultValue("RANDOM")
            .create());
    public final Parameter<String> bfsMode = addParameter(new StringParameterBuilder("mode-for-bfs")
            .withDescription("split mode, one of:  ALL, BEST, DEEP")
            .withShortOpt("bfs_mode")
            .withDefaultValue("ALL")
            .create());
    public final Parameter<String> resSizeMode = addParameter(new StringParameterBuilder("mode-for-result")
            .withDescription("sizes in result  SMALL, GOOD, BIG, ALL")
            .withShortOpt("res_mode")
            .withDefaultValue("GOOD")
            .create());
    public final Parameter<Integer> componentsForColor = addParameter(new IntParameterBuilder("components-for-color")
            .withDescription("max count of components for each color")
            .withShortOpt("cfc")
            .withDefaultValue(1000)
            .create());
    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("colored-cut-components"))
            .withDescription("Destination of resulting Color Kmers components")
            .create());
    private final InMemoryValue<File> componentsStatPr = new InMemoryValue<>();

    public ComponentColoredCutter() {
        super(NAME, DESCRIPTION);
    }

    protected ComponentColoredCutter(String name, String description) {
        super(name, description);
    }

    public static void main(String[] args) {
        new ComponentColoredCutter().mainImpl(args);
    }

    @Override
    protected void cleanImpl() {

    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        ColoredKmers coloredKmers = new ColoredKmers(coloredKmersFile.get(), availableProcessors.get());
        if (IS_TEST) {
            System.out.println("test");
            System.out.println(coloredKmers.size);
            System.out.println(coloredKmers.colorsCNT);
            System.out.println(coloredKmers.kmersColors.size());
            Integer[] colorsStat = new Integer[coloredKmers.colorsCNT + 1];
            Arrays.fill(colorsStat, 0);
            int cnt = 0;
            Iterator<MutableLongLongEntry> iterator = coloredKmers.kmersColors.entryIterator();
            while (iterator.hasNext()) {
                long kmer = iterator.next().getKey();
                cnt += 1;
                int c = coloredKmers.getColor(kmer);
                colorsStat[c] += 1;
            }
            System.out.println(cnt);
            for (int i = 0; i < coloredKmers.colorsCNT + 1; i++) {
                System.out.println(i + " = " + colorsStat[i]);
            }
        } else {
            BigLong2ShortHashMap hm = IOUtils.loadKmers(kmersFiles.get(), 0, availableProcessors.get(), logger);
            info("Searching for components...");
            List<ConnectedSetComponent> components;
            String statFP = workDir + File.separator + "components-stat-" +
                    minComponentSize.get() + "-" + maxComponentSize.get() + ".txt";
            SPLIT_MODE splitMode = SPLIT_MODE.valueOf(this.splitMode.get());
            START_KMER_MODE startMode = START_KMER_MODE.valueOf(this.startMode.get());
            BFS_MODE bfsMode = BFS_MODE.valueOf(this.bfsMode.get());
            COMPONENT_SIZES_MODE resMode = COMPONENT_SIZES_MODE.valueOf(resSizeMode.get());
            components = ColoredComponentBuilder.splitStrategy(hm, coloredKmers, k.get(), minComponentSize.get(),
                    maxComponentSize.get(), statFP, logger, availableProcessors.get(), splitMode, startMode, bfsMode, componentsForColor.get(), minForGreedStart.get(), resMode);
            componentsStatPr.set(new File(statFP));
            info("Total " + NumUtils.groupDigits(components.size()) + " components were found");
            if (components.size() == 0) {
                warn("No components were extracted! Perhaps you should decrease --min-component-size value");
            }
            try {
                ConnectedSetComponent.saveSetComponents(components, componentsFile.get().getAbsolutePath());
                info("Components saved to " + componentsFile.get());
            } catch (IOException e) {
                e.printStackTrace();
            }
            debug("Components-cutter has finished! Time = " + t);
        }
    }


    public enum COMPONENT_SIZES_MODE {
        SMALL, GOOD, BIG, ALL
    }

    public enum SPLIT_MODE {
        SEPARATE, COMMON
    }

    public enum BFS_MODE {
        ALL, BEST, DEEP
    }

    public enum START_KMER_MODE {
        RANDOM, BEST, GREED
    }
}
