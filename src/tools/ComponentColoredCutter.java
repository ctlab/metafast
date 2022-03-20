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
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
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

    private final InMemoryValue<File> componentsStatPr = new InMemoryValue<File>();


    private static final boolean IS_TEST = true;

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
            components = ColoredComponentBuilder.splitStrategy(hm, coloredKmers, k.get(), minComponentSize.get(),
                    maxComponentSize.get(), statFP, logger, availableProcessors.get());
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
