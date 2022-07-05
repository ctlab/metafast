package tools;

import algo.ColoredComponentsBuilder;
import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2LongHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import structures.SequenceComponent;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by -- on 28.06.2022.
 */
public class ColoredComponentMain extends Tool {
    public static final String NAME = "component-colored";
    public static final String DESCRIPTION = "Extract graph components from tangled graph based on k-mers coloring";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());


    public final Parameter<File[]> kmersFile = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("i")
            .withDescription("input file with colored k-mers in binary format")
            .create());

    public final Parameter<Integer> n_groups = addParameter(new IntParameterBuilder("n_groups")
            .optional()
            .withShortOpt("group")
            .withDefaultValue(3)
            .withDescription("number of classes (default: 3)")
            .create());

    public final Parameter<Boolean> separate = addParameter(new BoolParameterBuilder("separate")
            .optional()
            .withDefaultValue(false)
            .withDescription("use only color-specific k-mers in components (does not work in linear mode)")
            .create());

    public final Parameter<Boolean> linear = addParameter(new BoolParameterBuilder("linear")
            .optional()
            .withDefaultValue(false)
            .withDescription("choose best path on fork to create linear components")
            .create());

    public final Parameter<Integer> n_comps = addParameter(new IntParameterBuilder("n_comps")
            .optional()
            .withShortOpt("comp")
            .withDefaultValue(-1)
            .withDescription("select not more than X components for each class (default: -1, means all components)")
            .create());

    public final Parameter<Double> perc = addParameter(new DoubleParameterBuilder("perc")
            .optional()
            .withDefaultValue(0.9)
            .withDescription("relative abundance of k-mer in group to become color-specific")
            .create());



    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("colored-components"))
            .withDescription("Output directory")
            .create());

    @Override
    protected void runImpl() throws ExecutionFailedException {
        Timer t = new Timer();
        debug("Loading colored k-mers...");
        BigLong2LongHashMap hm = IOUtils.loadLongKmers(kmersFile.get(), k.get(), availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        info("Searching for colored components...");
        HashMap<Integer, List<SequenceComponent>> components;
        String statFP = workDir + File.separator + "components-stat.txt";
        try {
            components = ColoredComponentsBuilder.splitStrategy(hm, k.get(), statFP, logger, n_groups.get(), separate.get(), linear.get(), n_comps.get(), perc.get());
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }

        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }

        int total = 0;
        for (Map.Entry<Integer, List<SequenceComponent>> comp: components.entrySet()) {
            int sz = comp.getValue().size();
            total += sz;
            info(NumUtils.groupDigits(sz) + " components were found for class " + comp.getKey());
            File outFile = new File(outDir, "components_color_" + comp.getKey() + ".bin");
            try {
                SequenceComponent.saveComponents(comp.getValue(), outFile.getAbsolutePath());
            } catch (IOException e) {
                throw new ExecutionFailedException("Error while writing components to file: " + outFile);
            }
        }
        info("Total " + NumUtils.groupDigits(total) + " components were found");
        debug("Component-colored has finished! Time = " + t);
    }


    @Override
    protected void cleanImpl() {
    }


    public static void main(String[] args) {
        new ColoredComponentMain().mainImpl(args);
    }

    public ColoredComponentMain() {
        super(NAME, DESCRIPTION);
    }
}
