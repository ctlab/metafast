package tools;

import algo.ComponentsBuilder;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import structures.ConnectedComponent;
import io.IOUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.NumUtils;

import java.io.*;
import java.util.List;

public class ComponentCutterMain extends Tool {
    public static final String NAME = "component-cutter";
    public static final String DESCRIPTION = "Build graph components from tangled graph";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<Integer> minLen = addParameter(new IntParameterBuilder("min-seq-len")
            .important()
            .withShortOpt("l")
            .withDescription("minimum sequence length to be added")
            .withDefaultValue(100)
            .create());

    public final Parameter<Integer> minComponentSize = addParameter(new IntParameterBuilder("min-component-size")
            .important()
            .withShortOpt("b1")
            .withDescription("minimum component size in component-cutter (in k-mers)")
            .withDefaultValue(1000)
            .create());

    public final Parameter<Integer> maxComponentSize = addParameter(new IntParameterBuilder("max-component-size")
            .important()
            .withShortOpt("b2")
            .withDescription("maximum component size in component-cutter (in k-mers)")
            .withDefaultValue(10000)
            .create());

    public final Parameter<File[]> sequencesFiles = addParameter(new FileMVParameterBuilder("sequences")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .withDescription("file to write found components to")
            .withDefaultValue(workDir.append("components.bin"))
            .create());

    public File[] outputDescFiles = null;


    // output values
    public final InValue<File> componentsFileOut = addOutput("components-file", componentsFile, File.class);
    private final InMemoryValue<File> componentsStatPr = new InMemoryValue<File>();
    public final InValue<File> componentsStatOut = addOutput("components-stat", componentsStatPr, File.class);


    @Override
    protected void runImpl() throws ExecutionFailedException {
        Timer t = new Timer();
        debug("Loading sequences from files...");
        BigLong2ShortHashMap hm = IOUtils.loadReads(sequencesFiles.get(), k.get(), minLen.get(),
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        info("Searching for components...");
        List<ConnectedComponent> components;
        try {
            String statFP = workDir + File.separator + "components-stat-" +
                    minComponentSize.get() + "-" + maxComponentSize.get() + ".txt";
            components = ComponentsBuilder.splitStrategy(hm, k.get(), minComponentSize.get(),
                    maxComponentSize.get(), statFP, logger, availableProcessors.get());

            componentsStatPr.set(new File(statFP));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }
        info("Total " + NumUtils.groupDigits(components.size()) + " components were found");

        try {
            ConnectedComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
            info("Components printed to " + componentsFile.get());
        } catch (IOException e) {
            e.printStackTrace();
        }
        debug("Components-cutter has finished! Time = " + t);
    }

    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void postprocessing() {
        if (outputDescFiles != null) {
            for (File f : outputDescFiles) {
                try {
                    PrintWriter out = new PrintWriter(new FileWriter(f, true));
                    out.println();
                    out.println(componentsStatOut.get());
                    out.println("   File with components' statistics (in text format)");
                    out.println();
                    out.println(componentsFileOut.get());
                    out.println("   File with extracted components (in binary format)");
                    out.close();
                } catch (IOException e) {
                    // does not matter
                }
            }
        }
    }


    public static void main(String[] args) {
        new ComponentCutterMain().mainImpl(args);
    }

    public ComponentCutterMain() {
        super(NAME, DESCRIPTION);
    }
}
