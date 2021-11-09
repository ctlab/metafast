package tools;

import algo.ComponentsBuilderAroundPivot;
import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import structures.ConnectedComponent;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

/**
 * Created by -- on 05.02.2021.
 */
public class ComponentExtractorMain extends Tool {
    public static final String NAME = "component-extractor";
    public static final String DESCRIPTION = "Extract graph components from tangled graph based on pivot k-mers";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());


    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("i")
            .withDescription("list of input files with graph k-mers in binary format")
            .create());

    public final Parameter<File[]> pivotFiles = addParameter(new FileMVParameterBuilder("pivot")
            .mandatory()
            .withDescription("list of input files with pivot k-mers in binary format")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .withDescription("file to write found components to")
            .withDefaultValue(workDir.append("components.bin"))
            .create());

    /*public final Parameter<Integer> depth = addParameter(new IntParameterBuilder("depth")
            .optional()
            .withDescription("Depth of traversal from pivot k-mers")
            .withDefaultValue(1)
            .create());
    */

    private final InMemoryValue<File> componentsStatPr = new InMemoryValue<File>();


    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        debug("Loading graph from files...");
        BigLong2ShortHashMap hm = IOUtils.loadKmers(inputFiles.get(), 0, availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        debug("Loading graph from files...");
        BigLong2ShortHashMap pivot = IOUtils.loadKmers(pivotFiles.get(), 0, availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);



        info("Searching for components...");
        List<ConnectedComponent> components;
        try {
            String statFP = workDir + File.separator + "components-stat.txt";
            //components = IntelligentComponentsBuilderAroundPivot.splitStrategy(hm, pivot, k.get(), depth.get(), statFP, logger);
            components = ComponentsBuilderAroundPivot.splitStrategy(hm, k.get(), pivot, statFP, logger);

            componentsStatPr.set(new File(statFP));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        }
        info("Total " + NumUtils.groupDigits(components.size()) + " components were found");
        if (components.size() == 0) {
            warn("No components were extracted!");
        }

        try {
            ConnectedComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
            info("Components saved to " + componentsFile.get());
        } catch (IOException e) {
            e.printStackTrace();
        }
        debug("Components-extractor has finished! Time = " + t);
    }

    @Override
    protected void cleanImpl() {
    }


    public static void main(String[] args) {
        new ComponentExtractorMain().mainImpl(args);
    }

    public ComponentExtractorMain() {
        super(NAME, DESCRIPTION);
    }
}