package tools;

import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
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
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Created by -- on 03.02.2020.
 */
public class SequencesToComponents extends Tool {

    public static final String NAME = "seqs-to-comp";
    public static final String DESCRIPTION = "Transforms sequences to components";

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

    public final Parameter<File[]> sequencesFiles = addParameter(new FileMVParameterBuilder("sequences")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files")
            .create());

    public final Parameter<File> componentsFile = addParameter(new FileParameterBuilder("components-file")
            .withDescription("file to write found components to")
            .withDefaultValue(workDir.append("components.bin"))
            .create());


    // output values
    private final InMemoryValue<File> componentsStatPr = new InMemoryValue<File>();

    @Override
    protected void cleanImpl() {
    }

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        Timer t = new Timer();
        debug("Loading sequences from files...");
        List<ConnectedComponent> components = new ArrayList<ConnectedComponent>();
        for (File f : sequencesFiles.get()) {
            ConnectedComponent comp = new ConnectedComponent();
            BigLong2ShortHashMap hm = IOUtils.loadReads(new File[]{f}, k.get(), minLen.get(),
                    availableProcessors.get(), logger);
            Iterator<MutableLongShortEntry> it = hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                short value = entry.getValue();
                comp.add(key, value);
            }
            components.add(comp);
        }

        Collections.sort(components);

        // post processing...
        Tool.debug(logger, "ans.size = " + components.size());
        String statFP = workDir + File.separator + "components-stat.txt";
        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.no\tcomponent.size\tcomponent.weight\tusedFreqThreshold");
        for (int i = 0; i < components.size(); i++) {
            ConnectedComponent comp = components.get(i);
            statPW.println((i + 1) + "\t" + comp.size + "\t" + comp.weight + "\t" + comp.usedFreqThreshold);
        }
        statPW.close();
        componentsStatPr.set(new File(statFP));

        info("Total " + NumUtils.groupDigits(components.size()) + " components were found");

        try {
            ConnectedComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
            info("Components saved to " + componentsFile.get());
        } catch (IOException e) {
            e.printStackTrace();
        }
        debug("Components-cutter has finished! Time = " + t);
    }

    public static void main(String[] args) {
        new SequencesToComponents().mainImpl(args);
    }

    public SequencesToComponents() {
        super(NAME, DESCRIPTION);
    }
}
