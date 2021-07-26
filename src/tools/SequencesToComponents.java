package tools;

import algo.ComponentFromSequence;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import structures.SequenceComponent;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by -- on 03.02.2020.
 */
public class SequencesToComponents extends Tool {

    public static final String NAME = "seq2comp";
    public static final String DESCRIPTION = "Transforms sequences to components";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
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
        List<SequenceComponent> components = Collections.synchronizedList(new ArrayList<SequenceComponent>());
        for (File f : sequencesFiles.get()) {
            info("Loading file " + f.getName() + "...");
            List<Dna> sequences = ReadersUtils.loadDnas(f);
            int comps = components.size();
            ExecutorService execService = Executors.newFixedThreadPool(availableProcessors.get());
            for (Dna dna : sequences) {
                execService.execute(new ComponentFromSequence(components, dna, k.get()));
            }
            execService.shutdown();
            try {
                while (!execService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS)){}
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            info(components.size() - comps + " components added");
        }

        // post processing...
        Tool.debug(logger, "ans.size = " + components.size());
        String statFP = workDir + File.separator + "components-stat.txt";
        PrintWriter statPW = new PrintWriter(statFP);
        statPW.println("# component.no\tcomponent.size\tcomponent.weight");
        for (int i = 0; i < components.size(); i++) {
            SequenceComponent comp = components.get(i);
            statPW.println((i + 1) + "\t" + comp.size + "\t" + comp.weight);
        }
        statPW.close();
        componentsStatPr.set(new File(statFP));

        info("Total " + NumUtils.groupDigits(components.size()) + " components were found");

        try {
            SequenceComponent.saveComponents(components, componentsFile.get().getAbsolutePath());
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