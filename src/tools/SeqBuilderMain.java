package tools;

import algo.SequencesFinders;
import io.IOUtils;
import ru.ifmo.genetics.statistics.*;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;
import ru.ifmo.genetics.utils.FileUtils;
import ru.ifmo.genetics.utils.NumUtils;
import structures.Sequence;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;

public class SeqBuilderMain extends Tool {
    public static final String NAME = "seq-builder";
    public static final String DESCRIPTION = "Metagenome De Bruijn graph analysis and sequences building";

    static final int STAT_LEN = 1024;

    static final String SEQUENCES_FILENAME = "sequences.fasta";
    static final String DISTRIBUTION_FILENAME = "distribution";

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("i")
            .withDescription("list of input files with k-mers in binary format")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequency")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> bottomCutPercent = addParameter(new IntParameterBuilder("bottom-cut-percent")
            .optional()
            .withShortOpt("bp")
            .withDescription("k-mers percent to be assumed erroneous (if set, maximal-bad-frequency value isn't used)")
            .create());

    public final Parameter<Integer> sequenceLen = addParameter(new IntParameterBuilder("sequence-len")
            .mandatory()
            .withShortOpt("l")
            .withDescription("sequence minimal length to be written to " + SEQUENCES_FILENAME)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withShortOpt("o")
            .withDefaultValue(workDir.append("sequences"))
            .withDescription("Destination of resulting FASTA sequences")
            .create());


    // output values
    private final InMemoryValue<File> outputFilePr = new InMemoryValue<File>();
    public final InValue<File> outputFileOut = addOutput("output-file", outputFilePr, File.class);


    @Override
    protected void runImpl() throws ExecutionFailedException {
        Timer t = new Timer();
        BigLong2ShortHashMap hm =
                IOUtils.loadKmers(inputFiles.get(), maximalBadFrequency.get(), availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

        long totalKmers = 0;
        int[] stat = new int[STAT_LEN];
        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            int value = entry.getValue();
            totalKmers += value;
            if (value >= stat.length) {
                value = stat.length - 1;
            }
            ++stat[value];
        }

        try {
            dumpStat(stat, workDir + File.separator + DISTRIBUTION_FILENAME);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException(e);
        }

        if (bottomCutPercent.get() != null) {
            info("Using bottom cut percent = " + bottomCutPercent.get());
            long kmersToCut = totalKmers * bottomCutPercent.get() / 100;
            debug("K-mers under given threshold = " + NumUtils.groupDigits(kmersToCut));
            long currentKmersCount = 0;
            for (int i = 0; i < stat.length - 1; i++) {
                if (currentKmersCount >= kmersToCut) {
                    maximalBadFrequency.set(i);
                    break;
                }
                currentKmersCount += (long) i * stat[i];
            }
        }

        // currently not used
//        if (maximalBadFrequency.get() == null) {
//            int threshold = 1;
//            long currentSum = 0;
//            while (stat[threshold] * (long) threshold > stat[threshold + 1] * (long) (threshold + 1)) {
//                currentSum += stat[threshold];
//                if (currentSum * 2 > totalKmers) {
//                    debug("Threshold search stopped at 50 %");
//                    break;
//                }
//                threshold++;
//            }
//            maximalBadFrequency.set(threshold);
//        }

        info("Using maximal bad frequency = " + maximalBadFrequency.get());

        File dir = outputDir.get();
        if (!dir.isDirectory()) {
            dir.mkdir();
        }
        String basename = FileUtils.removeExtension(inputFiles.get()[0].getName(), ".kmers.bin");
        String fp = dir + File.separator + basename;
        fp += (inputFiles.get().length > 1 ? "+" : "") + ".seq.fasta";

        File destination = new File(fp);
        outputFilePr.set(destination);

        Deque<Sequence> sequences;
        try {
            sequences = SequencesFinders.thresholdStrategy(hm, availableProcessors.get(),
                    maximalBadFrequency.get(), sequenceLen.get(), k.get());
        } catch (InterruptedException e) {
            e.printStackTrace();
            return;
        }
        info(NumUtils.groupDigits(sequences.size()) + " sequences found");
        if (sequences.size() == 0) {
            warn("No sequences were found! Perhaps you should decrease --min-seq-len or --maximal-bad-frequency values");
        }
        debug("Memory used (without running GC) = " + Misc.usedMemoryWithoutRunningGCAsString());

        try {
            Sequence.printSequences(sequences, destination);
        } catch (IOException e) {
            throw new RuntimeException("Can't write sequences to file", e);
        }
        info("Sequences printed to " + destination);

        //info("N50 value of sequences = " + getN50(sequenceLen));
        //dumpSeqInfo(sequenceLen, sequenceWeight, workDir + File.separator + "seq-info");
    }

    void dumpStat(int[] stat, String filename) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(filename);
        for (int i = 1; i < stat.length; ++i) {
            pw.println(i + " " + stat[i]);
        }
        pw.close();
    }

    void dumpSeqInfo(ArrayList<Integer> lens, ArrayList<Long> weights, String filename) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(filename);
        for (int i = 1; i < lens.size(); ++i) {
            pw.println(lens.get(i) + " " + weights.get(i));
        }
        pw.close();
    }

    int getN50(ArrayList<Integer> lens) {
        ArrayList<Integer> sorted = new ArrayList<Integer>(lens);
        Collections.sort(sorted);
        long sum = 0;
        for (int x : sorted) {
            sum += x;
        }
        long topSum = 0;
        for (int i = sorted.size() - 1; i >= 0; i--) {
            topSum += sorted.get(i);
            if (topSum * 2 >= sum) {
                return sorted.get(i);
            }
        }
        return -1;
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new SeqBuilderMain().mainImpl(args);
    }

    public SeqBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
