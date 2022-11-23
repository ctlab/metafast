package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Iterator;

/**
 * Created by -- on 22.11.2022.
 */
public class KmersPerSampleCounter  extends Tool {

    public static final String NAME = "kmers-per-sample";

    public static final String DESCRIPTION = "Counts the abundance of frequent k-mers from dataset in each sample";


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size (maximum 31 due to realization details)")
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("k-mers")
            .mandatory()
            .withShortOpt("i")
            .withDescription("list of input files with k-mers in binary format")
            .create());

    public final Parameter<Integer> percentPresent = addParameter(new IntParameterBuilder("percent-present")
            .optional()
            .withShortOpt("perc")
            .withDescription("Output only k-mers present in at least '-perc'% of samples (Default: 20)")
            .withDefaultValue(20)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers_per_samples"))
            .create());


    @Override
    protected void runImpl() throws ExecutionFailedException {
        if (k.get() <= 0) {
            error("The size of k-mer must be at least 1.");
            System.exit(1);
        }
        if (k.get() > 31) {
            error("The size of k-mer must be no more than 31.");
            System.exit(1);
        }


        Timer t = new Timer();
        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }


        info("Loading all k-mers...");
        BigLong2ShortHashMap hm = null;
        for (File file : inputFiles.get()) {
            BigLong2ShortHashMap filt_hm = IOUtils.loadKmers(new File[]{file}, 0,
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            if (hm == null) {
                hm = filt_hm;
                hm.resetValues();
            }

            Iterator<MutableLongShortEntry> it = filt_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                short value = entry.getValue();

                if (value > 0) {
                    hm.put(key, (short) (hm.getWithZero(key) + 1));
                }
            }
        }
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

        info("Selecting  k-mers...");
        int thresh = inputFiles.get().length * percentPresent.get()/100;
        BigLong2ShortHashMap hm_filtered = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 12);

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            short val = entry.getValue();

            if (val >= thresh) {
                hm_filtered.put(key, val);
            }
        }
        hm.reset();
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        info("Calculating presence of selected k-mers...");
        File outFile = new File(outDir, "selected_kmers_"+ percentPresent.get() +".txt");
        PrintWriter out;
        try {
            out = new PrintWriter(outFile);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Couldn't open output file", e);
        }
        debug("Starting to print k-mers to " + outFile.getPath());

        {
            Iterator<MutableLongShortEntry> it_filtered = hm_filtered.entryIterator();
            while (it_filtered.hasNext()) {
                MutableLongShortEntry entry = it_filtered.next();
                long key = entry.getKey();
                out.print("\t" + new ShortKmer(key, k.get()));
            }
            out.println();
        }


        for (File file: inputFiles.get()) {
            BigLong2ShortHashMap file_hm = IOUtils.loadKmers(new File[]{file}, 0,
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            out.print(file.getName().replace(".kmers.bin", ""));

            Iterator<MutableLongShortEntry> it_filtered = hm_filtered.entryIterator();
            while (it_filtered.hasNext()) {
                MutableLongShortEntry entry = it_filtered.next();
                long key = entry.getKey();
                out.print("\t" + file_hm.getWithZero(key));
            }
            out.println();
        }
        out.close();

        info("K-mers printed to " + outFile.getPath());
    }

    @Override
    protected void cleanImpl() {
    }

    public KmersPerSampleCounter() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new KmersPerSampleCounter().mainImpl(args);
    }
}