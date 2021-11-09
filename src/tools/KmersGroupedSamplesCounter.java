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
 * Created by -- on 08.11.2021.
 */
public class KmersGroupedSamplesCounter  extends Tool {

    public static final String NAME = "kmers-grouped-counter";

    public static final String DESCRIPTION = "Count number of samples from 3 groups containing specified k-mers";


    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size (maximum 31 due to realization details)")
            .create());

    public final Parameter<File[]> kmersFile = addParameter(new FileMVParameterBuilder("kmers-file")
            .important()
            .withShortOpt("kf")
            .withDescription("binary file with kmers")
            .create());

    public final Parameter<File[]> CDFiles = addParameter(new FileMVParameterBuilder("cd-kmers")
            .mandatory()
            .withShortOpt("cd")
            .withDescription("list of input files with k-mers in binary format for CD group")
            .create());

    public final Parameter<File[]> UCFiles = addParameter(new FileMVParameterBuilder("uc-kmers")
            .mandatory()
            .withShortOpt("uc")
            .withDescription("list of input files with k-mers in binary format for UC group")
            .create());

    public final Parameter<File[]> nonIBDFiles = addParameter(new FileMVParameterBuilder("nonibd-kmers")
            .mandatory()
            .withShortOpt("nonibd")
            .withDescription("list of input files with k-mers in binary format for nonIBD group")
            .create());


    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous")
            .withDefaultValue(1)
            .create());

    public final Parameter<File> outputDir = addParameter(new FileParameterBuilder("output-dir")
            .withDescription("Output directory")
            .withDefaultValue(workDir.append("kmers"))
            .create());

    public final Parameter<File> statsDir = addParameter(new FileParameterBuilder("stats-dir")
            .withDescription("Directory with statistics")
            .withDefaultValue(workDir.append("stats"))
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
        File stDir = statsDir.get();
        if (!stDir.exists()) {
            stDir.mkdirs();
        }

        BigLong2ShortHashMap cd_hm = IOUtils.loadKmers(kmersFile.get(), 0,
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        BigLong2ShortHashMap uc_hm = IOUtils.loadKmers(kmersFile.get(), 0,
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);
        BigLong2ShortHashMap nonibd_hm = IOUtils.loadKmers(kmersFile.get(), 0,
                availableProcessors.get(), logger);
        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

        cd_hm.resetValues();
        uc_hm.resetValues();
        nonibd_hm.resetValues();

        debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);


        for (File file : CDFiles.get()) {
            BigLong2ShortHashMap filt_hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            Iterator<MutableLongShortEntry> it = cd_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();

                if (filt_hm.get(key) > maximalBadFrequency.get()) {
                    cd_hm.put(key, (short) (cd_hm.getWithZero(key) + 1));
                }
            }
        }

        for (File file : UCFiles.get()) {
            BigLong2ShortHashMap filt_hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            Iterator<MutableLongShortEntry> it = uc_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();

                if (filt_hm.get(key) > maximalBadFrequency.get()) {
                    uc_hm.put(key, (short) (uc_hm.getWithZero(key) + 1));
                }
            }
        }

        for (File file : nonIBDFiles.get()) {
            BigLong2ShortHashMap filt_hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            Iterator<MutableLongShortEntry> it = nonibd_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();

                if (filt_hm.get(key) > maximalBadFrequency.get()) {
                    nonibd_hm.put(key, (short) (nonibd_hm.getWithZero(key) + 1));
                }
            }
        }


        File outFile = new File(outDir, "kmers.groups.txt");
        PrintWriter out;
        try {
            out = new PrintWriter(outFile);
        } catch (FileNotFoundException e) {
            throw new ExecutionFailedException("Couldn't open output file", e);
        }

        debug("Starting to print k-mers to " + outFile.getPath());

        out.println("Kmer\tcd_count\tuc_count\tnonibd_count");
        Iterator<MutableLongShortEntry> it = cd_hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            out.println(new ShortKmer(key, k.get()) + "\t" + entry.getValue() +
                    "\t" + uc_hm.getWithZero(key) + "\t" + nonibd_hm.getWithZero(key));
        }
        out.close();
        info("K-mers printed to " + outFile.getPath());

    }

    @Override
    protected void cleanImpl() {
    }

    public KmersGroupedSamplesCounter() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new KmersGroupedSamplesCounter().mainImpl(args);
    }
}
