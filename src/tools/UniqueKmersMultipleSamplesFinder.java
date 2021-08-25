package tools;

import io.IOUtils;
import ru.ifmo.genetics.statistics.Timer;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;
import ru.ifmo.genetics.utils.tool.values.InMemoryValue;
import ru.ifmo.genetics.utils.tool.values.InValue;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Created by -- on 16.11.2020.
 */
public class UniqueKmersMultipleSamplesFinder extends Tool {

    public static final String NAME = "unique-kmers-multi";

    public static final String DESCRIPTION = "Output k-mers present in one dataset in fixed number of samples and missing in other";


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

    public final Parameter<File[]> filterFiles = addParameter(new FileMVParameterBuilder("filter-kmers")
            .mandatory()
            .withDescription("list of input files with k-mers in binary format used for filtering")
            .create());

    public final Parameter<Integer> minSamples = addParameter(new IntParameterBuilder("min-samples")
            .optional()
            .withDescription("minimal number of samples k-mer to be present in")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> maxSamples = addParameter(new IntParameterBuilder("max-samples")
            .optional()
            .withDescription("maximal number of samples k-mer to be present in")
            .withDefaultValue(1)
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a k-mer to be assumed erroneous (while saving kmers from sequences)")
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


    private final InMemoryValue<File[]> resultingKmerFilesPr = new InMemoryValue<File[]>();
    public final InValue<File[]> resultingKmerFiles =
            addOutput("resulting-kmers-file", resultingKmerFilesPr, File[].class);


    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        if (k.get() <= 0) {
            error("The size of k-mer must be at least 1.");
            System.exit(1);
        }
        if (k.get() > 31) {
            error("The size of k-mer must be no more than 31.");
            System.exit(1);
        }


        Timer t = new Timer();

        BigLong2ShortHashMap hm = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 8);
        BigLong2ShortHashMap hm_cnt = new BigLong2ShortHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 8);
        for (File file : inputFiles.get()) {
            BigLong2ShortHashMap tmp_hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            Iterator<MutableLongShortEntry> it = tmp_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                short value = entry.getValue();

                if (value > maximalBadFrequency.get()) {
                    hm.put(key, (short) (hm.getWithZero(key) + value));
                    hm_cnt.put(key, (short) (hm_cnt.getWithZero(key) + 1));
                }
            }
        }


        for (File file : filterFiles.get()) {
            BigLong2ShortHashMap filt_hm = IOUtils.loadKmers(new File[]{file}, maximalBadFrequency.get(),
                    availableProcessors.get(), logger);
            debug("Memory used = " + Misc.usedMemoryAsString() + ", time = " + t);

            Iterator<MutableLongShortEntry> it = filt_hm.entryIterator();
            while (it.hasNext()) {
                MutableLongShortEntry entry = it.next();
                long key = entry.getKey();
                short value = entry.getValue();

                if (value > maximalBadFrequency.get() && hm.get(key) > maximalBadFrequency.get()) {
                    hm.put(key, (short) 0);
                }
            }
        }



        File outDir = outputDir.get();
        if (!outDir.exists()) {
            outDir.mkdirs();
        }
        File stDir = statsDir.get();
        if (!stDir.exists()) {
            stDir.mkdirs();
        }

        for (int i = minSamples.get(); i < maxSamples.get() + 1; i++) {
            File outFile = new File(outDir, "filtered_" + i + ".kmers.bin");


            debug("Starting to print k-mers to " + outFile.getPath());
            long c = 0;
            try {
                c = IOUtils.filterAndPrintKmers(hm, hm_cnt, maximalBadFrequency.get(), i-1, outFile);
            } catch (IOException e) {
                e.printStackTrace();
            }
            info(NumUtils.groupDigits(hm.size()) + " k-mers found, "
                    + NumUtils.groupDigits(c) + " (" + String.format("%.1f", c * 100.0 / hm.size())
                    + "%) of them is good (present in one dataset and missing in other)");

            if (hm.size() == 0) {
                warn("No k-mers found in reads! Perhaps you reads file is empty or k-mer size is too big");
            } else if (c == 0 || c < (long) (hm.size() * 0.03)) {
                warn("Too few good k-mers were found! Perhaps you should decrease k-mer size or --maximal-bad-frequency value");
            }
            long allKmersNumber = (1L << (2 * k.get())) / 2;  // (4^k)/2
            if (hm.size() == allKmersNumber) {
                warn("All possible k-mers were found in reads! Perhaps you should increase k-mer size");
            } else if (hm.size() >= (long) (allKmersNumber * 0.99)) {
                warn("Almost all possible k-mers were found in reads! Perhaps you should increase k-mer size");
            }

            info("Good k-mers printed to " + outFile.getPath());
        }
    }

    @Override
    protected void cleanImpl() {
        resultingKmerFilesPr.set(new File[]{new File(outputDir.get(), "filtered_" + minSamples.get() + ".kmers.bin")});
    }

    public UniqueKmersMultipleSamplesFinder() {
        super(NAME, DESCRIPTION);
    }

    public static void main(String[] args) {
        new UniqueKmersMultipleSamplesFinder().mainImpl(args);
    }


}
