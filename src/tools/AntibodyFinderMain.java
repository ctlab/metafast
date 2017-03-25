package tools;

import algo.KmerOperations;
import io.IOUtils;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileMVParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.IntParameterBuilder;

import java.io.*;
import java.util.ArrayDeque;
import java.util.Queue;

public class AntibodyFinderMain extends Tool {
    public static final String NAME = "antibody-sequences-finder";
    public static final String DESCRIPTION = "Antibody sequences finder in De Bruijn graph";

    static final int LOAD_TASK_SIZE = 1 << 15;
    static final long MAX_SIZE = 10000000000L;

    static final int FRAGMENT_KMER_HM_VALUE = 30;

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<Integer> shift = addParameter(new IntParameterBuilder("shift")
            .withDescription("shift from the start")
            .withDefaultValue(50)
            .create());

    public final Parameter<Integer> maxDistance = addParameter(new IntParameterBuilder("max-distance")
            .mandatory()
            .withShortOpt("d")
            .withDescription("distance from constant fragment")
            .create());

    public final Parameter<File> constantFragmentFile = addParameter(new FileParameterBuilder("fragment-file")
            .mandatory()
            .withShortOpt("ff")
            .withDescription("file with constant fragment in FASTA")
            .create());

    public final Parameter<File[]> readsFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files in BINQ")
            .create());

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .mandatory()
            .withShortOpt("b")
            .withDescription("maximal frequency for a kmer to be assumed erroneous")
            .create());


    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {
        int LEN = k.get();

        String seq = loadConstantFragment();
        info("Constant fragment length = " + seq.length());

        BigLong2ShortHashMap readsHM =
                IOUtils.loadReads(readsFiles.get(), LEN, 0, availableProcessors.get(), logger);

        ShortKmer kmer = new ShortKmer(seq.substring(0, LEN));
        for (int pos = LEN; pos < seq.length(); pos++) {
            kmer.shiftRight(DnaTools.fromChar(seq.charAt(pos)));

            String deb = kmer.toString() + " " + readsHM.get(kmer.toLong());
            for (long neighbour : KmerOperations.leftNeighbours(kmer.fwKmer(), LEN)) {
//                deb += " " + new ShortKmer(neighbour, k.get()).toString() + " " + readsHM.get(neighbour);
                deb += " " + readsHM.get(neighbour);
            }
//            for (long neighbour : algo.KmerOperations.leftNeighbours(kmer.toLong(), k.get())) {
//                deb += " " + readsHM.get(algo.KmerOperations.rc(neighbour, k.get()));
//            }


            readsHM.addAndBound(kmer.toLong(), (short) (maximalBadFrequency.get() + 1));
//            deb += " " + kmer.toLong() + " " + readsHM.get(kmer.toLong());
            debug(deb);
        }

//        ShortKmer startKmer = new ShortKmer(seq.substring(seq.length() - k.get() - shift.get(),
//                seq.length() - shift.get()));
        ShortKmer startKmer = new ShortKmer(seq.substring(shift.get(), shift.get() + LEN));
        info("Start k-mer = " + startKmer.toString() + " " + readsHM.get(startKmer.toLong()));

        int depth = maxDistance.get() + shift.get();
        info("Depth from start k-mer = " + depth);

        ArrayLong2IntHashMap distFromStart =
                new ArrayLong2IntHashMap((int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4);
        distFromStart.add(startKmer.fwKmer(), 1);

        int freqThreshold = maximalBadFrequency.get();
        int[] uniqueKmers = new int[depth + 1];
        long[] totalKmers = new long[depth + 1];

        Queue<Long> queue = new ArrayDeque<Long>();
        queue.add(startKmer.fwKmer());

        PrintWriter kmersPW = null;
        try {
            kmersPW = new PrintWriter(workDir + File.separator + "kmers");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        int prevDepth = 1;
        while (queue.size() > 0) {
            ShortKmer currentKmer = new ShortKmer(queue.poll(), LEN);
            int kmerDepth = distFromStart.get(currentKmer.fwKmer());
            if (kmerDepth > depth) {
                break;
            }

            if (kmerDepth > prevDepth) {
                prevDepth++;
                kmersPW.println();
            }
            kmersPW.print(currentKmer.toString() + " ");

            uniqueKmers[kmerDepth]++;
            totalKmers[kmerDepth] += readsHM.getWithZero(currentKmer.toLong());

            byte rightNuc = currentKmer.nucAt(LEN - 1);
            for (byte nuc = 0; nuc <= 3; nuc++) {
                currentKmer.shiftLeft(nuc);
                long fwKmer = currentKmer.fwKmer(), kmerRepr = currentKmer.toLong();
                currentKmer.shiftRight(rightNuc);

                if (distFromStart.get(fwKmer) == 0 && readsHM.get(kmerRepr) > freqThreshold) {
                    distFromStart.add(fwKmer, kmerDepth + 1);
                    queue.add(fwKmer);
                }
            }
        }

        kmersPW.close();

        for (int i = 0; i < uniqueKmers.length; i++) {
            debug(i + " " + uniqueKmers[i] + " " + totalKmers[i]);
        }

        try {
            PrintWriter statPW = new PrintWriter(workDir + File.separator + "stat-b" + freqThreshold);
            for (int i = 0; i < uniqueKmers.length; i++) {
                statPW.println(i + " " + uniqueKmers[i] + " " + totalKmers[i]);
            }
            statPW.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }

    private String loadConstantFragment() throws ExecutionFailedException {
        String seq = "";
        try {
            BufferedReader reader = new BufferedReader(new FileReader(constantFragmentFile.get()));
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.charAt(0) != '>') {
                    seq += line.trim();
                }
            }
            reader.close();
        } catch (IOException e) {
            throw new ExecutionFailedException(
                    "Couldn't load constant fragment from " + constantFragmentFile.get(), e);
        }
        return seq;
    }

    @Override
    protected void cleanImpl() {
    }

    public static void main(String[] args) {
        new AntibodyFinderMain().mainImpl(args);
    }

    public AntibodyFinderMain() {
        super(NAME, DESCRIPTION);
    }
}
