package tools;

import io.IOUtils;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.KmerIteratorFactory;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.dna.kmers.ShortKmerIteratorFactory;
import ru.ifmo.genetics.structures.map.*;
import ru.ifmo.genetics.utils.Misc;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.Tool;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

public class SupergraphSeqBuilderMain extends Tool {
    public static final String NAME = "supergraph-sequence-builder";
    public static final String DESCRIPTION = "NOT COMPLETED";


    private final int STAT_LEN = 1024;

    public final Parameter<Integer> maximalBadFrequency = addParameter(new IntParameterBuilder("maximal-bad-frequence")
            .optional()
            .withShortOpt("b")
            .withDescription("maximal frequency for a kmer to be assumed erroneous in single library")
            .create());

    public final Parameter<Integer> bottomCutPercent = addParameter(new IntParameterBuilder("bottom-cut-percent")
            .optional()
            .withShortOpt("bp")
            .withDescription("k-mers percent to be assumed erroneous")
            .create());

    public final Parameter<Integer> supergraphFreq = addParameter(new IntParameterBuilder("supergraph-frequence")
            .mandatory()
            .withShortOpt("sb")
            .withDescription("maximal frequency for a kmer to be assumed erroneous in supergraph")
            .create());

    public final Parameter<Integer> k = addParameter(new IntParameterBuilder("k")
            .mandatory()
            .withShortOpt("k")
            .withDescription("k-mer size")
            .create());

    public final Parameter<Integer> sequenceLen = addParameter(new IntParameterBuilder("sequence-len")
            .mandatory()
            .withShortOpt("l")
            .withDescription("sequence minimal length to be written")
            .create());

    public final Parameter<Long> maxSize = addParameter(new LongParameterBuilder("max-size")
            .optional()
            .withDescription("maximal hashset size")
            .withDefaultValue(NumUtils.highestBits(Misc.availableMemory() / 42, 3))
            .memoryParameter()
            .create());

    public final Parameter<File[]> inputFiles = addParameter(new FileMVParameterBuilder("reads")
            .withShortOpt("i")
            .mandatory()
            .withDescription("list of input files")
            .create());

    public final Parameter<KmerIteratorFactory> kmerIteratorFactory = Parameter.createParameter(
            new KmerIteratorFactoryParameterBuilder("kmer-iterator-factory")
                    .optional()
                    .withDescription("factory used for iterating through kmers")
                    .withDefaultValue(new ShortKmerIteratorFactory())
                    .create());

    private long MAX_SIZE;

    @Override
    protected void runImpl() throws ExecutionFailedException, IOException {

        MAX_SIZE = maxSize.get();

        if (maximalBadFrequency.get() != null && bottomCutPercent.get() != null) {
            throw new IllegalArgumentException("-b and -bp can not be set both");
        }

        debug("MAXIMAL_SIZE = " + MAX_SIZE);

        BigLong2ShortHashMap superHM = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors.get()) / Math.log(2)) + 4, 8);

        for (File file : inputFiles.get()) {
            addToSupergraph(superHM, file);
        }

        calcSequences(superHM, workDir + File.separator + "sequences.fasta");
    }

    private void addToSupergraph(BigLong2ShortHashMap superHM, File readsFile) throws ExecutionFailedException, IOException {
        BigLong2ShortHashMap hm =
                IOUtils.loadReads(new File[]{readsFile}, k.get(), 0, availableProcessors.get(), logger);

        int freqThreshold = getThreshold(hm);

        long uniqueKmers = 0, uniqueAdded = 0, newKmers = 0;

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            uniqueKmers++;
            int value = entry.getValue();
            if (value > freqThreshold) {
                uniqueAdded++;
                int old = superHM.addAndBound(entry.getKey(), (short)1);
                if (old == 0) {
                    newKmers++;
                }
            }
        }

        debug("Reads from " + readsFile.getPath() + " added with threshold " + freqThreshold);
        debug("Unique k-mers count = " + uniqueKmers + ", unique k-mers added " + uniqueAdded +
                ", new k-mers added = " + newKmers);
    }

    private int getThreshold(BigLong2ShortHashMap hm) {
        if (maximalBadFrequency.get() != null) {
            return maximalBadFrequency.get();
        }

        long totalKmers = 0;
        int[] stat = new int[STAT_LEN];
        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            int v = entry.getValue();
            totalKmers += v;
            if (v >= stat.length) {
                v = stat.length - 1;
            }
            ++stat[v];
        }

        if (bottomCutPercent.get() != null) {
            long kmersToCut = totalKmers * bottomCutPercent.get() / 100;
            long currentKmersCount = 0;
            for (int i = 0; i < stat.length - 1; i++) {
                if (currentKmersCount >= kmersToCut) {
                    return i;
                }
                currentKmersCount += (long) i * stat[i];
            }
        }

        if (maximalBadFrequency.get() == null) {
            int threshold = 1;
            long currentSum = 0;
            while (stat[threshold] * (long) threshold > stat[threshold + 1] * (long) (threshold + 1)) {
                currentSum += stat[threshold];
                if (currentSum * 2 > totalKmers) {
                    debug("Threshold search stopped at 50 %");
                    break;
                }
                threshold++;
            }
            return threshold;
        }

        debug("Something strange about threshold");
        return 0;
    }

    public void calcSequences(BigLong2ShortHashMap hm, String fastaFP) throws FileNotFoundException {
        int freqThreshold = supergraphFreq.get();
        int lenThreshold = sequenceLen.get();
        int kValue = k.get();

        banKmers(hm, freqThreshold, kValue);

        int sequenceId = 0;

        ArrayList<Integer> sequenceLen = new ArrayList<Integer>();
        ArrayList<Long> sequenceWeight = new ArrayList<Long>();

        long kmersInSeq = 0;
        long totalKmersInSequences = 0;

        PrintWriter fastaPW = new PrintWriter(fastaFP);

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            int value = entry.getValue();
            if (value <= freqThreshold) {
                continue;
            }
            ShortKmer kmer = new ShortKmer(key, kValue);

            if (getLeftNucleotide(kmer, hm, freqThreshold) >= 0) {
                continue;
            }

            StringBuilder sequenceSB = new StringBuilder(kmer.toString());
            long seqWeight = 0, minWeight = value, maxWeight = value;

            while (true) {
                long kmerRepr = kmer.toLong();
                value = hm.getWithZero(kmerRepr);
                seqWeight += value;
                minWeight = Math.min(minWeight, value);
                maxWeight = Math.max(maxWeight, value);

                hm.put(kmerRepr, BAN_VALUE);

                byte rightNuc = getRightNucleotide(kmer, hm, freqThreshold);
                if (rightNuc < 0) {
                    break;
                }
                sequenceSB.append(DnaTools.toChar(rightNuc));
                kmer.shiftRight(rightNuc);
            }

            if (sequenceSB.length() >= lenThreshold) {
                sequenceId++;
                String sequenceStr = sequenceSB.toString();

                sequenceLen.add(sequenceStr.length());
                sequenceWeight.add(seqWeight);

                totalKmersInSequences += seqWeight;
                kmersInSeq += sequenceStr.length() - kValue + 1;

                String seqInfo = String.format(">%d length=%d sum_weight=%d min_weight=%d max_weight=%d",
                        sequenceId, sequenceStr.length(), seqWeight, minWeight, maxWeight);
                fastaPW.println(seqInfo);
                fastaPW.println(sequenceStr);

                if (sequenceId % 10000 == 0) {
                    debug("sequenceId = " + sequenceId + ", last len = " + sequenceStr.length());
                }
            }

        }
        info(sequenceId + " sequences found");
        info(kmersInSeq + " unique k-mers out of " + hm.size() + " in sequences");
        info("Total k-mers in sequences = " + totalKmersInSequences);
        info("N50 value of sequences = " + getN50(sequenceLen));

        dumpSeqInfo(sequenceLen, sequenceWeight, workDir + File.separator + "seq-info");

        fastaPW.close();
    }

    final short BAN_VALUE = -1;
    private void banKmers(BigLong2ShortHashMap hm, int freqThreshold, int k) {
        long totalKmers = 0, uniqueKmers = 0,
                totalBanned = 0, uniqueBanned = 0,
                totalUnderThreshold = 0, uniqueUnderThreshold = 0;

        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            long key = entry.getKey();
            int value = entry.getValue();
            totalKmers += value;
            uniqueKmers++;
            if (value <= freqThreshold) {
                totalUnderThreshold += value;
                uniqueUnderThreshold++;
                continue;
            }

            ShortKmer kmer = new ShortKmer(key, k);
            if (getLeftNucleotide(kmer, hm, freqThreshold) == -2 ||
                    getRightNucleotide(kmer, hm, freqThreshold) == -2) {
                hm.put(key, BAN_VALUE);
                totalBanned += value;
                uniqueBanned++;
            }
        }

        info("Total k-mers = " + totalKmers + ", unique k-mers = " + uniqueKmers);
        info("Total k-mers [<=] threshold = " + totalUnderThreshold + ", unique = " + uniqueUnderThreshold);
        info("Total k-mers banned = " + totalBanned + ", unique = " + uniqueBanned);
    }

    private static byte getLeftNucleotide(ShortKmer kmer, BigLong2ShortHashMap hm, int freqThreshold) {
        byte rightNuc = kmer.nucAt(kmer.length() - 1);
        byte ansNuc = -1;
        for (byte nuc = 0; nuc <= 3; nuc++) {
            kmer.shiftLeft(nuc);
            if (hm.get(kmer.toLong()) > freqThreshold) {
                if (ansNuc > -1) {
                    return -2;
                }
                ansNuc = nuc;
            }
            kmer.shiftRight(rightNuc);
        }
        return ansNuc;
    }

    private static byte getRightNucleotide(ShortKmer kmer, BigLong2ShortHashMap hm, int freqThreshold) {
        byte leftNuc = kmer.nucAt(0);
        byte ansNuc = -1;
        for (byte nuc = 0; nuc <= 3; nuc++) {
            kmer.shiftRight(nuc);
            if (hm.get(kmer.toLong()) > freqThreshold) {
                if (ansNuc > -1) {
                    return -2;
                }
                ansNuc = nuc;
            }
            kmer.shiftLeft(leftNuc);
        }
        return ansNuc;
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
        new SupergraphSeqBuilderMain().mainImpl(args);
    }

    public SupergraphSeqBuilderMain() {
        super(NAME, DESCRIPTION);
    }
}
