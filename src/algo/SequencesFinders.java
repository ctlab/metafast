package algo;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.executors.BlockingThreadPoolExecutor;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import structures.Sequence;

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

public class SequencesFinders {

    public static List<Sequence> thresholdStrategy(ArrayLong2IntHashMap hm,
                                                   int availableProcessors,
                                                   int freqThreshold,
                                                   int lenThreshold,
                                                   int k) throws InterruptedException {
        //List<Sequence> ans = Collections.synchronizedList(new ArrayList<Sequence>());
        List<Sequence> ans = new CopyOnWriteArrayList<Sequence>();

        BlockingThreadPoolExecutor executor = new BlockingThreadPoolExecutor(availableProcessors);

        for (int i = 0; i < hm.hm.length; ++i) {
            executor.blockingExecute(new
                    AddSequencesShiftingRightTask(hm, hm.hm[i], freqThreshold, ans, lenThreshold, k));
        }

        System.out.println(executor.getTaskCount());
        executor.shutdownAndAwaitTermination();
        return ans;
    }

    public static Sequence getSequenceShiftingRight(ArrayLong2IntHashMap hm,
                                                    ShortKmer kmer,
                                                    int freqThreshold) {
        int value = hm.get(kmer.toLong());

        StringBuilder sequenceSB = new StringBuilder(kmer.toString());
        long seqWeight = value, minWeight = value, maxWeight = value;

        byte rightNuc = HashMapOperations.getRightNucleotide(hm, kmer, freqThreshold);
        if (rightNuc < 0) {
            return null;
        }

        while (true) {
            kmer.shiftRight(rightNuc);
            byte nextRightNuc = HashMapOperations.getRightNucleotide(hm, kmer, freqThreshold);

            if (nextRightNuc < 0 || HashMapOperations.getLeftNucleotide(hm, kmer, freqThreshold) < 0) {
                break;
            }

            sequenceSB.append(DnaTools.toChar(rightNuc));

            value = hm.get(kmer.toLong());
            seqWeight += value;
            minWeight = Math.min(minWeight, value);
            maxWeight = Math.max(maxWeight, value);

            rightNuc = nextRightNuc;
        }

        return new Sequence(sequenceSB.toString(), seqWeight, minWeight, maxWeight);
    }

}
