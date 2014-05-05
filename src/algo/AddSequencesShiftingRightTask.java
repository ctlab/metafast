package algo;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import structures.Sequence;

import java.util.List;

/**
 * Task class to add sequences from single <tt>Long2IntOpenHashMap</tt>
 * taken from <tt>ArrayLong2IntHashMap</tt>
 * @author Vladimir Ulyantsev
 */
public class AddSequencesShiftingRightTask implements Runnable {

    ArrayLong2IntHashMap hm;
    Long2IntOpenHashMap openHM;
    int freqThreshold;
    List<Sequence> sequenceList;
    int lenThreshold;
    int k;

    public AddSequencesShiftingRightTask(ArrayLong2IntHashMap hm,
                                         Long2IntOpenHashMap openHM,
                                         int freqThreshold,
                                         List<Sequence> sequenceList,
                                         int lenThreshold,
                                         int k) {
        this.hm = hm;
        this.openHM = openHM;
        this.freqThreshold = freqThreshold;
        this.sequenceList = sequenceList;
        this.lenThreshold = lenThreshold;
        this.k = k;
    }

    @Override
    public void run() {
        for (Long2IntMap.Entry entry : openHM.long2IntEntrySet()) {
            int value = entry.getIntValue();
            if (value <= freqThreshold) {
                continue;
            }

            long key = entry.getLongKey();
            ShortKmer kmer = new ShortKmer(key, k);

            if (HashMapOperations.getRightNucleotide(hm, kmer, freqThreshold) < 0 ||
                    HashMapOperations.getLeftNucleotide(hm, kmer, freqThreshold) >= 0) {
                continue;
            }

            processSequence(kmer);
        }
    }

    private void processSequence(ShortKmer startKmer) {
        int value = hm.get(startKmer.toLong());

        StringBuilder sequenceSB = new StringBuilder(startKmer.toString());
        long seqWeight = value, minWeight = value, maxWeight = value;

        byte rightNuc = HashMapOperations.getRightNucleotide(hm, startKmer, freqThreshold);
        if (rightNuc < 0) {
            return;
        }

        ShortKmer kmer = new ShortKmer(startKmer);
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

        if (sequenceSB.length() >= lenThreshold) {
            sequenceList.add(new Sequence(sequenceSB.toString(), seqWeight, minWeight, maxWeight));
        }
    }
}
