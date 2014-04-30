package algo;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.DnaTools;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import structures.Sequence;

import java.util.ArrayList;
import java.util.List;

public class SequencesFinders {
    /*
    TODO: Parallel it
     */
    public static List<Sequence> thresholdStrategy(ArrayLong2IntHashMap hm,
                                                   int freqThreshold,
                                                   int lenThreshold,
                                                   int k,
                                                   Logger logger) {
        int LOG_GAP = 10000;

        List<Sequence> ans = new ArrayList<Sequence>();

        for (int i = 0; i < hm.hm.length; ++i) {
            for (Long2IntMap.Entry entry : hm.hm[i].long2IntEntrySet()) {
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

                Sequence sequence = getSequenceShiftingRight(hm, kmer, freqThreshold);

                if (sequence.length() >= lenThreshold) {
                    ans.add(sequence);

                    if (ans.size() % LOG_GAP == 0) {
                        logger.debug("sequenceId = " + ans.size() + ", last len = " + sequence.length());
                    }
                }

            }
        }

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
