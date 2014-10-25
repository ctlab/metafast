package algo;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.ArrayLong2IntHashMap;
import ru.ifmo.genetics.structures.map.BigLong2IntHashMap;
import ru.ifmo.genetics.structures.map.MutableLongIntEntry;

import java.util.HashMap;
import java.util.Iterator;

public class HashMapOperations {

    public static byte getLeftNucleotide(BigLong2IntHashMap hm, ShortKmer kmer, int freqThreshold) {
        byte rightNuc = kmer.nucAt(kmer.length() - 1);
        byte ansNuc = -1;
        for (byte nuc = 0; nuc <= 3; nuc++) {
            kmer.shiftLeft(nuc);
            long neighbourRepr = kmer.toLong();
            kmer.shiftRight(rightNuc);

            if (hm.get(neighbourRepr) > freqThreshold) {
                if (ansNuc > -1) {
                    return -2;
                }
                ansNuc = nuc;
            }
        }
        return ansNuc;
    }

    public static byte getRightNucleotide(BigLong2IntHashMap hm, ShortKmer kmer, int freqThreshold) {
        byte leftNuc = kmer.nucAt(0);
        byte ansNuc = -1;
        for (byte nuc = 0; nuc <= 3; nuc++) {
            kmer.shiftRight(nuc);
            long neighbourRepr = kmer.toLong();
            kmer.shiftLeft(leftNuc);

            if (hm.get(neighbourRepr) > freqThreshold) {
                if (ansNuc > -1) {
                    return -2;
                }
                ansNuc = nuc;
            }
        }
        return ansNuc;
    }

    public static HashMap<ShortKmer, Integer> getNeighbours(BigLong2IntHashMap hm,
                                                            ShortKmer kmer,
                                                            int depth) {

        return null;
    }

    public static void banBranchingKmers(BigLong2IntHashMap hm,
                                         int freqThreshold,
                                         int k,
                                         Logger logger) {
        int BAN_VALUE = 1000000000;
        long totalKmers = 0, uniqueKmers = 0,
                totalBanned = 0, uniqueBanned = 0,
                totalUnderThreshold = 0, uniqueUnderThreshold = 0;

        Iterator<MutableLongIntEntry> it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongIntEntry entry = it.next();
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
            if (getLeftNucleotide(hm, kmer, freqThreshold) == -2 ||
                    getRightNucleotide(hm, kmer, freqThreshold) == -2) {
                hm.addAndBound(key, BAN_VALUE);
                totalBanned += value;
                uniqueBanned++;
            }
        }

        logger.info("Total k-mers = " + totalKmers + ", unique k-mers = " + uniqueKmers);
        logger.info("Total k-mers [<=] threshold = " + totalUnderThreshold + ", unique = " + uniqueUnderThreshold);
        logger.info("Total k-mers banned = " + totalBanned + ", unique = " + uniqueBanned);

        it = hm.entryIterator();
        while (it.hasNext()) {
            MutableLongIntEntry entry = it.next();
            int value = entry.getValue();
            if (value >= BAN_VALUE) {
                long key = entry.getKey();
                hm.addAndBound(key, -(value + 1));
            }
        }
    }
}
