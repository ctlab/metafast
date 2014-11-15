package algo;

import it.unimi.dsi.fastutil.longs.Long2IntMap;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.kmers.ShortKmer;
import ru.ifmo.genetics.structures.map.*;

import java.util.HashMap;
import java.util.Iterator;

public class HashMapOperations {

    public static byte getLeftNucleotide(BigLong2ShortHashMap hm, ShortKmer kmer, int freqThreshold) {
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

    public static byte getRightNucleotide(BigLong2ShortHashMap hm, ShortKmer kmer, int freqThreshold) {
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

    public static void banBranchingKmers(BigLong2ShortHashMap hm,
                                         int freqThreshold,
                                         int k,
                                         Logger logger) {
        final short BAN_VALUE = -1;
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
            if (getLeftNucleotide(hm, kmer, freqThreshold) == -2 ||
                    getRightNucleotide(hm, kmer, freqThreshold) == -2) {
                hm.put(key, BAN_VALUE);
                totalBanned += value;
                uniqueBanned++;
            }
        }

        logger.info("Total k-mers = " + totalKmers + ", unique k-mers = " + uniqueKmers);
        logger.info("Total k-mers under threshold = " + totalUnderThreshold + ", unique = " + uniqueUnderThreshold);
        logger.info("Total k-mers banned = " + totalBanned + ", unique = " + uniqueBanned);
    }
}
