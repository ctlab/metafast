package algo;

import org.apache.commons.lang.ArrayUtils;
import ru.ifmo.genetics.dna.kmers.ShortKmer;

import java.util.Arrays;

public class KmerOperations {
    public static long[] possibleNeighbours(long kmerRepr, int k) {
        long[] ans = new long[8];

        ShortKmer goRight = new ShortKmer(kmerRepr, k);
        goRight.shiftRight((byte) 0);
        ans[0] = goRight.toLong();
        ShortKmer goLeft = new ShortKmer(kmerRepr, k);
        goLeft.shiftLeft((byte) 0);
        ans[1] = goLeft.toLong();

        for (byte nuc = 1; nuc <= 3; nuc++) {
            goRight.updateAt(k - 1, nuc);
            ans[nuc*2] = goRight.toLong();
            goLeft.updateAt(0, nuc);
            ans[nuc*2 + 1] = goLeft.toLong();
        }
        return ans;
    }

    public static long[] calcNeighbours(long kmer, int k) {
        return ArrayUtils.addAll(leftNeighbours(kmer, k), rightNeighbours(kmer, k));
    }

    public static long[] rightNeighbours(long kmer, int k) {
        long mask = (1L << (2 * k)) - 1;
        long[] ans = new long[] {(kmer << 2) & mask,
                ((kmer << 2) & mask) | 1,
                ((kmer << 2) & mask) | 2,
                ((kmer << 2) & mask) | 3};
        for (int i = 0; i < ans.length; i++) {
            long rc = rc(ans[i], k);
            if (rc < ans[i]) {
                ans[i] = rc;
            }
        }
        return ans;
    }

    public static long[] leftNeighbours(long kmerRepr, int k) {
        long[] ans = new long[4];

        ShortKmer kmer = new ShortKmer(kmerRepr, k);
        byte rightNuc = kmer.nucAt(k - 1);

        for (byte nuc = 0; nuc <= 3; nuc++) {
            kmer.shiftLeft(nuc);
            ans[nuc] = kmer.toLong();
            kmer.shiftRight(rightNuc);
            assert kmer.toLong() == kmerRepr;
        }
        return ans;
    }

    public static long rc(long kmer, long k) {
        kmer = ((kmer & 0x3333333333333333L) << 2) | ((kmer & 0xccccccccccccccccL) >>> 2);
        kmer = ((kmer & 0x0f0f0f0f0f0f0f0fL) << 4) | ((kmer & 0xf0f0f0f0f0f0f0f0L) >>> 4);
        kmer = ((kmer & 0x00ff00ff00ff00ffL) << 8) | ((kmer & 0xff00ff00ff00ff00L) >>> 8);
        kmer = ((kmer & 0x0000ffff0000ffffL) << 16) | ((kmer & 0xffff0000ffff0000L) >>> 16);
        kmer = ((kmer & 0x00000000ffffffffL) << 32) | ((kmer & 0xffffffff00000000L) >>> 32);

        kmer = ~kmer;

        return kmer >>> (64 - 2 * k);
    }

}
