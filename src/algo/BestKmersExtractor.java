package algo;

import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.Iterator;

/**
 * Created by -- on 13.09.2022.
 */

public class BestKmersExtractor {

    final BigLong2ShortHashMap hm;
    final int[] ranks;

    public BestKmersExtractor(BigLong2ShortHashMap hm, DataInputStream streamIn) {
        this.hm = hm;
        ranks = new int[(int) Math.min(hm.capacity(), Integer.MAX_VALUE)];
        try {
            for (int i = 0; i < Integer.MAX_VALUE; i++) {
                ranks[i] = streamIn.readInt();
            }
        } catch (IOException ex) {
            ex.getMessage();
        }
    }

    public void outTopNKmers(int n, DataOutputStream outStream) throws IOException {
        Iterator<MutableLongShortEntry> it = this.hm.entryIterator();
        int i = 0;
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            int value = ranks[i];
            if (value < n) {
                long key = entry.getKey();
                short coverage = this.hm.get(key);
                outStream.writeLong(key);
                outStream.writeShort(coverage);
            }
            i++;
        }
    }
}

