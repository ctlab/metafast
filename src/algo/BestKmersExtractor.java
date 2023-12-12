package algo;

import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.structures.map.MutableLongShortEntry;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.Iterator;


public class BestKmersExtractor {

    final BigLong2ShortHashMap hm;
    final int[] ranks;

    public BestKmersExtractor(BigLong2ShortHashMap hm, DataInputStream ranksStream) throws IOException {
        this.hm = hm;
        ranks = new int[(int) Math.min(hm.capacity(), Integer.MAX_VALUE)];
        for (int i = 0; i < Integer.MAX_VALUE; i++) {
            ranks[i] = ranksStream.readInt();
        }
    }

    public void outTopNKmers(int n, DataOutputStream outStream) throws IOException {
        Iterator<MutableLongShortEntry> it = hm.entryIterator();
        int i = 0;
        while (it.hasNext()) {
            MutableLongShortEntry entry = it.next();
            int value = ranks[i];
            if (value < n) {
                long key = entry.getKey();
                outStream.writeLong(key);
                outStream.writeShort(hm.get(key));
            }
            i++;
        }
    }
}

