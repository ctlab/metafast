package structures.map;

import it.unimi.dsi.fastutil.HashCommon;
import org.apache.commons.lang.mutable.MutableLong;
import org.apache.log4j.Logger;
import ru.ifmo.genetics.structures.set.BigLongHashSet;
import ru.ifmo.genetics.structures.set.LongHashSet;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.Tool;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Memory-efficient map based on many hash tables with open addressing.<br></br>
 * It is resizable (not full support!). Map is synchronized.<br></br>
 * It can contain up to 2^60 (~10^18) elements.<br></br>
 * <br></br>
 */
public class BigLong2BitShortaHashMap implements Long2BitShortaHashMapInterface {
    private static final Logger logger = Logger.getLogger("BigLong2BitLongaHashMap");

    public Long2BitShortaHashMap[] maps;
    protected int mask;


    public BigLong2BitShortaHashMap(long capacity) {
        this(
                NumUtils.getPowerOf2(capacity >> 20 +
                        ((capacity & ((1 << 20) - 1)) == 0 ? 0 : 1)
                ),
                20    // 1 M elements per small map
        );
    }

    public BigLong2BitShortaHashMap(int logSmallMapNumber, int logSmallCapacity) {
        this(logSmallMapNumber, logSmallCapacity, false);
    }

    public BigLong2BitShortaHashMap(int logSmallMapNumber, int logSmallCapacity, boolean debugInfo) {
        this(logSmallMapNumber, logSmallCapacity, false, 0);
    }

    public BigLong2BitShortaHashMap(int logSmallMapNumber, int logSmallCapacity, boolean debugInfo, int sizeBitSet) {
        if (logSmallMapNumber > 30) {
            throw new IllegalArgumentException("logSmallMapNumber > 30!");
        }

        int smallMapNumber = 1 << logSmallMapNumber;
        mask = smallMapNumber - 1;

        maps = new Long2BitShortaHashMap[smallMapNumber];
        for (int i = 0; i < smallMapNumber; i++) {
            maps[i] = new Long2BitShortaHashMap(logSmallCapacity, LongHashSet.DEFAULT_MAX_LOAD_FACTOR, sizeBitSet);
        }
        if (debugInfo) {
            Tool.debug(logger, "Created " + NumUtils.groupDigits(smallMapNumber) + " small Long2BitSetHashMaps");
        }
    }


    @Override
    public short[] set(long key, int bitIndex) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].set(key, bitIndex);
    }

    @Override
    public boolean get(long key, int bitIndex) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].get(key, bitIndex);
    }

    @Override
    public short[] get(long key) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].get(key);
    }

    @Override
    public short[] getWithEmpty(long key) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].getWithEmpty(key);
    }

    public int getCardinality(long key) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].getCardinality(key);
    }

    public int getCardinality(long key, int from, int to) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].getCardinality(key, from, to);
    }

    @Override
    public boolean contains(long key) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        return maps[n].contains(key);
    }

    @Override
    public long size() {
        long size = 0;
        for (Long2BitShortaHashMap map : maps) {
            size += map.size();
        }
        return size;
    }

    @Override
    public long capacity() {
        long capacity = 0;
        for (Long2BitShortaHashMap map : maps) {
            capacity += map.capacity();
        }
        return capacity;
    }



    // --------------  Other methods from interface Long2BitLongaHashMapInterface  ---------------

    @Override
    public void reset() {
        for (Long2BitShortaHashMap map : maps) {
            map.reset();
        }
    }
    @Override
    public void resetValues() {
        for (Long2BitShortaHashMap map : maps) {
            map.resetValues();
        }
    }


    long[] off;

    @Override
    public void prepare() {
        off = new long[maps.length];
        off[0] = 0;
        for (int i = 1; i < maps.length; i++) {
            off[i] = off[i - 1] + maps[i - 1].maxPosition() + 1;
        }
    }

    @Override
    public long maxPosition() {
        return maps.length == 0 ? -1 :
                off[maps.length - 1] + maps[maps.length - 1].maxPosition();
    }

    /**
     * Working ONLY if small sets stay unchanged!
     * Call prepare() before using.
     */
    @Override
    public long getPosition(long key) {
        int n = HashCommon.murmurHash3((int) key) & mask;
        long pos = maps[n].getPosition(key);
        return pos == -1 ? -1 : (off[n] + pos);
    }

    @Override
    public long keyAt(long pos) {
        int n = Arrays.binarySearch(off, pos);
        if (n < 0) {
            n = (-n - 1) - 1;
        }
        return maps[n].keyAt(pos - off[n]);
    }
    @Override
    public short[] valueAt(long pos) {
        int n = Arrays.binarySearch(off, pos);
        if (n < 0) {
            n = (-n - 1) - 1;
        }
        return maps[n].valueAt(pos - off[n]);
    }
    @Override
    public boolean containsAt(long pos) {
        int n = Arrays.binarySearch(off, pos);
        if (n < 0) {
            n = (-n - 1) - 1;
        }
        return maps[n].containsAt(pos - off[n]);
    }


    @Override
    public void write(DataOutput out) throws IOException {
        out.writeInt(maps.length);

        for (Long2BitShortaHashMap map : maps) {
            map.write(out);
        }
    }

    @Override
    public void readFields(DataInput in) throws IOException {
        int len = in.readInt();
        if (Integer.bitCount(len) != 1) {
            throw new RuntimeException("Length is not a power of two!");
        }
        maps = new Long2BitShortaHashMap[len];
        mask = maps.length - 1;

        for (int i = 0; i < len; i++) {
            maps[i] = new Long2BitShortaHashMap(2);
            maps[i].readFields(in);
        }
    }


    @Override
    public Iterator<MutableLong> iterator() {
        return new BigLongHashSet.MyIterator(maps);
    }

    @Override
    public Iterator<MutableLongBitShortaEntry> entryIterator() {
        return new MyIterator();
    }

    class MyIterator implements Iterator<MutableLongBitShortaEntry> {
        private int index;
        private Iterator<MutableLongBitShortaEntry> it = null;

        public MyIterator() {
            index = 0;
            if (index < maps.length) {
                it = maps[index].entryIterator();
            }
        }

        @Override
        public boolean hasNext() {
            while (index < maps.length) {
                if (it.hasNext()) {
                    return true;
                }
                index++;
                if (index < maps.length) {
                    it = maps[index].entryIterator();
                }
            }
            return false;
        }

        @Override
        public MutableLongBitShortaEntry next() {
            if (hasNext()){
                return it.next();
            }
            throw new NoSuchElementException();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
